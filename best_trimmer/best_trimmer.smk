# Configuration
configfile: "config.yaml"

# Define global variables
INPUT_DIR = config["input_dir"]
OUTPUT_DIR = config["output_dir"]
ITERATIONS = config["iterations"]
SAMPLES = glob_wildcards(os.path.join(INPUT_DIR, "{sample}_R1.fastq.gz")).sample
LOG_DIR = "logs"

# Ensure output directories exist
onstart:
    shell("""
        mkdir -p {OUTPUT_DIR}/params
        mkdir -p {OUTPUT_DIR}/results
        mkdir -p {OUTPUT_DIR}/combined
        mkdir -p {OUTPUT_DIR}/plots
        mkdir -p {LOG_DIR}
        mkdir -p {INPUT_DIR}
    """)

rule all:
    input:
        os.path.join(OUTPUT_DIR, "report/best_trimmer_report.html"),
        os.path.join(OUTPUT_DIR, "optimal_parameters.json"),
        os.path.join(OUTPUT_DIR, "best_trimmer_summary.tsv"),
        expand("{output_dir}/optimized_trimming/{sample}_R{read}.fastq.gz", 
               output_dir=OUTPUT_DIR, sample=SAMPLES, read=[1,2])

# Add checkpoint to handle dynamic sample discovery
checkpoint find_samples:
    output:
        directory(os.path.join(INPUT_DIR))
    shell:
        """
        if [ -d "{output}" ]; then
            # Check if directory has fastq files
            if [ -z "$(ls {output}/*_R[12].fastq.gz 2>/dev/null)" ]; then
                echo "Creating test FASTQ files in {output}"
                mkdir -p {output}
                # Create more realistic test data with randomized sequences
                python3 - <<'EOF'
import gzip
import random
import os

bases = ['A', 'C', 'G', 'T']
quals = ''.join([chr(x) for x in range(33, 74)])  # Phred+33 quality scores

def make_read(read_num):
    seq = ''.join(random.choices(bases, k=100))  # 100bp reads
    qual = ''.join(random.choices(quals, k=100))  # Good quality scores
    return f"@test_read{{read_num}}\\n{{seq}}\\n+\\n{{qual}}\\n".format(
        read_num=read_num, seq=seq, qual=qual)

output_dir = "{output}"  # Snakemake output directory

# Create R1 file
with gzip.open(os.path.join(output_dir, "test_R1.fastq.gz"), 'wt') as f:
    for i in range(10000):  # 10k reads
        f.write(make_read(i))

# Create R2 file
with gzip.open(os.path.join(output_dir, "test_R2.fastq.gz"), 'wt') as f:
    for i in range(10000):  # 10k reads
        f.write(make_read(i))
EOF
            else
                echo "Using existing FASTQ files in {output}"
            fi
        else
            echo "Creating input directory and test files in {output}"
            mkdir -p {output}
            # Use the same Python script as above with slight modifications
            python3 - <<'EOF'
import gzip
import random
import os

bases = ['A', 'C', 'G', 'T']
quals = ''.join([chr(x) for x in range(33, 74)])  # Phred+33 quality scores

def make_read(read_num):
    seq = ''.join(random.choices(bases, k=100))  # 100bp reads
    qual = ''.join(random.choices(quals, k=100))  # Good quality scores
    return f"@test_read{{read_num}}\\n{{seq}}\\n+\\n{{qual}}\\n".format(
        read_num=read_num, seq=seq, qual=qual)

output_dir = "{output}"  # Snakemake output directory

# Create R1 file
with gzip.open(os.path.join(output_dir, "test_R1.fastq.gz"), 'wt') as f:
    for i in range(10000):  # 10k reads
        f.write(make_read(i))

# Create R2 file
with gzip.open(os.path.join(output_dir, "test_R2.fastq.gz"), 'wt') as f:
    for i in range(10000):  # 10k reads
        f.write(make_read(i))
EOF
        fi
        """

rule run_random_iterations:
    input:
        os.path.join(INPUT_DIR),
        r1=os.path.join(INPUT_DIR, "{sample}_R1.fastq.gz"),
        r2=os.path.join(INPUT_DIR, "{sample}_R2.fastq.gz")
    output:
        outdir=directory(os.path.join(OUTPUT_DIR, "random_iterations/{sample}")),
        quality=os.path.join(OUTPUT_DIR, "results/quality_results_{sample}.tsv"),
        combination=os.path.join(OUTPUT_DIR, "results/combination_results_{sample}.tsv"),
        trimmomatic_params=os.path.join(OUTPUT_DIR, "params/trimmomatic_params_{sample}.tsv"),
        fastp_params=os.path.join(OUTPUT_DIR, "params/fastp_params_{sample}.tsv"),
        cutadapt_params=os.path.join(OUTPUT_DIR, "params/cutadapt_params_{sample}.tsv"),
        bbduk_params=os.path.join(OUTPUT_DIR, "params/bbduk_params_{sample}.tsv"),
        sickle_params=os.path.join(OUTPUT_DIR, "params/sickle_params_{sample}.tsv")
    conda:
        "trimming_env"
    params:
        iterations=ITERATIONS,
        output_dir=OUTPUT_DIR,
        params_dir=os.path.join(OUTPUT_DIR, "params"),
        results_dir=os.path.join(OUTPUT_DIR, "results")
    shell:
        """
        mkdir -p {params.params_dir}
        mkdir -p {params.results_dir}
        mkdir -p {output.outdir}
        python trim_randomizer.py \
            --input-dir $(dirname {input.r1}) \
            --output-dir {output.outdir} \
            --iterations {params.iterations} \
            --quality-output {output.quality} \
            --combination-output {output.combination} \
            --params-dir {params.params_dir}
        """

def aggregate_input(wildcards):
    checkpoint_output = checkpoints.find_samples.get(**wildcards).output[0]
    samples = glob_wildcards(os.path.join(checkpoint_output, "{sample}_R1.fastq.gz")).sample
    quality_files = expand(os.path.join(OUTPUT_DIR, "results/quality_results_{sample}.tsv"), sample=samples)
    combination_files = expand(os.path.join(OUTPUT_DIR, "results/combination_results_{sample}.tsv"), sample=samples)
    return {"quality": quality_files, "combination": combination_files}

rule combine_results:
    input:
        unpack(aggregate_input)
    output:
        quality=os.path.join(OUTPUT_DIR, "combined/quality_results.tsv"),
        combination=os.path.join(OUTPUT_DIR, "combined/combination_results.tsv")
    params:
        output_dir=os.path.join(OUTPUT_DIR, "combined")
    shell:
        """
        mkdir -p {params.output_dir}
        if [ -n "$(ls -A {input.quality} 2>/dev/null)" ]; then
            head -n 1 "$(ls {input.quality[0]})" > {output.quality}
            tail -n +2 -q {input.quality} >> {output.quality}
        else
            echo "No quality results found"
            touch {output.quality}
        fi
        
        if [ -n "$(ls -A {input.combination} 2>/dev/null)" ]; then
            head -n 1 "$(ls {input.combination[0]})" > {output.combination}
            tail -n +2 -q {input.combination} >> {output.combination}
        else
            echo "No combination results found"
            touch {output.combination}
        fi
        """

rule combine_params:
    input:
        trimmomatic=expand(os.path.join(OUTPUT_DIR, "params/trimmomatic_params_{sample}.tsv"), sample=SAMPLES),
        fastp=expand(os.path.join(OUTPUT_DIR, "params/fastp_params_{sample}.tsv"), sample=SAMPLES),
        cutadapt=expand(os.path.join(OUTPUT_DIR, "params/cutadapt_params_{sample}.tsv"), sample=SAMPLES),
        bbduk=expand(os.path.join(OUTPUT_DIR, "params/bbduk_params_{sample}.tsv"), sample=SAMPLES),
        sickle=expand(os.path.join(OUTPUT_DIR, "params/sickle_params_{sample}.tsv"), sample=SAMPLES)
    output:
        trimmomatic=os.path.join(OUTPUT_DIR, "params/trimmomatic_params.tsv"),
        fastp=os.path.join(OUTPUT_DIR, "params/fastp_params.tsv"),
        cutadapt=os.path.join(OUTPUT_DIR, "params/cutadapt_params.tsv"),
        bbduk=os.path.join(OUTPUT_DIR, "params/bbduk_params.tsv"),
        sickle=os.path.join(OUTPUT_DIR, "params/sickle_params.tsv")
    shell:
        """
        # Combine trimmomatic parameters
        head -n 1 {input.trimmomatic[0]} > {output.trimmomatic}
        tail -n +2 -q {input.trimmomatic} >> {output.trimmomatic}
        
        # Combine fastp parameters
        head -n 1 {input.fastp[0]} > {output.fastp}
        tail -n +2 -q {input.fastp} >> {output.fastp}
        
        # Combine cutadapt parameters
        head -n 1 {input.cutadapt[0]} > {output.cutadapt}
        tail -n +2 -q {input.cutadapt} >> {output.cutadapt}
        
        # Combine bbduk parameters
        head -n 1 {input.bbduk[0]} > {output.bbduk}
        tail -n +2 -q {input.bbduk} >> {output.bbduk}
        
        # Combine sickle parameters
        head -n 1 {input.sickle[0]} > {output.sickle}
        tail -n +2 -q {input.sickle} >> {output.sickle}
        """

rule train_ml_model:
    input:
        iterations=expand(os.path.join(OUTPUT_DIR, "random_iterations/{sample}"), sample=SAMPLES),
        quality=os.path.join(OUTPUT_DIR, "combined/quality_results.tsv"),
        combination=os.path.join(OUTPUT_DIR, "combined/combination_results.tsv"),
        params_trimmomatic=os.path.join(OUTPUT_DIR, "params/trimmomatic_params.tsv"),
        params_fastp=os.path.join(OUTPUT_DIR, "params/fastp_params.tsv"),
        params_cutadapt=os.path.join(OUTPUT_DIR, "params/cutadapt_params.tsv"),
        params_bbduk=os.path.join(OUTPUT_DIR, "params/bbduk_params.tsv"),
        params_sickle=os.path.join(OUTPUT_DIR, "params/sickle_params.tsv")
    output:
        models=directory(os.path.join(OUTPUT_DIR, "trained_models")),
        optimal_params=os.path.join(OUTPUT_DIR, "optimal_parameters.json"),
        combination=os.path.join(OUTPUT_DIR, "optimal_combination.txt")
    shell:
        """
        python ml_optimizer.py \
            --output-dir {output.models} \
            --params-output {output.optimal_params} \
            --combination-output {output.combination} \
            --quality-input {input.quality} \
            --combination-input {input.combination} \
            --params-dir {OUTPUT_DIR}/params
        """

rule apply_optimal_trimming:
    input:
        r1=os.path.join(INPUT_DIR, "{sample}_R1.fastq.gz"),
        r2=os.path.join(INPUT_DIR, "{sample}_R2.fastq.gz"),
        params=os.path.join(OUTPUT_DIR, "optimal_parameters.json"),
        combination=os.path.join(OUTPUT_DIR, "optimal_combination.txt")
    output:
        r1=os.path.join(OUTPUT_DIR, "optimized_trimming/{sample}_R1.fastq.gz"),
        r2=os.path.join(OUTPUT_DIR, "optimized_trimming/{sample}_R2.fastq.gz")
    conda:
        "trimming_env"
    shell:
        """
        python trim_randomizer.py \
            --input-dir $(dirname {input.r1}) \
            --output-dir $(dirname {output.r1}) \
            --iterations 1 \
            --use-ml
        """

rule generate_best_trimmer_summary:
    input:
        params=os.path.join(OUTPUT_DIR, "optimal_parameters.json"),
        combination=os.path.join(OUTPUT_DIR, "optimal_combination.txt"),
        quality_results=os.path.join(OUTPUT_DIR, "combined/quality_results.tsv"),
        combination_results=os.path.join(OUTPUT_DIR, "combined/combination_results.tsv")
    output:
        summary=os.path.join(OUTPUT_DIR, "best_trimmer_summary.tsv"),
        plots=directory(os.path.join(OUTPUT_DIR, "plots"))
    script:
        "scripts/summarize_best_trimmer.py"

rule create_parameter_visualization:
    input:
        os.path.join(OUTPUT_DIR, "best_trimmer_summary.tsv")
    output:
        os.path.join(OUTPUT_DIR, "report/parameter_visualization.html")
    shell:
        """
        mkdir -p $(dirname {output})
        python scripts/visualize_parameters.py \
            --input {input} \
            --output {output}
        """

rule export_best_parameters:
    input:
        os.path.join(OUTPUT_DIR, "optimal_parameters.json")
    output:
        config=os.path.join(OUTPUT_DIR, "config/best_trimmer_config.yaml"),
        script=os.path.join(OUTPUT_DIR, "scripts/run_best_trimmer.sh")
    shell:
        """
        mkdir -p $(dirname {output.config})
        mkdir -p $(dirname {output.script})
        python scripts/export_best_parameters.py \
            --input {input} \
            --config-output {output.config} \
            --script-output {output.script}
        chmod +x {output.script}
        """

rule generate_report:
    input:
        params=os.path.join(OUTPUT_DIR, "optimal_parameters.json"),
        combination=os.path.join(OUTPUT_DIR, "optimal_combination.txt"),
        quality_results=os.path.join(OUTPUT_DIR, "combined/quality_results.tsv"),
        combination_results=os.path.join(OUTPUT_DIR, "combined/combination_results.tsv"),
        summary=os.path.join(OUTPUT_DIR, "best_trimmer_summary.tsv"),
        parameter_viz=os.path.join(OUTPUT_DIR, "report/parameter_visualization.html")
    output:
        report=os.path.join(OUTPUT_DIR, "report/best_trimmer_report.html")
    script:
        "scripts/generate_report.py"
