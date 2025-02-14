import os
import random
import subprocess
import argparse

generated_parameters = {}

def random_trimmomatic_parameters(iteration, sample):
    seed_mismatches = random.randint(0, 2)
    palindrome_clip_threshold = random.randint(20, 40)
    simple_clip_threshold = random.randint(7, 15)
    leading = random.randint(3, 20)
    trailing = random.randint(3, 20)
    slidingwindow_size = random.randint(4, 10)
    slidingwindow_quality = random.randint(15, 30)
    minlen = random.randint(36, 100)
    
    param_key = f"{sample}_trimmomatic_{iteration}"
    generated_parameters[param_key] = (seed_mismatches, palindrome_clip_threshold, simple_clip_threshold, 
                                       leading, trailing, slidingwindow_size, slidingwindow_quality, minlen)
    
    return generated_parameters[param_key]

def random_fastp_parameters(iteration, sample):
    trim_front1 = random.randint(3, 20)
    trim_tail1 = random.randint(3, 20)
    cut_window_size = random.randint(4, 10)
    cut_mean_quality = random.randint(15, 30)
    length_required = random.randint(36, 100)
    
    param_key = f"{sample}_fastp_{iteration}"
    generated_parameters[param_key] = (trim_front1, trim_tail1, cut_window_size, cut_mean_quality, length_required)
    
    return generated_parameters[param_key]

def random_cutadapt_parameters(iteration, sample):
    adapter_fwd = 'AGATCGGAAGAGCACACGTCTGAACTCCAGTCA'
    adapter_rev = 'AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT'
    minlen = random.randint(36, 100)
    quality = random.randint(15, 30)
    
    param_key = f"{sample}_cutadapt_{iteration}"
    generated_parameters[param_key] = (adapter_fwd, adapter_rev, quality, minlen)
    
    return generated_parameters[param_key]

def random_bbduk_parameters(iteration, sample):
    ktrim = 'r'
    k = random.randint(20, 40)
    mink = random.randint(7, 15)
    hdist = random.randint(0, 2)
    trimq = random.randint(15, 30)
    qtrim = 'rl'
    minlen = random.randint(36, 100)
    
    param_key = f"{sample}_bbduk_{iteration}"
    generated_parameters[param_key] = (ktrim, k, mink, hdist, trimq, qtrim, minlen)
    
    return generated_parameters[param_key]

def random_sickle_parameters(iteration, sample):
    quality = random.randint(15, 30)
    length_threshold = random.randint(36, 100)
    
    param_key = f"{sample}_sickle_{iteration}"
    generated_parameters[param_key] = (quality, length_threshold)
    
    return generated_parameters[param_key]

def remove_singles_and_unpaired(output_dir):
    """Remove all singles and unpaired samples from the output directory."""
    for file in os.listdir(output_dir):
        if 'singles' in file or 'unpaired' in file:
            file_path = os.path.join(output_dir, file)
            print(f"Removing file: {file_path}")  # Debug print
            os.remove(file_path)

def process_directory(input_dir, output_dir, iterations, log_file):
    with open(log_file, 'w') as f:
        f.write("Tool\tIteration\tSample\tParameters\n")
        f.flush()  # Ensure header is written immediately
        for file in os.listdir(input_dir):
            if file.endswith('_1.fastq.gz'):
                input_file_1 = os.path.join(input_dir, file)
                input_file_2 = input_file_1.replace('_1.fastq.gz', '_2.fastq.gz')
                sample_base = os.path.splitext(file)[0].replace('_1', '')
                for i in range(iterations):
                    sample = os.path.basename(input_file_1)
                    
                    # Run Trimmomatic
                    params = random_trimmomatic_parameters(i, sample)
                    f.write(f"Trimmomatic\t{i}\t{sample}\t{params}\n")
                    f.flush()  # Ensure log is written immediately
                    command = [
                        'trimmomatic', 'PE', '-phred33',
                        input_file_1, input_file_2,
                        os.path.join(output_dir, f'{sample_base}_trimmomatic_{i}_paired_R1.fastq'), 
                        os.path.join(output_dir, f'{sample_base}_trimmomatic_{i}_unpaired_R1.fastq'),
                        os.path.join(output_dir, f'{sample_base}_trimmomatic_{i}_paired_R2.fastq'), 
                        os.path.join(output_dir, f'{sample_base}_trimmomatic_{i}_unpaired_R2.fastq'),
                        'ILLUMINACLIP:adapters.fa:{0}:{1}:{2}'.format(params[0], params[1], params[2]),
                        f'LEADING:{params[3]}', f'TRAILING:{params[4]}',
                        f'SLIDINGWINDOW:{params[5]}:{params[6]}', f'MINLEN:{params[7]}'
                    ]
                    print(f"Running Trimmomatic iteration {i+1} with parameters: {params}")
                    subprocess.run(command)
                    
                    # Run Cutadapt
                    params = random_cutadapt_parameters(i, sample)
                    f.write(f"Cutadapt\t{i}\t{sample}\t{params}\n")
                    f.flush()  # Ensure log is written immediately
                    command = [
                        'cutadapt',
                        '-a', params[0], '-A', params[1],
                        '-m', str(params[3]),
                        '-q', f'{params[2]}',
                        '-o', os.path.join(output_dir, f'{sample_base}_cutadapt_{i}_R1.fastq'),
                        '-p', os.path.join(output_dir, f'{sample_base}_cutadapt_{i}_R2.fastq'),
                        input_file_1, input_file_2
                    ]
                    print(f"Running Cutadapt iteration {i+1} with parameters: {params}")
                    subprocess.run(command)

                    # Run BBDuk
                    params = random_bbduk_parameters(i, sample)
                    f.write(f"BBDuk\t{i}\t{sample}\t{params}\n")
                    f.flush()  # Ensure log is written immediately
                    command = [
                        'bbduk.sh',
                        f'in1={input_file_1}', f'in2={input_file_2}',
                        f'out1={os.path.join(output_dir, f"{sample_base}_bbduk_{i}_R1.fastq")}', 
                        f'out2={os.path.join(output_dir, f"{sample_base}_bbduk_{i}_R2.fastq")}',
                        f'ref=adapters.fa', f'ktrim={params[0]}', f'k={params[1]}', f'mink={params[2]}',
                        f'hdist={params[3]}', f'trimq={params[4]}', f'qtrim={params[5]}',
                        f'minlen={params[6]}'
                    ]
                    print(f"Running BBDuk iteration {i+1} with parameters: {params}")
                    subprocess.run(command)

                    # Run Fastp
                    params = random_fastp_parameters(i, sample)
                    f.write(f"Fastp\t{i}\t{sample}\t{params}\n")
                    f.flush()  # Ensure log is written immediately
                    command = [
                        'fastp',
                        '-i', input_file_1, '-I', input_file_2,
                        '-o', os.path.join(output_dir, f'{sample_base}_fastp_{i}_R1.fastq'), 
                        '-O', os.path.join(output_dir, f'{sample_base}_fastp_{i}_R2.fastq'),
                        '--detect_adapter_for_pe', '--trim_front1', str(params[0]), '--trim_tail1', str(params[1]),
                        '--cut_front', '--cut_tail', '--cut_window_size', str(params[2]), '--cut_mean_quality', str(params[3]),
                        '--length_required', str(params[4])
                    ]
                    print(f"Running Fastp iteration {i+1} with parameters: {params}")
                    subprocess.run(command)

                    # Run Sickle
                    params = random_sickle_parameters(i, sample)
                    f.write(f"Sickle\t{i}\t{sample}\t{params}\n")
                    f.flush()  # Ensure log is written immediately
                    # Create temporary directory for fastq files
                    temp_dir = os.path.join(output_dir, 'temp_fastq')
                    os.makedirs(temp_dir, exist_ok=True)
                    input_file_1_fastq = os.path.join(temp_dir, os.path.basename(input_file_1).replace('.fastq.gz', '.fastq'))
                    input_file_2_fastq = os.path.join(temp_dir, os.path.basename(input_file_2).replace('.fastq.gz', '.fastq'))
                    # Stream the outputs
                    subprocess.run(f'gunzip -c {input_file_1} > {input_file_1_fastq}', shell=True)
                    subprocess.run(f'gunzip -c {input_file_2} > {input_file_2_fastq}', shell=True)
                    command = [
                        'sickle', 'pe',
                        '-f', input_file_1_fastq, '-r', input_file_2_fastq,
                        '-t', 'sanger',
                        '-o', os.path.join(output_dir, f'{sample_base}_sickle_{i}_R1.fastq'), 
                        '-p', os.path.join(output_dir, f'{sample_base}_sickle_{i}_R2.fastq'),
                        '-s', os.path.join(output_dir, f'{sample_base}_sickle_{i}_singles.fastq'),
                        '-q', str(params[0]), '-l', str(params[1])
                    ]
                    print(f"Running Sickle iteration {i+1} with parameters: {params}")
                    subprocess.run(command)
                    # Remove temporary fastq files
                    os.remove(input_file_1_fastq)
                    os.remove(input_file_2_fastq)
                    os.rmdir(temp_dir)
    
    # Remove singles and unpaired files after processing
    remove_singles_and_unpaired(output_dir)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Process paired-end reads with various trimming tools.')
    parser.add_argument('-i', '--input_dir', type=str, required=True, help='Directory containing input paired-end reads')
    parser.add_argument('-o', '--output_dir', type=str, required=True, help='Directory to save output files')
    parser.add_argument('-n', '--iterations', type=int, default=10, help='Number of iterations to run for each tool')
    parser.add_argument('-l', '--log_file', type=str, required=True, help='File to log the parameters used')
    
    args = parser.parse_args()
    
    process_directory(args.input_dir, args.output_dir, args.iterations, args.log_file)