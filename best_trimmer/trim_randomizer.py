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
    adapter_fwd = 'ADAPTER_FWD'
    adapter_rev = 'ADAPTER_REV'
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

def run_trimmomatic(input_file, output_file, iterations):
    for i in range(iterations):
        params = random_trimmomatic_parameters(i, os.path.basename(input_file))
        command = [
            'trimmomatic', 'PE', '-phred33',
            input_file, input_file.replace('_R1_', '_R2_'),
            output_file + f'_paired_{i}_R1.fastq', output_file + f'_unpaired_{i}_R1.fastq',
            output_file + f'_paired_{i}_R2.fastq', output_file + f'_unpaired_{i}_R2.fastq',
            f'adapters.fa:{params[0]}:{params[1]}:{params[2]}',
            f'LEADING:{params[3]}', f'TRAILING:{params[4]}',
            f'SLIDINGWINDOW:{params[5]}:{params[6]}', f'MINLEN:{params[7]}'
        ]
        print(f"Running iteration {i+1} with parameters: {params}")
        subprocess.run(command)

def process_directory(input_dir, output_dir, iterations):
    for file in os.listdir(input_dir):
        if file.endswith('_1.fastq.gz'):
            input_file_1 = os.path.join(input_dir, file)
            input_file_2 = input_file_1.replace('_1.fastq.gz', '_2.fastq.gz')
            output_file_base = os.path.join(output_dir, os.path.splitext(file)[0].replace('_1', ''))
            run_trimmomatic(input_file_1, output_file_base, iterations)
            for i in range(iterations):
                sample = os.path.basename(input_file_1)
                params = random_cutadapt_parameters(i, sample)
                # Run Cutadapt
                command = [
                    'cutadapt',
                    '-a', params[0], '-A', params[1],
                    '-m', str(params[3]),
                    '-q', f'{params[2]}',
                    '-o', output_file_base + f'_cutadapt_{i}_R1.fastq',
                    '-p', output_file_base + f'_cutadapt_{i}_R2.fastq',
                    input_file_1, input_file_2
                ]
                print(f"Running Cutadapt iteration {i+1} with parameters: {params}")
                subprocess.run(command)

                params = random_bbduk_parameters(i, sample)
                # Run BBDuk
                command = [
                    'bbduk.sh',
                    f'in1={input_file_1}', f'in2={input_file_2}',
                    f'out1={output_file_base}_bbduk_{i}_R1.fastq', f'out2={output_file_base}_bbduk_{i}_R2.fastq',
                    f'ref=adapters.fa', f'ktrim={params[0]}', f'k={params[1]}', f'mink={params[2]}',
                    f'hdist={params[3]}', f'trimq={params[4]}', f'qtrim={params[5]}',
                    f'minlen={params[6]}'
                ]
                print(f"Running BBDuk iteration {i+1} with parameters: {params}")
                subprocess.run(command)

                params = random_fastp_parameters(i, sample)
                # Run Fastp
                command = [
                    'fastp',
                    '-i', input_file_1, '-I', input_file_2,
                    '-o', output_file_base + f'_fastp_{i}_R1.fastq', '-O', output_file_base + f'_fastp_{i}_R2.fastq',
                    '--detect_adapter_for_pe', '--trim_front1', str(params[0]), '--trim_tail1', str(params[1]),
                    '--cut_front', '--cut_tail', '--cut_window_size', str(params[2]), '--cut_mean_quality', str(params[3]),
                    '--length_required', str(params[4])
                ]
                print(f"Running Fastp iteration {i+1} with parameters: {params}")
                subprocess.run(command)

                params = random_sickle_parameters(i, sample)
                # Run Sickle
                command = [
                    'sickle', 'pe',
                    '-f', input_file_1, '-r', input_file_2,
                    '-t', 'sanger',
                    '-o', output_file_base + f'_sickle_{i}_R1.fastq', '-p', output_file_base + f'_sickle_{i}_R2.fastq',
                    '-s', output_file_base + f'_sickle_{i}_singles.fastq',
                    '-q', str(params[0]), '-l', str(params[1])
                ]
                print(f"Running Sickle iteration {i+1} with parameters: {params}")
                subprocess.run(command)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Process paired-end reads with various trimming tools.')
    parser.add_argument('-i', '--input_dir', type=str, required=True, help='Directory containing input paired-end reads')
    parser.add_argument('-o', '--output_dir', type=str, required=True, help='Directory to save output files')
    parser.add_argument('-n', '--iterations', type=int, default=10, help='Number of iterations to run for each tool')
    
    args = parser.parse_args()
    
    process_directory(args.input_dir, args.output_dir, args.iterations)