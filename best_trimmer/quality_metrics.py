from Bio import SeqIO
from Bio.SeqUtils import gc_fraction
import os
import numpy as np
import math
import csv
import argparse

def calculate_shannon_entropy(sequence):
    """Calculate the Shannon entropy of a sequence."""
    bases = ['A', 'C', 'G', 'T']
    entropy = 0
    for base in bases:
        p = sequence.count(base) / len(sequence)
        if p > 0:
            entropy -= p * math.log2(p)
    return entropy

def calculate_n_content(sequence):
    """Calculate the N content of a sequence."""
    return sequence.count('N') / len(sequence)

def calculate_base_composition(records):
    """Calculate base composition per position."""
    base_counts = {'A': [], 'C': [], 'G': [], 'T': [], 'N': []}
    for record in records:
        for i, base in enumerate(record.seq):
            if len(base_counts['A']) <= i:
                for key in base_counts:
                    base_counts[key].append(0)
            base_counts[base][i] += 1
    return base_counts

def calculate_per_sequence_gc_content(sequence):
    """Calculate the GC content for a sequence."""
    return gc_fraction(sequence)

def calculate_per_base_n_content(records):
    """Calculate the proportion of 'N' bases at each position."""
    n_content = []
    for record in records:
        for i, base in enumerate(record.seq):
            if len(n_content) <= i:
                n_content.append(0)
            if base == 'N':
                n_content[i] += 1
    return [count / len(records) for count in n_content]

def calculate_per_sequence_n_content(sequence):
    """Calculate the proportion of 'N' bases in a sequence."""
    return sequence.count('N') / len(sequence)

def identify_overrepresented_sequences(sequence_counts, threshold=0.01):
    """Identify sequences that are overrepresented in the dataset."""
    total_sequences = sum(sequence_counts.values())
    return {seq: count for seq, count in sequence_counts.items() if count / total_sequences > threshold}

def detect_adapter_content(sequence, adapters):
    """Detect and quantify the presence of adapter sequences."""
    adapter_counts = {adapter: 0 for adapter in adapters}
    for adapter in adapters:
        adapter_counts[adapter] += sequence.count(adapter)
    return adapter_counts

def calculate_base_percentages(sequence):
    """Calculate the percentages of A, C, G, and T in a sequence."""
    length = len(sequence)
    return {
        'A': (sequence.count('A') / length) * 100,
        'C': (sequence.count('C') / length) * 100,
        'G': (sequence.count('G') / length) * 100,
        'T': (sequence.count('T') / length) * 100
    }

def calculate_at_content(sequence):
    """Calculate the AT content of a sequence."""
    return (sequence.count('A') + sequence.count('T')) / len(sequence)

def calculate_sequence_complexity(sequence):
    """Calculate the complexity of a sequence."""
    return len(set(sequence)) / len(sequence)

def calculate_dinucleotide_frequencies(sequence):
    """Calculate the frequencies of dinucleotides in a sequence."""
    dinucleotides = ['AA', 'AC', 'AG', 'AT', 'CA', 'CC', 'CG', 'CT', 'GA', 'GC', 'GG', 'GT', 'TA', 'TC', 'TG', 'TT']
    frequencies = {dinucleotide: 0 for dinucleotide in dinucleotides}
    for i in range(len(sequence) - 1):
        dinucleotide = sequence[i:i+2]
        if dinucleotide in frequencies:
            frequencies[dinucleotide] += 1
    total_dinucleotides = len(sequence) - 1
    return {dinucleotide: count / total_dinucleotides for dinucleotide, count in frequencies.items()}

def process_fasta(file_path, adapters=[]):
    """Process a FASTA file and gather quality metrics."""
    metrics_list = []

    try:
        with open(file_path, "rt") as handle:
            records = list(SeqIO.parse(handle, "fasta"))
    except Exception as e:
        print(f"Error processing file {file_path}: {e}")
        return []

    print(f"Number of records in {file_path}: {len(records)}")  # Debug print
    if len(records) > 0:
        print(f"First record in {file_path}: {records[0].format('fasta')}")  # Debug print

    for record in records:
        seq_str = str(record.seq)
        metrics = {
            "entropy": calculate_shannon_entropy(seq_str),
            "gc_content": calculate_per_sequence_gc_content(record.seq),
            "at_content": calculate_at_content(seq_str),
            "sequence_length": len(record.seq),
            "n_content": calculate_per_sequence_n_content(seq_str),
            "n_percentage": calculate_n_content(seq_str) * 100,
            "base_percentages": calculate_base_percentages(seq_str),
            "sequence_complexity": calculate_sequence_complexity(seq_str),
            "dinucleotide_frequencies": calculate_dinucleotide_frequencies(seq_str),
            "adapter_content": detect_adapter_content(seq_str, adapters)
        }
        metrics_list.append((record.id, metrics))

    return metrics_list

def write_metrics_to_csv(metrics_list, output_file, sample_name):
    """Write the metrics to a CSV file."""
    file_exists = os.path.isfile(output_file)
    if not file_exists:
        with open(output_file, mode='w', newline='') as file:
            writer = csv.writer(file)
            header = [
                "Sample Name", "Read ID", "Entropy", "GC Content", 
                "AT Content", "Sequence Length", "N Content", "N Percentage", 
                "Base Percentages", "Sequence Complexity", "Dinucleotide Frequencies", 
                "Adapter Content"
            ]
            writer.writerow(header)
    with open(output_file, mode='a', newline='') as file:
        writer = csv.writer(file)
        for read_id, metrics in metrics_list:
            print(f"Writing data for {sample_name}, read {read_id}: {metrics}")  # Debug print
            writer.writerow([
                sample_name, read_id, metrics['entropy'], metrics['gc_content'], 
                metrics['at_content'], metrics['sequence_length'], metrics['n_content'], 
                metrics['n_percentage'], metrics['base_percentages'], 
                metrics['sequence_complexity'], metrics['dinucleotide_frequencies'], 
                metrics['adapter_content']
            ])

def read_adapters(adapter_file):
    """Read adapter sequences from a file."""
    adapters = []
    with open(adapter_file, 'r') as file:
        for line in file:
            if not line.startswith('>'):
                adapters.append(line.strip())
    return adapters

def main(input_directory, adapter_file):
    """Main function to process all FASTA files in a directory."""
    output_file = "read_metrics.csv"  # Hardcoded output file name
    # Ensure the output file is created if it doesn't exist
    if not os.path.isfile(output_file):
        with open(output_file, mode='w', newline='') as file:
            writer = csv.writer(file)
            header = [
                "Sample Name", "Read ID", "Entropy", "GC Content", 
                "AT Content", "Sequence Length", "N Content", "N Percentage", 
                "Base Percentages", "Sequence Complexity", "Dinucleotide Frequencies", 
                "Adapter Content"
            ]
            writer.writerow(header)
    
    adapters = read_adapters(adapter_file)
    print(f"Adapters: {adapters}")  # Debug print
    for filename in os.listdir(input_directory):
        # Skip FASTQ files and only process FASTA files
        if filename.endswith(('.fastq', '.fq')):
            print(f"Skipping FASTQ file: {filename}")
            continue
        if filename.endswith(('.fasta', '.fa')):
            file_path = os.path.join(input_directory, filename)
            print(f"Processing FASTA file: {file_path}")  # Debug print
            metrics_list = process_fasta(file_path, adapters)
            print(f"Results for {filename}: {metrics_list}")  # Debug print
            write_metrics_to_csv(metrics_list, output_file, filename)
            print(f"Metrics for {filename} have been written to {output_file}")  # Debug print

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process FASTA files and gather quality metrics.")
    parser.add_argument("-i", "--input", required=True, help="Input directory containing FASTA files")
    parser.add_argument("-a", "--adapters", required=True, help="File containing adapter sequences")
    
    args = parser.parse_args()
    main(args.input, args.adapters)
    print("Metrics have been written to read_metrics.csv")