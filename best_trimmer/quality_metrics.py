import numpy as np
from scipy.stats import entropy
from Bio import SeqIO
from Bio.SeqUtils import gc_fraction  # Updated import
import glob
import os
from collections import Counter
from itertools import product
import statistics
import gzip

def calculate_entropy(sequence):
    """Calculate Shannon entropy of a sequence."""
    _, counts = np.unique(list(sequence), return_counts=True)
    return entropy(counts)

def calculate_quality_score(fastq_file):
    """Calculate quality metrics for a FASTQ file (handles gzipped files)"""
    quality_scores = []
    try:
        # Open gzipped files
        if str(fastq_file).endswith('.gz'):
            with gzip.open(fastq_file, 'rt') as handle:
                for record in SeqIO.parse(handle, "fastq"):
                    quality_scores.append(sum(record.letter_annotations["phred_quality"]) / len(record.seq))
        # Open regular files
        else:
            with open(fastq_file, 'r') as handle:
                for record in SeqIO.parse(handle, "fastq"):
                    quality_scores.append(sum(record.letter_annotations["phred_quality"]) / len(record.seq))
    except Exception as e:
        print(f"Error processing {fastq_file}: {str(e)}")
        raise

    return {
        'quality_scores': quality_scores,
        'avg_quality': sum(quality_scores) / len(quality_scores) if quality_scores else 0
    }

def calculate_gc_content(sequence):
    """Calculate GC content of a sequence."""
    return gc_fraction(sequence) * 100  # Updated GC calculation

def calculate_sequence_complexity(sequence):
    """Calculate linguistic sequence complexity."""
    total_kmers = len(sequence)
    unique_kmers = len(set(sequence[i:i+3] for i in range(len(sequence)-2)))
    return unique_kmers / total_kmers if total_kmers > 0 else 0

def calculate_n50(lengths):
    """Calculate N50 of sequence lengths."""
    sorted_lengths = sorted(lengths, reverse=True)
    total = sum(sorted_lengths)
    running_sum = 0
    for length in sorted_lengths:
        running_sum += length
        if running_sum >= total/2:
            return length
    return 0

def calculate_dinucleotide_freq(sequence):
    """Calculate dinucleotide frequencies."""
    if len(sequence) < 2:
        return 0
    dinucs = [''.join(p) for p in product('ATGC', repeat=2)]
    freq = Counter([sequence[i:i+2] for i in range(len(sequence)-1)])
    return sum(abs(freq[d]/sum(freq.values()) - 0.0625) for d in dinucs)

def calculate_quality_distribution(record):
    """Analyze quality score distribution."""
    quals = record.letter_annotations["phred_quality"]
    return {
        'mean': statistics.mean(quals),
        'std': statistics.stdev(quals) if len(quals) > 1 else 0,
        'median': statistics.median(quals),
        'q1': np.percentile(quals, 25),
        'q3': np.percentile(quals, 75)
    }

def calculate_base_balance(sequence):
    """Calculate base composition balance."""
    counts = Counter(sequence)
    total = sum(counts.values())
    expected = total/4  # Expected count for perfect balance
    return 1 - sum(abs(counts[base]/total - 0.25) for base in 'ATGC')/2

def calculate_extended_biological_metrics(record):
    """Calculate comprehensive biological metrics."""
    sequence = str(record.seq)
    qual_dist = calculate_quality_distribution(record)
    
    metrics = {
        'gc_content': calculate_gc_content(sequence),
        'complexity': calculate_sequence_complexity(sequence),
        'length': len(sequence),
        'dinuc_bias': calculate_dinucleotide_freq(sequence),
        'base_balance': calculate_base_balance(sequence),
        'quality_mean': qual_dist['mean'],
        'quality_std': qual_dist['std'],
        'quality_consistency': 1 - (qual_dist['std'] / qual_dist['mean'] if qual_dist['mean'] > 0 else 1)
    }
    return metrics

def parse_trim_results(results_dir):
    """Parse trim_randomizer outputs with extended metrics."""
    results = []
    for output_file in glob.glob(os.path.join(results_dir, "*.fastq")):
        params = os.path.basename(output_file).split('_')
        
        # Initialize aggregate metrics
        metrics_sum = {
            'entropy': 0,
            'gc_content': 0,
            'complexity': 0,
            'dinuc_bias': 0,
            'base_balance': 0,
            'quality_mean': 0,
            'quality_consistency': 0
        }
        lengths = []
        read_count = 0
        
        reads = list(SeqIO.parse(output_file, "fastq"))
        for record in reads:
            bio_metrics = calculate_extended_biological_metrics(record)
            metrics_sum['entropy'] += calculate_entropy(str(record.seq))
            lengths.append(len(record.seq))
            
            for key in ['gc_content', 'complexity', 'dinuc_bias', 
                       'base_balance', 'quality_mean', 'quality_consistency']:
                metrics_sum[key] += bio_metrics[key] if key in bio_metrics else 0
            
            read_count += 1
            
        # Calculate averages
        if read_count > 0:
            for key in metrics_sum:
                metrics_sum[key] /= read_count
                
        metrics_sum['n50'] = calculate_n50(lengths)
        
        results.append({
            'file': output_file,
            'params': params,
            **metrics_sum
        })
    
    return results

def find_optimal_trim(results_dir, weights=None):
    """Find optimal parameters with extended metrics."""
    if weights is None:
        weights = {
            'entropy': 0.15,
            'gc_content': 0.15,
            'complexity': 0.15,
            'dinuc_bias': 0.1,
            'base_balance': 0.15,
            'quality_mean': 0.15,
            'quality_consistency': 0.1,
            'n50': 0.05
        }
    
    results = parse_trim_results(results_dir)
    
    # Normalize scores
    max_scores = {
        metric: max(r[metric] for r in results)
        for metric in ['entropy', 'gc_content', 'complexity', 'dinuc_bias', 'base_balance', 'quality_mean', 'quality_consistency', 'n50']
    }
    
    for r in results:
        r['total_score'] = sum(
            weights[metric] * (r[metric] / max_scores[metric])
            for metric in weights.keys()
        )
    
    return max(results, key=lambda x: x['total_score'])

if __name__ == "__main__":
    import argparse
    
    parser = argparse.ArgumentParser(description='Find optimal trimming parameters')
    parser.add_argument('results_dir', help='Directory containing trim_randomizer outputs')
    parser.add_argument('--weights', type=str, default=None,
                       help='Comma-separated weights for metrics (entropy,gc,complexity,dinuc,balance,qual_mean,qual_cons,n50)')
    
    args = parser.parse_args()
    
    weights = None
    if args.weights:
        weight_values = [float(x) for x in args.weights.split(',')]
        weights = {
            'entropy': weight_values[0],
            'gc_content': weight_values[1],
            'complexity': weight_values[2],
            'dinuc_bias': weight_values[3],
            'base_balance': weight_values[4],
            'quality_mean': weight_values[5],
            'quality_consistency': weight_values[6],
            'n50': weight_values[7]
        }
    
    best = find_optimal_trim(args.results_dir, weights)
    
    print(f"Best trimming parameters found:")
    print(f"File: {best['file']}")
    print(f"Parameters: {best['params']}")
    for metric in ['entropy', 'gc_content', 'complexity', 'dinuc_bias', 'base_balance', 'quality_mean', 'quality_consistency', 'n50', 'total_score']:
        print(f"{metric.replace('_', ' ').title()}: {best[metric]:.3f}")
