import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from scipy import stats
import os
from ast import literal_eval
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from sklearn.manifold import TSNE

def load_and_clean_data(file_path):
    """Load and prepare the data for analysis, with improved error handling."""
    print(f"Loading data from {file_path}")
    df = pd.read_csv(file_path)
    print(f"Loaded {len(df)} rows from CSV")
    
    # Debug: Print first few sample names
    print("First few sample names:")
    print(df['Sample Name'].head())
    
    print("Extracting trimmer names...")
    
    def extract_trimmer(sample_name):
        try:
            if pd.isna(sample_name) or not isinstance(sample_name, str):
                return None
                
            # Convert to uppercase for case-insensitive matching
            name_upper = sample_name.upper()
            
            # First try direct pattern matching from filename structure
            if '_' in sample_name:
                parts = sample_name.split('_')
                for part in parts:
                    part_upper = part.upper()
                    if part_upper in ['BBDUK', 'CUTADAPT', 'FASTP', 'SICKLE', 'TRIMMOMATIC']:
                        return part_upper
                    elif 'BB' in part_upper and 'DUK' in part_upper:  # Handle variations of BBDUK
                        return 'BBDUK'
            
            # Fallback to more complex pattern matching if needed
            trimmer_patterns = {
                'SICKLE': ['SICKLE'],
                'FASTP': ['FASTP'],
                'TRIMMOMATIC': ['TRIMMOMATIC'],
                'CUTADAPT': ['CUTADAPT'],
                'BBDUK': ['BBDUK', 'BB_DUK', 'BB-DUK', 'BBDUK']
            }
            
            for trimmer, patterns in trimmer_patterns.items():
                if any(pattern in name_upper for pattern in patterns):
                    return trimmer
            
            return None
        except Exception as e:
            print(f"Error processing sample name {sample_name}: {str(e)}")
            return None
    
    # Extract trimmer names
    df['Trimmer'] = df['Sample Name'].apply(extract_trimmer)
    
    # Debug: Print unique trimmer names
    print("\nUnique trimmer names found:")
    print(df['Trimmer'].value_counts())
    
    # Parse dictionary strings
    print("\nParsing data structures...")
    try:
        df['Base Percentages'] = df['Base Percentages'].apply(lambda x: literal_eval(x) if pd.notnull(x) else {})
        df['Dinucleotide Frequencies'] = df['Dinucleotide Frequencies'].apply(lambda x: literal_eval(x) if pd.notnull(x) else {})
        df['Adapter Content'] = df['Adapter Content'].apply(lambda x: literal_eval(x) if pd.notnull(x) else {})
    except Exception as e:
        print(f"Error parsing dictionaries: {str(e)}")
        raise
    
    # Remove rows with invalid trimmer names
    df = df.dropna(subset=['Trimmer'])
    print(f"\nFinal number of valid rows: {len(df)}")
    
    return df

def generate_basic_stats(df):
    numeric_cols = ['Entropy', 'GC.Content', 'AT.Content', 'Sequence.Length', 
                    'N.Content', 'N.Percentage']
    stats = df[numeric_cols].describe()
    return stats

def generate_trimmer_stats(df):
    trimmer_stats = df.groupby('Trimmer').agg(
        Mean_Length=('Sequence.Length', 'mean'),
        Mean_GC=('GC.Content', 'mean'),
        Mean_Entropy=('Entropy', 'mean'),
        Mean_N_Percent=('N.Percentage', 'mean'),
        Sequences_Count=('Sequence.Length', 'count')
    ).round(2)
    
    return trimmer_stats

def setup_plotting_style():
    sns.set(style="whitegrid")
    plt.rcParams.update({
        'figure.figsize': (12, 8),
        'axes.titlesize': 14,
        'axes.labelsize': 12,
        'xtick.labelsize': 10,
        'ytick.labelsize': 10,
        'legend.title_fontsize': 12,
        'legend.fontsize': 10
    })

def save_plot(plot, filename, output_dir):
    filepath = os.path.join(output_dir, f"{filename}.png")
    plot.figure.savefig(filepath, dpi=600)
    plt.close(plot.figure)

def plot_length_distribution(df, output_dir):
    plt.figure()
    sns.violinplot(y='Sequence.Length', data=df, inner=None, color="lightgray")
    sns.boxplot(y='Sequence.Length', data=df, whis=np.inf, color="white")
    plt.axhline(df['Sequence.Length'].mean(), color='red', linestyle='dashed')
    plt.axhline(df['Sequence.Length'].median(), color='blue', linestyle='dashed')
    plt.title("Sequence Length Distribution")
    plt.ylabel("Length (bp)")
    save_plot(plt, "length_distribution", output_dir)

def plot_gc_content(df, output_dir):
    plt.figure()
    sns.histplot(df['GC.Content'], kde=True, color="lightblue")
    plt.axvline(50, color='red', linestyle='dashed')
    plt.axvline(df['GC.Content'].mean(), color='green', linestyle='solid')
    plt.title("GC Content Distribution")
    plt.xlabel("GC Content (%)")
    save_plot(plt, "gc_content_distribution", output_dir)

def plot_dimensionality_reduction(df, output_dir):
    features = ['Entropy', 'GC.Content', 'AT.Content', 'Sequence.Length', 'N.Percentage']
    X = StandardScaler().fit_transform(df[features])
    
    pca = PCA(n_components=2)
    pca_result = pca.fit_transform(X)
    df['PC1'] = pca_result[:, 0]
    df['PC2'] = pca_result[:, 1]
    
    plt.figure()
    sns.scatterplot(x='PC1', y='PC2', hue='Trimmer', data=df, alpha=0.6)
    plt.title("PCA of Trimming Results")
    save_plot(plt, "pca_trimming_results", output_dir)
    
    tsne = TSNE(n_components=2, perplexity=30)
    tsne_result = tsne.fit_transform(X)
    df['TSNE1'] = tsne_result[:, 0]
    df['TSNE2'] = tsne_result[:, 1]
    
    plt.figure()
    sns.scatterplot(x='TSNE1', y='TSNE2', hue='Trimmer', data=df, alpha=0.6)
    plt.title("t-SNE of Trimming Results")
    save_plot(plt, "tsne_trimming_results", output_dir)

def generate_html_report(stats, trimmer_stats, output_dir):
    html_content = f'''
    <!DOCTYPE html>
    <html>
    <head>
        <title>Trimmer Comparison Report</title>
        <style>
            body {{ font-family: Arial, sans-serif; margin: 20px; }}
            table {{ border-collapse: collapse; width: 100%; }}
            th, td {{ border: 1px solid #ddd; padding: 8px; }}
            th {{ background-color: #f2f2f2; }}
        </style>
    </head>
    <body>
        <h1>Trimmer Comparison Report</h1>
        <h2>Overall Statistics</h2>
        {stats.to_html()}
        <h2>Trimmer Statistics</h2>
        {trimmer_stats.to_html()}
        <div class="plots">
            <h2>Analysis Plots</h2>
            <img src="length_distribution.png" alt="Length Distribution">
            <img src="gc_content_distribution.png" alt="GC Content Distribution">
            <img src="pca_trimming_results.png" alt="PCA of Trimming Results">
            <img src="tsne_trimming_results.png" alt="t-SNE of Trimming Results">
        </div>
    </body>
    </html>
    '''
    
    with open(os.path.join(output_dir, "trimmer_comparison_report.html"), 'w') as f:
        f.write(html_content)

def main():
    # Create output directory
    output_dir = 'quality_reports'
    os.makedirs(output_dir, exist_ok=True)
    
    try:
        # Load and process data with improved error handling
        df = load_and_clean_data("read_metrics.csv")
        
        if len(df) == 0:
            raise ValueError("No valid data rows after processing")
        
        # Generate statistics and plots
        stats = generate_basic_stats(df)
        trimmer_stats = generate_trimmer_stats(df)
        
        # Setup plotting
        setup_plotting_style()
        
        # Generate plots
        plot_length_distribution(df, output_dir)
        plot_gc_content(df, output_dir)
        plot_dimensionality_reduction(df, output_dir)
        
        # Generate HTML report
        generate_html_report(stats, trimmer_stats, output_dir)
        print(f"Report generated successfully in {output_dir}/trimmer_comparison_report.html")
        
    except Exception as e:
        print(f"Error generating report: {str(e)}")
        raise

if __name__ == "__main__":
    main()
