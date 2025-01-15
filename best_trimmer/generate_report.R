import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
from pathlib import Path
import json
from datetime import datetime

def load_data(snakemake):
    """Load all required input data files"""
    with open(snakemake.input.params) as f:
        optimal_params = json.load(f)
    
    with open(snakemake.input.combination) as f:
        optimal_combination = f.read().strip().split(',')
    
    quality_results = pd.read_csv(snakemake.input.quality_results, sep='\t')
    combination_results = pd.read_csv(snakemake.input.combination_results, sep='\t')
    summary = pd.read_csv(snakemake.input.summary, sep='\t')
    
    return optimal_params, optimal_combination, quality_results, combination_results, summary

def create_quality_plot(quality_results):
    """Create quality comparison plot"""
    fig = px.box(quality_results, 
                 x='Tool', 
                 y=['Quality_Gain', 'Entropy_Loss', 'Richness_Loss'],
                 title='Quality Metrics by Trimming Tool')
    return fig

def create_combination_plot(combination_results):
    """Create tool combination performance plot"""
    combo_summary = combination_results.groupby('Tool_Combination').agg({
        'Quality_Gain': 'mean',
        'Entropy_Loss': 'mean',
        'Richness_Loss': 'mean'
    }).reset_index()
    
    fig = px.parallel_coordinates(combo_summary,
                                dimensions=['Quality_Gain', 'Entropy_Loss', 'Richness_Loss'],
                                color='Tool_Combination',
                                title='Tool Combination Performance')
    return fig

def generate_parameter_summary(optimal_params):
    """Generate HTML table of optimal parameters"""
    html = "<h3>Optimal Parameters</h3>"
    html += "<table border='1'>"
    html += "<tr><th>Tool</th><th>Parameter</th><th>Value</th></tr>"
    
    for tool, params in optimal_params.items():
        for param, value in params.items():
            html += f"<tr><td>{tool}</td><td>{param}</td><td>{value}</td></tr>"
    
    html += "</table>"
    return html

def generate_html_report(snakemake):
    """Generate the final HTML report"""
    # Load all data
    optimal_params, optimal_combination, quality_results, combination_results, summary = load_data(snakemake)
    
    # Create plots
    quality_plot = create_quality_plot(quality_results)
    combination_plot = create_combination_plot(combination_results)
    
    # Generate HTML content
    html_content = f"""
    <!DOCTYPE html>
    <html>
    <head>
        <title>Best Trimmer Analysis Report</title>
        <script src="https://cdn.plot.ly/plotly-latest.min.js"></script>
        <style>
            body {{ font-family: Arial, sans-serif; margin: 40px; }}
            .container {{ max-width: 1200px; margin: auto; }}
            .section {{ margin-bottom: 40px; }}
        </style>
    </head>
    <body>
        <div class="container">
            <h1>Best Trimmer Analysis Report</h1>
            <p>Generated on: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}</p>
            
            <div class="section">
                <h2>Optimal Tool Combination</h2>
                <p>The best performing combination of tools is:</p>
                <ul>
                    {''.join(f'<li>{tool}</li>' for tool in optimal_combination)}
                </ul>
            </div>
            
            <div class="section">
                {generate_parameter_summary(optimal_params)}
            </div>
            
            <div class="section">
                <h2>Quality Metrics Comparison</h2>
                {quality_plot.to_html(full_html=False)}
            </div>
            
            <div class="section">
                <h2>Tool Combination Performance</h2>
                {combination_plot.to_html(full_html=False)}
            </div>
            
            <div class="section">
                <h2>Performance Summary</h2>
                {summary.to_html(index=False)}
            </div>
        </div>
    </body>
    </html>
    """
    
    # Write the report
    with open(snakemake.output.report, 'w') as f:
        f.write(html_content)

if __name__ == "__main__":
    generate_html_report(snakemake)
