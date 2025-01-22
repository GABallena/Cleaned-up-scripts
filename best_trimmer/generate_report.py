import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from jinja2 import Template
import os
import numpy as np
from scipy import stats

class QualityReport:
    def __init__(self, csv_file):
        self.df = pd.read_csv(csv_file)
        self.output_dir = "report_output"
        os.makedirs(self.output_dir, exist_ok=True)
        
    def generate_summary_stats(self):
        """Generate summary statistics for each metric by tool"""
        metrics = ['Entropy', 'GC Content', 'Sequence Length', 'N Content']
        summary = self.df.groupby('Tool')[metrics].agg(['mean', 'std', 'min', 'max'])
        return summary
    
    def plot_metric_distributions(self):
        """Create violin plots with jittered points for key metrics across different tools"""
        metrics = ['Entropy', 'GC Content', 'Sequence Length', 'N Content']
        figs = []
        
        for metric in metrics:
            plt.figure(figsize=(12, 7))
            
            # Create violin plot
            sns.violinplot(x='Tool', y=metric, data=self.df, 
                         inner=None, color='lightgray')
            
            # Add jittered points
            sns.stripplot(x='Tool', y=metric, data=self.df,
                        size=4, alpha=0.3, jitter=0.2)
            
            plt.xticks(rotation=45)
            plt.title(f'{metric} Distribution by Tool')
            plt.xlabel('Tool')
            plt.ylabel(metric)
            plt.tight_layout()
            
            # Save plot
            filepath = os.path.join(self.output_dir, f'{metric.lower().replace(" ", "_")}_violin.png')
            plt.savefig(filepath, dpi=300, bbox_inches='tight')
            plt.close()
            figs.append(filepath)
        
        return figs

    def analyze_base_composition(self):
        """Analyze and plot base composition across tools"""
        base_cols = [col for col in self.df.columns if col.startswith('Base_') and col.endswith('_Percentage')]
        base_data = self.df.groupby('Tool')[base_cols].mean()
        
        plt.figure(figsize=(12, 6))
        base_data.plot(kind='bar', stacked=True)
        plt.title('Average Base Composition by Tool')
        plt.xlabel('Tool')
        plt.ylabel('Percentage')
        plt.legend(title='Base', bbox_to_anchor=(1.05, 1))
        plt.tight_layout()
        
        filepath = os.path.join(self.output_dir, 'base_composition.png')
        plt.savefig(filepath)
        plt.close()
        
        return filepath, base_data

    def statistical_analysis(self):
        """Perform statistical tests between tools"""
        results = {}
        metrics = ['Entropy', 'GC Content', 'Sequence Length', 'N Content']
        
        for metric in metrics:
            # Perform Kruskal-Wallis H-test
            tools = self.df['Tool'].unique()
            groups = [self.df[self.df['Tool'] == tool][metric] for tool in tools]
            h_stat, p_value = stats.kruskal(*groups)
            
            results[metric] = {
                'h_statistic': h_stat,
                'p_value': p_value,
                'significant': p_value < 0.05
            }
        
        return results

    def generate_html_report(self):
        """Generate an HTML report with all analyses"""
        # Get all analysis results
        summary_stats = self.generate_summary_stats()
        plot_files = self.plot_metric_distributions()
        base_plot, base_data = self.analyze_base_composition()
        stat_results = self.statistical_analysis()
        
        # HTML template
        template_str = """
        <html>
        <head>
            <title>Quality Metrics Report</title>
            <style>
                body { font-family: Arial, sans-serif; margin: 20px; }
                h1, h2 { color: #333; }
                .plot { margin: 20px 0; }
                table { border-collapse: collapse; width: 100%; }
                th, td { border: 1px solid #ddd; padding: 8px; text-align: left; }
                th { background-color: #f2f2f2; }
            </style>
        </head>
        <body>
            <h1>Quality Metrics Report</h1>
            
            <h2>Summary Statistics</h2>
            {{ summary_stats.to_html() }}
            
            <h2>Metric Distributions</h2>
            {% for plot in plot_files %}
            <div class="plot">
                <img src="{{ plot }}" alt="Distribution Plot">
            </div>
            {% endfor %}
            
            <h2>Base Composition</h2>
            <div class="plot">
                <img src="{{ base_plot }}" alt="Base Composition">
            </div>
            
            <h2>Statistical Analysis</h2>
            <table>
                <tr>
                    <th>Metric</th>
                    <th>H-statistic</th>
                    <th>p-value</th>
                    <th>Significant Difference</th>
                </tr>
                {% for metric, results in stat_results.items() %}
                <tr>
                    <td>{{ metric }}</td>
                    <td>{{ "%.4f"|format(results.h_statistic) }}</td>
                    <td>{{ "%.4f"|format(results.p_value) }}</td>
                    <td>{{ "Yes" if results.significant else "No" }}</td>
                </tr>
                {% endfor %}
            </table>
        </body>
        </html>
        """
        
        # Render template
        template = Template(template_str)
        html_content = template.render(
            summary_stats=summary_stats,
            plot_files=[os.path.basename(p) for p in plot_files],
            base_plot=os.path.basename(base_plot),
            stat_results=stat_results
        )
        
        # Save HTML report
        report_path = os.path.join(self.output_dir, 'quality_report.html')
        with open(report_path, 'w') as f:
            f.write(html_content)
        
        return report_path

def main():
    # Create report
    report = QualityReport('read_metrics.csv')
    
    # Generate report
    report_path = report.generate_html_report()
    print(f"Report generated: {report_path}")

if __name__ == "__main__":
    main()
