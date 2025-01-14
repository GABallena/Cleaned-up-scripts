import pandas as pd
import numpy as np
from sklearn.ensemble import RandomForestRegressor
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import train_test_split
import joblib

class TrimmingOptimizer:
    def __init__(self):
        self.models = {}
        self.scalers = {}
        self.best_params = {}
        self.best_combination = None
    
    def load_parameter_data(self):
        # Load parameters for each tool
        params_data = {
            'trimmomatic': pd.read_csv('trimmomatic_params.tsv', sep='\t'),
            'fastp': pd.read_csv('fastp_params.tsv', sep='\t'),
            'cutadapt': pd.read_csv('cutadapt_params.tsv', sep='\t'),
            'bbduk': pd.read_csv('bbduk_params.tsv', sep='\t'),
            'sickle': pd.read_csv('sickle_params.tsv', sep='\t')
        }
        return params_data

    def load_results(self):
        return pd.read_csv('quality_results.tsv', sep='\t')

    def load_combination_data(self):
        """Load the combination results"""
        return pd.read_csv('combination_results.tsv', sep='\t')

    def prepare_features(self, tool_params, tool):
        if tool == 'trimmomatic':
            return tool_params[['MinLen', 'SlidingWindow']].values
        elif tool == 'fastp':
            return tool_params[['QualityCut', 'LengthCut', 'NBaseLimit']].values
        elif tool == 'cutadapt':
            return tool_params[['ErrorRate', 'MinimumLength', 'Overlap']].values
        elif tool == 'bbduk':
            return tool_params[['K', 'MinK', 'HDist', 'MinLen', 'TrimQ', 'MinKmerHits', 'MinKmerFraction']].values
        elif tool == 'sickle':
            return tool_params[['QualityThreshold', 'LengthThreshold']].values

    def calculate_score(self, row):
        # Custom scoring function that prioritizes minimizing entropy/richness loss
        # while maximizing quality gain
        return (row['Quality_Gain'] * 0.4 - 
                row['Entropy_Loss'] * 0.3 - 
                row['Richness_Loss'] * 0.3)

    def calculate_combination_score(self, row):
        """Calculate score for a combination of tools"""
        return (row['Quality_Gain'] * 0.4 - 
                row['Entropy_Loss'] * 0.3 - 
                row['Richness_Loss'] * 0.3)

    def train_models(self):
        params_data = self.load_parameter_data()
        results = self.load_results()
        
        for tool in params_data.keys():
            # Merge parameters with results
            tool_data = pd.merge(params_data[tool], results[results['Tool'] == tool],
                               on=['Iteration', 'Sample'])
            
            # Prepare features and target
            X = self.prepare_features(tool_data, tool)
            tool_data['Score'] = tool_data.apply(self.calculate_score, axis=1)
            y = tool_data['Score'].values
            
            # Scale features
            scaler = StandardScaler()
            X_scaled = scaler.fit_transform(X)
            
            # Train model
            model = RandomForestRegressor(n_estimators=100, random_state=42)
            model.fit(X_scaled, y)
            
            # Store model and scaler
            self.models[tool] = model
            self.scalers[tool] = scaler
            
            # Find best parameters
            self.best_params[tool] = self.find_best_parameters(model, scaler, tool)
        
        # Also find the best combination
        self.find_best_combination()
            
    def find_best_parameters(self, model, scaler, tool):
        # Use random search to find best parameters
        n_trials = 1000
        best_score = float('-inf')
        best_params = None
        
        param_ranges = self.get_parameter_ranges(tool)
        
        for _ in range(n_trials):
            params = {k: np.random.uniform(v[0], v[1]) for k, v in param_ranges.items()}
            X_trial = np.array([[params[k] for k in param_ranges.keys()]])
            X_scaled = scaler.transform(X_trial)
            score = model.predict(X_scaled)[0]
            
            if score > best_score:
                best_score = score
                best_params = params
                
        return best_params

    def get_parameter_ranges(self, tool):
        if tool == 'trimmomatic':
            return {
                'MinLen': (30, 50),
                'SlidingWindow': (4, 10)
            }
        # Add ranges for other tools...
        return {}

    def find_best_combination(self):
        """Find the best performing combination of tools"""
        combo_data = self.load_combination_data()
        
        # Calculate scores for each combination
        combo_data['Score'] = combo_data.apply(self.calculate_combination_score, axis=1)
        
        # Group by combination and calculate mean score
        combo_scores = combo_data.groupby('Tool_Combination')['Score'].mean()
        
        # Get the best combination
        best_combo = combo_scores.idxmax()
        self.best_combination = best_combo.split(',')
        
        return self.best_combination

    def save_models(self):
        for tool in self.models:
            joblib.dump(self.models[tool], f'{tool}_model.joblib')
            joblib.dump(self.scalers[tool], f'{tool}_scaler.joblib')

    def get_optimal_parameters(self):
        return self.best_params

    def get_optimal_combination(self):
        """Return the optimal combination of tools"""
        if self.best_combination is None:
            self.find_best_combination()
        return self.best_combination

def main():
    optimizer = TrimmingOptimizer()
    optimizer.train_models()
    optimizer.save_models()
    
    optimal_params = optimizer.get_optimal_parameters()
    print("Optimal parameters for each tool:")
    for tool, params in optimal_params.items():
        print(f"\n{tool.upper()}:")
        for param, value in params.items():
            print(f"{param}: {value}")
    
    optimal_combination = optimizer.get_optimal_combination()
    print("\nOptimal combination of tools:")
    print(", ".join(optimal_combination))

if __name__ == "__main__":
    main()
