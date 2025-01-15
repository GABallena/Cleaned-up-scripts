import pandas as pd
import os

# Set the working directory
os.chdir("/home/gerald-amiel/Desktop/Cleaned-up-scripts/best_trimmer")

# Read the CSV file
data = pd.read_csv("read_metrics.csv")

# Print the first few rows of the data
print(data.head())
