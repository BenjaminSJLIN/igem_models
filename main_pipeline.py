import yaml
import os
import sys

# Ensure the parent directory is in the path so we can import from models
sys.path.append(os.path.dirname(os.path.abspath(__file__)))

from models.model_0_5 import process as process_0_5
from utils.plot_utils import Plotter

def load_config(filepath="config.yaml"):
    with open(filepath, 'r', encoding='utf-8') as f:
        return yaml.safe_load(f)

def run_pipeline():
    print("--- Loading Configuration ---")
    config = load_config()
    model_0_5_config = config.get('model_0_5', {})
    
    print("\n--- Running Stage 0.5 (Cas12a Activation) ---")
    # Execute Model 0.5
    result_0_5 = process_0_5.run(model_0_5_config)
    
    # Save intermediate data
    intermediate_dir = os.path.join("data", "02_intermediate")
    if not os.path.exists(intermediate_dir):
        os.makedirs(intermediate_dir)
        
    intermediate_path = os.path.join(intermediate_dir, "model_0_5_result.csv")
    result_0_5.to_csv(intermediate_path, index=False)
    print(f"Saved intermediate data to: {intermediate_path}")
    
    print("\n--- Generating Plots for Stage 0.5 ---")
    # Initialize plotter for Model 0.5
    plotter_0_5 = Plotter(output_dir=os.path.join("outputs", "model_0_5_plots"))
    
    # Time-series plot for Cas12a
    plotter_0_5.plot_time_series(
        data=result_0_5, 
        y_col="Cas12a", 
        title="Cas12a Expression Dynamics over Time", 
        filename="cas12a_time_series.png"
    )
    
    # Time-series plot for mRNA (bonus)
    plotter_0_5.plot_time_series(
        data=result_0_5, 
        y_col="mRNA", 
        title="Cas12a mRNA Dynamics over Time", 
        filename="mRNA_time_series.png"
    )
    
    # Dose-Response curve for Cas12a
    plotter_0_5.plot_dose_response(
        data=result_0_5, 
        y_col="Cas12a", 
        title="Steady State Cas12a vs Arabinose Concentration", 
        filename="cas12a_dose_response.png"
    )
    
    print("\n--- Pipeline Execution Completed! ---")

if __name__ == "__main__":
    run_pipeline()
