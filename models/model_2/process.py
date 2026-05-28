import numpy as np
import pandas as pd
from scipy.integrate import odeint
from .equations import rcr_dynamics
import yaml
import os

def run(config, input_data=None):
    """
    Run simulation for Model 2.
    input_data is expected to be the result dataframe from Model 1.0.
    """
    params = config['parameters']
    sim_config = config['simulation']
    
    t_end = sim_config['t_end']
    t_steps = sim_config['t_steps']
    
    t = np.linspace(0, t_end, t_steps)
    initial_conditions = [0.0] # [S_0]
    
    # Load k_cat and K_M from config.yaml since they are defined in model_1_0
    config_path = os.path.join(os.path.dirname(__file__), '..', '..', 'config.yaml')
    with open(config_path, 'r', encoding='utf-8') as f:
        full_config = yaml.safe_load(f)
    k_cat = full_config['model_1_0']['parameters']['k_cat']
    K_M = full_config['model_1_0']['parameters']['K_M']
    
    results = []
    
    # We want to simulate for different Target_DNA and Ara_Concentration
    # Filter for 'Activation' type to get Activated_Cas12a values from input_data
    if input_data is not None and 'Analysis_Type' in input_data.columns:
        activation_df = input_data[input_data['Analysis_Type'] == 'Activation']
    else:
        # Fallback dummy data if no valid input_data is provided
        activation_df = pd.DataFrame({
            'Ara_Concentration': [10.0],
            'Target_DNA': [10.0],
            'Activated_Cas12a': [0.5]
        })
        
    # To avoid too many plots, take the max Ara_concentration
    max_ara = activation_df['Ara_Concentration'].max()
    filtered_df = activation_df[activation_df['Ara_Concentration'] == max_ara]
    
    # Select a few representative Target DNA levels to simulate
    target_dna_levels = [0.1, 1.0, 10.0, 100.0]
    
    for t_dna in target_dna_levels:
        # Find the closest Target_DNA in the filtered_df
        closest_row = filtered_df.iloc[(filtered_df['Target_DNA'] - t_dna).abs().argsort()[:1]]
        if closest_row.empty:
            continue
        A = closest_row['Activated_Cas12a'].values[0]
        actual_t_dna = closest_row['Target_DNA'].values[0]
        
        solution = odeint(
            rcr_dynamics,
            initial_conditions,
            t,
            args=(A, k_cat, K_M, params)
        )
        
        S_sol = solution[:, 0]
        
        df_temp = pd.DataFrame({
            'Time': t,
            'Ara_Concentration': max_ara,
            'Target_DNA': actual_t_dna,
            'Activated_Cas12a': A,
            'ssDNA': S_sol
        })
        results.append(df_temp)
        
    if results:
        final_df = pd.concat(results, ignore_index=True)
    else:
        final_df = pd.DataFrame()
        
    return final_df
