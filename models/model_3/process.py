import numpy as np
import pandas as pd
from scipy.integrate import odeint
from .equations import plasmid_dynamics

def run(config, input_data=None):
    """
    Run simulation for Model 3 (Plasmid degradation).
    input_data is expected to be the result dataframe from Model 1.0.
    """
    params = config['parameters']
    sim_config = config['simulation']
    
    t_end = sim_config['t_end']
    t_steps = sim_config['t_steps']
    P_s_initial = sim_config['P_s_initial']
    
    t = np.linspace(0, t_end, t_steps)
    initial_conditions = [P_s_initial]
    
    results = []
    
    # We want to simulate for different Target_DNA
    # Filter for 'Activation' type to get Cleavage_Rate values from input_data
    if input_data is not None and 'Analysis_Type' in input_data.columns:
        activation_df = input_data[input_data['Analysis_Type'] == 'Activation']
    else:
        # Fallback dummy data if no valid input_data is provided
        activation_df = pd.DataFrame({
            'Ara_Concentration': [10.0],
            'Target_DNA': [10.0],
            'Cleavage_Rate': [50.0]
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
        v_cleavage = closest_row['Cleavage_Rate'].values[0]
        actual_t_dna = closest_row['Target_DNA'].values[0]
        
        solution = odeint(
            plasmid_dynamics,
            initial_conditions,
            t,
            args=(v_cleavage, params)
        )
        
        Ps_sol = solution[:, 0]
        
        df_temp = pd.DataFrame({
            'Time': t,
            'Ara_Concentration': max_ara,
            'Target_DNA': actual_t_dna,
            'v_cleavage': v_cleavage,
            'Plasmid_Concentration': Ps_sol
        })
        results.append(df_temp)
        
    if results:
        final_df = pd.concat(results, ignore_index=True)
    else:
        final_df = pd.DataFrame()
        
    return final_df
