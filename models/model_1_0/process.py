import numpy as np
import pandas as pd
from .equations import calculate_activated_cas12a, calculate_cleavage_rate, calculate_LOD

def run(config, input_data=None):
    """
    Run simulation for Model 1.0.
    input_data is expected to be the result dataframe from Model 0.5.
    """
    params = config['parameters']
    sim_config = config['simulation']
    
    K_D = params['K_D']
    k_cat = params['k_cat']
    K_M = params['K_M']
    P_threshold = params['P_threshold']
    t_incubate = params['t_incubate']
    
    # Generate log-spaced arrays for sweeping
    t_dna_vals = np.logspace(
        np.log10(sim_config['target_dna_min']), 
        np.log10(sim_config['target_dna_max']), 
        sim_config['target_dna_steps']
    )
    
    ssdna_vals = np.logspace(
        np.log10(sim_config['ssdna_min']), 
        np.log10(sim_config['ssdna_max']), 
        sim_config['ssdna_steps']
    )
    
    # 1. Extract steady-state Cas12a from input_data
    if input_data is not None:
        steady_state_df = input_data.groupby('Ara_Concentration').last().reset_index()
    else:
        # Fallback if no input data is provided (e.g. testing standalone)
        steady_state_df = pd.DataFrame({
            'Ara_Concentration': [0, 1.0, 10.0],
            'Cas12a': [0.01, 0.5, 1.0]
        })
        
    results = []
    
    # Analysis 1: Activation vs Target DNA (Fix ssDNA=max, sweep Target DNA)
    for _, row in steady_state_df.iterrows():
        ara = row['Ara_Concentration']
        c_total = row['Cas12a']
        
        for t_dna in t_dna_vals:
            A = calculate_activated_cas12a(c_total, t_dna, K_D)
            # We can also calculate a nominal cleavage rate at saturation ssDNA
            v = calculate_cleavage_rate(A, sim_config['ssdna_max'], k_cat, K_M)
            
            results.append({
                'Analysis_Type': 'Activation',
                'Ara_Concentration': ara,
                'Cas12a_Total': c_total,
                'Target_DNA': t_dna,
                'ssDNA': sim_config['ssdna_max'],
                'Activated_Cas12a': A,
                'Cleavage_Rate': v,
                'LOD': np.nan
            })
            
    # Analysis 2: Cleavage Rate vs ssDNA (Fix Cas12a at max Arabinose, sweep ssDNA for a few Target DNA levels)
    max_ara_row = steady_state_df.loc[steady_state_df['Ara_Concentration'].idxmax()]
    max_c_total = max_ara_row['Cas12a']
    max_ara = max_ara_row['Ara_Concentration']
    
    # Pick a few representative Target DNA concentrations
    t_dna_levels = [0.1, 1.0, 10.0, 100.0] 
    
    for t_dna in t_dna_levels:
        A = calculate_activated_cas12a(max_c_total, t_dna, K_D)
        for s in ssdna_vals:
            v = calculate_cleavage_rate(A, s, k_cat, K_M)
            results.append({
                'Analysis_Type': 'Cleavage',
                'Ara_Concentration': max_ara,
                'Cas12a_Total': max_c_total,
                'Target_DNA': t_dna,
                'ssDNA': s,
                'Activated_Cas12a': A,
                'Cleavage_Rate': v,
                'LOD': np.nan
            })
            
    # Analysis 3: LOD vs Arabinose Concentration
    for _, row in steady_state_df.iterrows():
        ara = row['Ara_Concentration']
        c_total = row['Cas12a']
        
        lod = calculate_LOD(c_total, sim_config['ssdna_max'], P_threshold, t_incubate, k_cat, K_M, K_D)
        
        results.append({
            'Analysis_Type': 'LOD',
            'Ara_Concentration': ara,
            'Cas12a_Total': c_total,
            'Target_DNA': np.nan,
            'ssDNA': sim_config['ssdna_max'],
            'Activated_Cas12a': np.nan,
            'Cleavage_Rate': np.nan,
            'LOD': lod
        })
        
    final_df = pd.DataFrame(results)
    return final_df
