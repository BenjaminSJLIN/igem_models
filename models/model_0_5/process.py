import numpy as np
import pandas as pd
from scipy.integrate import odeint
from .equations import cas12a_dynamics

def run(config):
    """
    Run simulation for Model 0.5 across different Arabinose concentrations.
    """
    params = config['parameters']
    sim_config = config['simulation']
    
    t_end = sim_config['t_end']
    t_steps = sim_config['t_steps']
    ara_concentrations = sim_config['ara_concentrations']
    
    t = np.linspace(0, t_end, t_steps)
    initial_conditions = [0.0, 0.0] # [mRNA_0, Cas12a_0]
    
    results = []
    
    for ara in ara_concentrations:
        # Solve ODE
        solution = odeint(
            cas12a_dynamics, 
            initial_conditions, 
            t, 
            args=(ara, params)
        )
        
        # Store results
        mRNA_sol = solution[:, 0]
        Cas12a_sol = solution[:, 1]
        
        df_temp = pd.DataFrame({
            'Time': t,
            'Ara_Concentration': ara,
            'mRNA': mRNA_sol,
            'Cas12a': Cas12a_sol
        })
        results.append(df_temp)
        
    final_df = pd.concat(results, ignore_index=True)
    return final_df
