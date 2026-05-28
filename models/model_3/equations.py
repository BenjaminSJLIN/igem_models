import numpy as np

def plasmid_dynamics(y, t, v_cleavage, params):
    """
    System of ODEs for Sensor Plasmid (Ps) degradation.
    y[0] = P_s (Sensor Plasmid concentration)
    """
    P_s = y[0]
    
    r = params['r']
    mu = params['mu']
    k_loss = params['k_loss']
    
    # dP_s/dt = r*P_s - mu*P_s - k_loss*v_cleavage*P_s
    dPs_dt = r * P_s - mu * P_s - k_loss * v_cleavage * P_s
    
    return [dPs_dt]
