import numpy as np

def rcr_dynamics(y, t, A, k_cat, K_M, params):
    """
    System of ODEs for Rolling Circle Replication (RCR).
    y[0] = S (ssDNA concentration)
    """
    S = y[0]
    
    alpha_RCR = params['alpha_RCR']
    P_s = params['P_s']
    G = params['G']
    K_G = params['K_G']
    delta_s = params['delta_s']
    
    # Calculate v_cleavage dynamically based on current S
    v_cleavage = (k_cat * A * S) / (K_M + S)
    
    # dS/dt = alpha_RCR * P_s * (G / (K_G + G)) - v_cleavage - delta_s * S
    dS_dt = alpha_RCR * P_s * (G / (K_G + G)) - v_cleavage - delta_s * S
    
    return [dS_dt]
