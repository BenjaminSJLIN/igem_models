import numpy as np

def calculate_activated_cas12a(C_total, T_total, K_D):
    """
    Calculate the concentration of Activated Cas12a (A) given total Cas12a (C_total), 
    total Target DNA (T_total), and binding affinity (K_D).
    Based on quadratic solution to: K_D = ([C_total] - [A])([T_total] - [A]) / [A]
    """
    b = C_total + T_total + K_D
    c = C_total * T_total
    
    # Quadratic formula: A^2 - bA + c = 0
    # We take the negative root because A cannot exceed C_total or T_total
    A = (b - np.sqrt(b**2 - 4 * c)) / 2.0
    return A

def calculate_cleavage_rate(A, S, k_cat, K_M):
    """
    Calculate the initial cleavage rate (v) using Michaelis-Menten kinetics.
    v = (k_cat * A * S) / (K_M + S)
    """
    return (k_cat * A * S) / (K_M + S)

def calculate_LOD(C_total, S, P_threshold, t_incubate, k_cat, K_M, K_D):
    """
    Calculate the Limit of Detection (LOD) which is the Target DNA concentration (T_total)
    required to achieve a specific cleaved product concentration (P_threshold) 
    after a given incubation time (t_incubate).
    """
    # Required average cleavage rate to reach P_threshold in t_incubate
    v_threshold = P_threshold / t_incubate
    
    # 1. Find the required amount of Activated Cas12a (A)
    # v_threshold = (k_cat * A * S) / (K_M + S)
    A_req = (v_threshold * (K_M + S)) / (k_cat * S)
    
    # If the required activated Cas12a is greater than the total available Cas12a,
    # the threshold is unachievable (LOD is infinite/NaN)
    if A_req >= C_total:
        return np.nan
        
    # 2. Find the required total Target DNA to achieve A_req
    # K_D = (C_total - A_req)(T_total - A_req) / A_req
    T_total = A_req + (K_D * A_req) / (C_total - A_req)
    
    return T_total
