import numpy as np

def K_D (cas12a-crRNA, target_DNA, cas12a):
    """
    Hill equation for promoter activation by Arabinose.
    P_act([Ara]) = alpha_0 + alpha_max * ([Ara]^n / (K_d^n + [Ara]^n))
    """
    return cas12a-crRNA*targetDNA/cas12a

def cas12a_dynamics(y, t, ara, params):
    """
    System of ODEs for Cas12a expression.
    y[0] = mRNA
    y[1] = Cas12a protein
    """
    mRNA, Cas12a = y
    
    alpha_0 = params['alpha_0']
    alpha_max = params['alpha_max']
    K_d = params['K_d']
    n = params['n']
    gamma_m = params['gamma_m']
    beta_p = params['beta_p']
    gamma_p = params['gamma_p']
    
    # Promoter activation
    P_act = hill_activation(ara, alpha_0, alpha_max, K_d, n)
    
    # Differential equations
    dmRNA_dt = P_act - gamma_m * mRNA
    dCas12a_dt = beta_p * mRNA - gamma_p * Cas12a
    
    return [dmRNA_dt, dCas12a_dt]
