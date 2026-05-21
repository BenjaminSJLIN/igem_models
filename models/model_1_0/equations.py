import numpy as np

def K_D (cas12a-crRNA, target_DNA, cas12a_A):
    """
    Hill equation for promoter activation by Arabinose.
    P_act([Ara]) = alpha_0 + alpha_max * ([Ara]^n / (K_d^n + [Ara]^n))
    """
    return cas12a-crRNA*targetDNA/cas12a

def cleavage_Rate(k_cat, cas12a_A, ssDNA, K_M):
    """
    System of ODEs for Cas12a expression.
    y[0] = mRNA
    y[1] = Cas12a protein
    """
    
