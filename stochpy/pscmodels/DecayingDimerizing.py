model = """
# PySCeS test input file
# Stochastic Simulation Algorithm input format
# Decay Dimerazing model

# S1  --> 0    
# 2S1 --> S2
# S2  --> 2S1
# S2  --> S3

R1:
    S1 > $pool
    S1*k1

R2: 
    S1 + S1 > S2
    0.5*k2*S1*(S1-1)

R3:
    S2 > S1 + S1
    k3*S2

R4:
    S2 > S3
    k4*S2

#InitPar
k1 = 1.0
k2 = 0.002
k3 = 0.5
k4 = 0.04

#InitVar
S1 = 100000
S2 = 0.0
S3 = 0.0
"""
