model = """
# PySCeS test input file
# Stochastic Simulation Algorithm input format
# BurstModel

R1:
    ONstate > OFFstate
    koff*ONstate

R2: 
    OFFstate > ONstate
    kon*OFFstate

R3:
    ONstate > mRNA + ONstate
    ksyn*ONstate

R4:
    mRNA > $pool
    kdeg*mRNA

# InitPar
kon  = 0.05
koff = 0.05
kdeg = 2.5
ksyn = 80

# InitVar
ONstate = 0
OFFstate = 1
mRNA = 0
"""

