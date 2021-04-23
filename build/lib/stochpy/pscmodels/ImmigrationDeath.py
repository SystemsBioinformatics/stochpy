model = """# Stochastic Simulation Algorithm input file
# --> mRNA --> 

# Reactions
R1:
    $pool > mRNA
    Ksyn

R2:
    mRNA > $pool
    Kdeg*mRNA
 
# Fixed species
 
# Variable species
mRNA = 50.0
 
# Parameters
Ksyn = 10
Kdeg = 0.2
"""
