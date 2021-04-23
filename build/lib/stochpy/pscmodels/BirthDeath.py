model = """# Stochastic Simulation Algorithm input file
# --> mRNA --> 

# Reactions
R1:
    mRNA > {2} mRNA
    Ksyn*mRNA

R2:
    mRNA > $pool
    Kdeg*mRNA
 
# Fixed species
 
# Variable species
mRNA = 100
 
# Parameters
Ksyn = 2.9
Kdeg = 3
"""
