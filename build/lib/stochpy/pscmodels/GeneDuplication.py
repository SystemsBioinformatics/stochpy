model = """

# Reactions
R1:
    G1 >  G1 + mRNA1
    Ksyn*G1

R2:
    mRNA1 > $pool
    Kdeg*mRNA1
R3:
    G2 >  G2 + mRNA2
    Ksyn*G2
R4:
    mRNA2 > $pool
    Kdeg*mRNA2 
 
# Fixed species
 
# Variable species
mRNA1 = 50.0
G1 = 1
mRNA2 = 50.0
G2 = 1
 
# Parameters
Ksyn = 10
Kdeg = 0.2
"""
