model = """# NonConsuming vs Consuming delayed reactions

# Reactions
R1:
    polymerase > mRNA + polymerase
    Ksyn*polymerase

R2:
    mRNA > $pool
    Kdeg*mRNA

# Variable species
polymerase  = 10
mRNA = 0

# Parameters
Ksyn = 0.5
Kdeg = 0.1
"""
