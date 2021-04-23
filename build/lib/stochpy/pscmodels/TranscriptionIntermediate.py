model = """# Stochastic Simulation Algorithm input format

R1:
    Polymerase          > PolymeraseMoving 
    kini*Polymerase
    
R2:
    PolymeraseMoving    > mRNA  + Polymerase
    0.1*PolymeraseMoving # SingleMolecule non-exponential reaction.
    
R3:
    mRNA > $pool
    kdeg*mRNA

#Species
mRNA = 0
Polymerase = 10
PolymeraseMoving = 0

# Rates
kini = 1/60 #1/min
kdeg = 1/600 #1/10 min
"""
