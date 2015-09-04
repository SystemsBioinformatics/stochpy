model = """# Stochastic Simulation Algorithm input format

R1:
    $pool > TF
    kTFsyn
R2:
    TF > $pool
    kTFdeg*TF
R3:
    TFactive > $pool
    kTFdeg*TFactive
R4:
    TF > TFactive
    kActivate*TF  
R5: 
    TFactive > TF
    kInactivate*TFactive
R6:
    TFactive > mRNA + TFactive
    kmRNAsyn*(TFactive/(TFactive+kX))
R7:
    mRNA > $pool
    kmRNAdeg*mRNA
R8:
    mRNA > Protein + mRNA
    kProteinsyn*mRNA
R9:
    Protein > $pool
    kProteindeg*Protein
    
# InitPar
kTFsyn = 200
kTFdeg = 20
kActivate = 2000
kInactivate = 200
kX = 5
kmRNAsyn = 240
kmRNAdeg = 20
kProteinsyn = 400
kProteindeg = 2

# InitVar
TF = 2
TFactive = 10
mRNA = 10
Protein = 220
"""
