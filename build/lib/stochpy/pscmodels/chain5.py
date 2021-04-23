model = """R1: 
    S1 > S2 
    k1*S1
R2: 
    S2 > S1 
    k1*S2
R3: 
    S2 > S3 
    k1*S2
R4: 
    S3 > S2 
    k1*S3
R5: 
    S3 > S4 
    k1*S3
R6: 
    S4 > S3 
    k1*S4
R7: 
    S4 > S5 
    k1*S4
R8: 
    S5 > S4 
    k1*S5

# InitPar
k1 = 1 # noise

# InitVar
S1 = 1
S2 = 0
S3 = 0
S4 = 0
S5 = 0
"""
