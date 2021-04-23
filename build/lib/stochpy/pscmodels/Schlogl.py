model = """
# Schlogl model (Schlogl 1972, Chemical reaction models for nonequilibrium phase transitions)

FIX: A B

# Reactions
R1:
    A + {2} X > {3} X
    1/2*c1*A*X*(X-1)

R2:
    {3} X > A + {2} X
    1/6*c2*X*(X-1)*(X-2)
R3: 
    B > X
    c3 * B
R4:
    X > B
    c4*X
 
# Fixed species
A = 100000
B = 200000 
 
# Variable species
X = 250

c1 = 3*10**-7
c2 = 10**-4
c3 = 10**-3
c4 = 3.5
"""
