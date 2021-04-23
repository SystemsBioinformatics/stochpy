model = """# signalling pathway with response times to stepwise L changes, 'up-down' staircase

FIX: L

#Reactions
R1f:
	S > SL 
	k1f*S*L
	
R2:
	SL > SLP 
	k2*SL
	
R3f:
	SLP > SP
	k3f*SLP
	
R4f:
	SLP + R > SLPR
	k4f*SLP*R 
	
R5f:
	SP + R > SPR
	k5f*SP*R
	
R6f:
	SLPR > SL + RP
	k6f*SLPR
	
R7f:
	SPR > S + RP
	k7f*SPR
	
R8f:
	RP + S > RPS
	k8f*RP*S
	
R9:
	RPS > R + S
	k9*RPS
	
R10f:
	SLPR > SPR
	k10f*SLPR
	
R11:
	RP > R
	k11*RP
	
R12:
	SP > S
	k12*SP
	
R13:
	SLP > SL
	k13*SLP
	
R1b:
	SL > S
	k1b*SL

R3b:
	SP > SLP
	k3b*SP*L

R4b:
	SLPR > SLP + R
	k4b*SLPR

R5b:
	SPR > SP + R
	k5b*SPR
	
R6b:
	SL + RP > SLPR
	k6b*SL*RP
	
R7b:
	S + RP > SPR
	k7b*S*RP
	
R8b:
	RPS > S + RP
	k8b*RPS
	
R10b:
	SPR > SLPR
	k10b*SPR*L
	
#Event Definitions

Event: jump1, operator.ge(_TIME_,200), 0.0 
{
L = 20
}

Event: jump2, operator.ge(_TIME_,400), 0.0 
{
L = 40
}

Event: jump3, operator.ge(_TIME_,600), 0.0 
{
L = 60
}

Event: jump4, operator.ge(_TIME_,800), 0.0 
{
L = 80
}

Event: jump5, operator.ge(_TIME_,1000), 0.0 
{
L = 100
}

Event: jump6, operator.ge(_TIME_,1400), 0.0 
{
L = 80
}

Event: jump7, operator.ge(_TIME_,1600), 0.0 
{
L = 60
}

Event: jump8, operator.ge(_TIME_,1800), 0.0 
{
L = 40
}

Event: jump9, operator.ge(_TIME_,2000), 0.0 
{
L = 20
}

Event: jump10, operator.ge(_TIME_,2200), 0.0 
{
L = 1
}

#Fixed species

#Variable species:
S = 100
SL = 0
SP = 0
SLP = 0
SLPR = 0
SPR = 0
L = 1
RPS = 0
RP = 0
R = 100

# Parameters
k1f = 0.01
k1b = 0.059
k2 = 0.1
k3f = 0.001
k3b =0.001
k4f = 0.0238
k4b = 0.001
k5f = 0.0238
k5b = 0.001
k6f = 19.8
k6b = 0.001
k7f = 19.8
k7b = 0.001
k8f = 0.0238
k8b = 0.2
k9 = 0.5
k10f = 0.001
k10b = 0.001
k11 = 0.0001
k12 = 0.0001
k13 = 0.0001
"""
