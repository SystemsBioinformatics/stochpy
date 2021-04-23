model = """# signalling pathway incl. feedback
# R and S on the same operon
# Delay can be applied such that mRNA_S is produced after mRNA_R via mRNA_prelim
# RP dimerises before activating genes
# Rate constants set for bistability in Cell Division module
# Set DNA and DNAa as 'Non-dividing' species in Stochpy


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

#Basal transcription/translation and degradation

#Reaction names correspond to relevant parameter names in Table S2 of the Supplementary Information
#Recation names are further appended by letters if multiple reactions are associated with the same parameter

#Note: mRNA_prelim is a 'trick' species to model production of mRNA_S strictly after delay following mRNA_R production
Rp1a:
	$pool > mRNA_R + mRNA_Sprelim 
	kr1*(DNA+DNAa)
	
Rp1b: 
	mRNA_Sprelim > mRNA_S #approximately instantaneous reaction (so mRNA_prelim doesn't really 'exist')
	10000*mRNA_Sprelim	
Rp3a:
	$pool > R
	kr2*mRNA_R	
Rd1a:
	mRNA_R > $pool
	kr3*mRNA_R	
Rp3b:
	$pool > S
	ks2*mRNA_S
	
Rd1b:
	mRNA_S > $pool
	ks3*mRNA_S	

	
#Feedback reactions
	
RDimerf:
	RP + RP > RPDimer
	kDimerf*RP*(RP-1)
	
RDimerb:
	RPDimer > RP + RP
	kDimerb*RPDimer
	
RActf:
	RPDimer + DNA > DNAa
	kActf*RPDimer*DNA
	
RActb:
	DNAa > RPDimer + DNA
	kActb*DNAa
	
Rp2:
	$pool > mRNA_R + mRNA_Sprelim
	kr6*DNAa
	
	
#Fixed Species

#Variable species:
DNA = 1
DNAa = 0
S = 0
SL = 0
SP = 0
SLP = 100
SLPR = 0
SPR = 0
RPS = 0
RP = 80
RPDimer = 0
R = 0
mRNA_R = 5
mRNA_S = 5
mRNA_Sprelim = 0

# Parameters
L=200
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

kr1 = 0.003
kr2 = 0.002
kr3 = 0.002
kr6 = 0.1

ks2 = 0.002
ks3 = 0.002

kDimerf = 0.001
kDimerb = 1
kActf = 0.001
kActb = 0.01



"""
