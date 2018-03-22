# -*- coding: utf-8 -*-

from z3 import *
from utils import buildONE, buildZERO, sub, extract1BV

## a refers to the bit-vector of the present activators ##
## of a given gene, size: |genes|                       ##

## r refers to the bit-vector of the present repressors ##
## of a given gene, size: |genes|                       ##

## Present regulators are regulators that CAN be made   ##
## active != active regulators in a given state         ##

## q refers to the bit-vector of the system state q     ##
## size: |genes|                                        ##

## Active regulators in state q are thus a & q, r & q   ##

#######################
#  SHORTCUTS          #
#######################

#_____________________#
#  SOLVER-SIDE        #
#_____________________#

## BIT-VECTOR OPERATIONS                                ##

## a refers to the bit-vector of the present activators ##
## of a given gene, size: |genes|                       ##

## Those bit-vectors encode conditions. They can be     ##
## casted to conditions using operation castBVinCond              ##
## where castBVinCond: b -> b == 1..1                             ##
## 1..1 is the full one vector                          ##

## Equality of all bit values                           ##
## <=> b1_i == b2_i, for all i                          ##
## <=> b1_i == b2_i == 1 OR b1_i == b2_i == 0, for all i##
nxor = lambda b1, b2 : b1&b2|~b1&~b2
xor = lambda b1, b2 : b1&~b2|b1&~b2
## Negation of a condition b encoded as a BV is 1..1    ##
## if ~b > 0, else ~b,                                  ##
## where ~ is the bit-wise negation                     ##
## SIGNED COMPARISON                                    ##
neg = lambda bv : If(UGT(~bv, buildZERO(bv.size())), buildONE(bv.size()), ~bv)
## Every bit of a = 0 <=> no present activators         ##
notInducible = lambda a : nxor(a, buildZERO(a.size())) 
## Every bit of r = 0 <=> no present repressors (same)  ##
notRepressible = notInducible

## CONDITION OPERATIONS                                 ##

## Caster bit-vector -> condition (BoolRef type)        ##
castBVinCond = lambda x : x == buildONE(x.size())

## Pre-computation of template functions knowing        ##
## regulators (activators and repressors) of a given    ##
## gene                                                 ##
def prepreCompute(a, r):
	n = a.size()
	zeroG = buildZERO(n)
	oneG = buildONE(n)
	notInducibleG = notInducible(a)
	notRepressibleG = notRepressible(r)
	negnotInducibleG = neg(notInducibleG)
	negnotRepressibleG = neg(notRepressibleG)
	## a != 0 and all (present) activators are active       ##
	## i.e. q & a == a                                      ##œ		
	allActivatorsG = lambda q : negnotInducibleG & nxor(q & a, a)
	## r != 0 and all (present) repressors are active       ##
	## i.e. q & r == r (same)                               ##
	allRepressorsG = lambda q : negnotRepressibleG & nxor(q & r, r)
	## a == 0 or all (present) activators are NOT active    ##
	## i.e. q & a == 0                                      ##
	noActivatorsG = lambda q : nxor(q & a, zeroG)
	## a == 0 or all (present) activators are NOT active    ##
	## i.e. q & a == 0 (same)                               ##
	noRepressorsG = lambda q : nxor(q & r, zeroG)
	indRegG = neg(negnotInducibleG & notRepressibleG)
	repRegG = negnotRepressibleG & notInducibleG
	return([allActivatorsG, allRepressorsG, noActivatorsG, noRepressorsG, indRegG, repRegG])

## Pre-computation of template functions knowing        ##
## regulators (activators and repressors) of a given    ##
## gene, AND global system state                        ##
def preCompute(q, prepreComputationG, test=False):
	[allActivatorsG, allRepressorsG, noActivatorsG, noRepressorsG, indRegG, repRegG] = prepreComputationG
	allActivatorsGQ = allActivatorsG(q)
	allRepressorsGQ = allRepressorsG(q)
	noActivatorsGQ = noActivatorsG(q)
	noRepressorsGQ = noRepressorsG(q)
	negallActivatorsGQ = neg(allActivatorsGQ)
	negallRepressorsGQ = neg(allRepressorsGQ)
	negnoActivatorsGQ = neg(noActivatorsGQ)
	negnoRepressorsGQ = neg(noRepressorsGQ)
	inducibleRegulationGQ = indRegG | negnoActivatorsGQ
	repressibleRegulationGQ = repRegG & noRepressorsGQ
	if (test):
		getRegFctGQ = lambda x : (x & inducibleRegulationGQ) | repressibleRegulationGQ
	else:
		getRegFctGQ = lambda x : castBVinCond((x & inducibleRegulationGQ) | repressibleRegulationGQ)
	return([allActivatorsGQ, allRepressorsGQ, noActivatorsGQ, noRepressorsGQ, negallActivatorsGQ, negallRepressorsGQ, negnoActivatorsGQ, negnoRepressorsGQ, getRegFctGQ])

## allActivatorsGQ, allRepressorsGQ, noActivatorsGQ, noRepressorsGQ, ##
## negallActivatorsGQ, negallRepressorsGQ, negnoActivatorsGQ,        ##
## negnoRepressorsGQ, getRegFctGQ                                    ##
## arg = aA, aR, nA, nR, naA, naR, nnA, nnR, g                       ##
diCompute = {  0: lambda aA, aR, nA, nR, naA, naR, nnA, nnR, g : g(aA & nR),
	1: lambda aA, aR, nA, nR, naA, naR, nnA, nnR, g : g(nnA & nR),
	2: lambda aA, aR, nA, nR, naA, naR, nnA, nnR, g : g(aA & naR),
	3 : lambda aA, aR, nA, nR, naA, naR, nnA, nnR, g : g((nR & nnA) | (naR & aA)),
	4: lambda aA, aR, nA, nR, naA, naR, nnA, nnR, g : g(aA),
	5: lambda aA, aR, nA, nR, naA, naR, nnA, nnR, g : g(aA | (nR & nnA)),
	6: lambda aA, aR, nA, nR, naA, naR, nnA, nnR, g : g(nnA & naR),
	7: lambda aA, aR, nA, nR, naA, naR, nnA, nnR, g : g((nnA & naR) | aA),
	8: lambda aA, aR, nA, nR, naA, naR, nnA, nnR, g : g(nnA),
	9: lambda aA, aR, nA, nR, naA, naR, nnA, nnR, g : g(nR),
	10: lambda aA, aR, nA, nR, naA, naR, nnA, nnR, g : g(nR | (naR & aA)),
	11: lambda aA, aR, nA, nR, naA, naR, nnA, nnR, g : g(nR | (nnA & naR)),
	12: lambda aA, aR, nA, nR, naA, naR, nnA, nnR, g : g(naR),
	13: lambda aA, aR, nA, nR, naA, naR, nnA, nnR, g : g(nR | aA),
	14: lambda aA, aR, nA, nR, naA, naR, nnA, nnR, g : g((nR | aA) | (naR & nnA)),
	15: lambda aA, aR, nA, nR, naA, naR, nnA, nnR, g : g(naR | aA),
	16: lambda aA, aR, nA, nR, naA, naR, nnA, nnR, g : g(nR | nnA),
	17: lambda aA, aR, nA, nR, naA, naR, nnA, nnR, g : g(naR | nnA)
	}

def rule18(q, a, r, c):
	nq = q.size()
	zeroQ = buildZERO(nq)
	oneQ = buildONE(nq)
	oneQ1 = buildONE(nq+1)
	sr = sub(q & r)
	sa = sub(q & a)
	return(castBVinCond(If(UGT(sa, sr), oneQ, If(nxor(Concat(extract1BV(q, c), sa), Concat(buildONE(1), sr)) == oneQ1, oneQ, zeroQ))))

def rule19(q, a, r, c):
	return(UGT(sub(q & a), sub(q & r)))

#_____________________#
#  SIMPLIFIER-SIDE    #
#_____________________#

## GRF GENERATORS                                       ##

def writeKO(pgN, x):
	return("(NOT (" + pgN + ") AND " + x + ")")

def writeFE(pgN, x):
	return("(" + pgN + " OR " + x + ")")

def writePerturbation(kog, feg, x, kogN, fegN):
	if (kog and feg):
		return(writeFE(fegN, writeKO(kogN, x)))
	elif (kog):
		return(writeKO(kogN, x))
	elif (feg):
		return(writeFE(fegN, x))
	else:
		return(x)

## RE:IN style                               ##
#def writePerturbation(kog, feg, x, kogN, fegN):
#	if (kog and feg):
#		return("((NOT (" + kogN + ") OR " + fegN + ") AND " + x + ")")
#	elif (kog):
#		return("(NOT (" + kogN + ") AND " + x + ")")
# 	elif (feg):
#		return("(" + fegN + " OR " + x + ")")
#	else:
#		return(x)

def writePerturbedGene(pg):
	return(pg)

def allActivatorsGRF(a):
	## allActivators needs the existence of at least ONE activator ##
	## and that any existent activator is active                   ##
	if (not a):
		return("False")
	return("(" + reduce(lambda x, y : x + " AND " + y, a) + ")")

def noActivatorsGRF(a):
	## There is no activator indeed                                ##
	if (not a):
		return("True")
	return("((NOT " + reduce(lambda x, y : x + ") AND (NOT " + y, a) + "))")

allRepressorsGRF = allActivatorsGRF
noRepressorsGRF = noActivatorsGRF

def wGRF(a, r, x):
	## means a == None ##
	notInducible = not a
	## means r == None ##
	notRepressible = not r
	inducibleRegulation = not (not notInducible and notRepressible)
	repressibleRegulation = not notRepressible and notInducible 
	if (inducibleRegulation and repressibleRegulation):
		return("(" + x + " OR " + noRepressorsGRF(r) + ")")
	elif (inducibleRegulation):
		return(x)
	elif (repressibleRegulation):
		return("((" + x  + " AND " + "(NOT " + noActivatorsGRF(a) + ")" + ")" + " OR " + noRepressorsGRF(r) + ")")
	else:
		return("(" + x + " AND " + "(NOT " + noActivatorsGRF(a) + ")" + ")")

def listActivators(a):
	if (not a):
		return("0")
	else:
		return("[" + reduce(lambda x,y : x + ", " + y, a) + "]")

listRepressors = listActivators	

diGRF = lambda i : {
       	 0: lambda a, r, c : wGRF(a, r, "(" + allActivatorsGRF(a) + " AND " + noRepressorsGRF(r) + ")"),
      	 1: lambda a, r, c : wGRF(a, r, "((NOT " + noActivatorsGRF(a) + ") AND " + noRepressorsGRF(r) + ")"),
         2: lambda a, r, c : wGRF(a, r, "(" + allActivatorsGRF(a) + " AND " + "(NOT " + allRepressorsGRF(r) + "))"),
         3: lambda a, r, c : wGRF(a, r, "((" + noRepressorsGRF(r) + " AND " + "(NOT " + noActivatorsGRF(a) + "))" + " OR " 
		+ "((NOT " + allRepressorsGRF(r) + ") AND " + allActivatorsGRF(a) + "))"),
         4: lambda a, r, c : wGRF(a, r, allActivatorsGRF(a)),
	 5: lambda a, r, c : wGRF(a, r, "(" + allActivatorsGRF(a) + " OR " + "(" + noRepressorsGRF(r) + " AND " 
		+ "(NOT " + noActivatorsGRF(a) + ")))"),
	 6: lambda a, r, c : wGRF(a, r, "((NOT " + noActivatorsGRF(a) + ")" + " AND " + "(NOT " + allRepressorsGRF(r) + "))"),
	 7: lambda a, r, c : wGRF(a, r, "(((NOT " + noActivatorsGRF(a) + ")" + " AND " + "(NOT " + allRepressorsGRF(r) + "))" 
		+ " OR " + allActivatorsGRF(a) + ")"),
	 8: lambda a, r, c : wGRF(a, r, "(NOT " + noActivatorsGRF(a) + ")"),
	 9: lambda a, r, c : wGRF(a, r, noRepressorsGRF(r)),
	 10: lambda a, r, c : wGRF(a, r, "(" + noRepressorsGRF(r) + " OR " + "((NOT " + allRepressorsGRF(r) + ")" 
		+ " AND " + allActivatorsGRF(a) + "))"),
	 11: lambda a, r, c : wGRF(a, r, "(" + noRepressorsGRF(r) + " OR " + "((NOT " + noActivatorsGRF(a) + ")" 
		+ " AND (NOT " + allRepressorsGRF(r) + "))" + ")"),
	 12: lambda a, r, c : wGRF(a, r, "(NOT " + allRepressorsGRF(r) + ")"),
	 13: lambda a, r, c : wGRF(a, r, "(" + noRepressorsGRF(r) + " OR " + allActivatorsGRF(a) + ")"),
	 14: lambda a, r, c : wGRF(a, r, "((" + noRepressorsGRF(r) + " OR " + allActivatorsGRF(a) + ")" + " OR " 
		+ "((NOT " + allRepressorsGRF(r) + ")" + " AND (NOT " + noActivatorsGRF(a) + ")))"),
	 15: lambda a, r, c : wGRF(a, r, "((NOT " + allRepressorsGRF(r) + ")" + " OR " + allActivatorsGRF(a) + ")"),
	 16: lambda a, r, c : wGRF(a, r, "(" + noRepressorsGRF(r) + " OR (NOT " + noActivatorsGRF(a) + "))"),
	 17: lambda a, r, c : wGRF(a, r, "((NOT " + allRepressorsGRF(r) + ")" + " OR (NOT " + noActivatorsGRF(a) + "))"),
	 ## Threshold rules   ##
	 18: lambda a, r, c : "(" + "(" + listActivators(a) + " > " + listRepressors(r) 
		+ ") OR ((" + listActivators(a) + " = " + listRepressors(r) + ") AND " + c + "))",
	 19: lambda a, r, c : "(" + listActivators(a) + " > " + listRepressors(r) + ")"
}

