# -*- coding: utf-8 -*-

########################################################
## UTILS FUNCTIONS TO PROCESS Z3 BIT-VECTORS          ##
########################################################

from z3 import *

#________#
# Debug  #
#________#

def testSolver(s, C):
    print(s.check())
    try:
    	print(s.model())
    	print(C)
    except:
	print("No model")
    return(None)

def checkSolver(s):
	print(s.sexpr())
	return(s)

def printPretty(x):
	ls = ["\t" + str(a[0]) + " = " + str(a[1]) for a in x]
	for a in ls:
		print(a)

def verboseIt(msg, verbose):
	if (verbose):
		print("MSG: " + msg)
	return(None)

def f(x):
	try:
		return(simplify(x))
	except:
		return(None)

def addArgValues(s, argnames, args):
	## length(argnames)==length(args) ##
	for i in range(len(args)):
		s += " " + argnames[i] + "=" + str(args[i]) + ", "
	return(s)

def noBV(adj="", argnames=None, args=None):
	s = "ERROR: A " + adj + " Bit-Vector could not be created: "
	return(addArgValues(s, argnames, args))

def warnBV(warning="", argnames=None, args=None):
	s = "WARNING: " + warning + ": "
	return(addArgValues(s, argnames, args))

def regularError(msg=""):
	return("ERROR: " + msg + ".")

def strList2Str(ls):
	res = ""
	for s in ls:
		res += s
	return(res)

#____________________________#
# Build generic Bit-Vectors  #
#____________________________#

#' Construct Bit-Vector encoding 2^i of given size
#'
#' @param i log2 of value 
#'          (integer)
#' @param s size of the vector
#'          (positive integer)
#' @return b Z3 Bit-Vector encoding 2^i of size s 
def bv(i, s):
	return(BitVecVal(2**i, s))

#' Construct Bit-Vector 0,0,0,0,...,0
#' of size l
#'
#' @param l size of the vector
#'          (positive integer)
#' @return b zero Z3 Bit-Vector of size l 
def buildZERO(l):
	if (not l):
		print(noBV("ZERO", ["l"], [l]))
		return(None)
	return(BitVecVal(0, l))

#' Construct Bit-Vector 1,1,1,1,...,1
#' of size l
#'
#' @param l size of the vector
#'          (integer)
#' @return b full-one Z3 Bit-Vector of size l 
def buildONE(l):
	if (not l):
		print(noBV("ONE", ["l"], [l]))
		return(None)
	return(BitVecVal(2**l-1, l))

#' Construct Bit-Vector with the 
#' coordinate i removed: 
#' b1, ..., bi-1, bi+1, ..., bn
#'
#' @param b Z3 Bit-Vector
#' @param i index of the coordinate to remove
#'          (integer)
#' @return res Z3 Bit-Vector minus the i^th coordinate
def trim1BV(b, i):
	n = b.size()
	if (i == None or i >= n or i < 0 or n==0):
		print(noBV("TRIMMED", ["i", "n"], [i, n]))
		return(None)
	if (i==0 and n==1):
		print(noBV("EMPTY", ["i", "n"], [i, n]))
		return(None)
	if (i==0):
		return(Extract(n-1, 1, b))
	if (i==n-1):
		return(Extract(n-2, 0, b))
	else:
		return(Concat(Extract(n-1, i+1, b), Extract(i-1, 0, b)))

#' Construct Bit-Vector with only the 
#' coordinate i: bi
#'
#' @param b Z3 Bit-Vector
#' @param i index of the coordinate
#'          (integer)
#' @return res Z3 Bit-Vector of size 1 with i^th coordinate
def extract1BV(b, i):
	n = b.size()
	if (i == None or i >= n or i < 0 or n==0):
		print(noBV("1-SIZED", ["i", "n", "b"], [i, n, b]))
		return(None)
	return(Extract(i, i, b))

#' Construct constant vector using the list 
#' of non-negative bit indices
#'
#' @param ls list of non-negative bit indices
#'           (integer list)
#' @param size size of the final Bit-Vector
#' @return res Z3 Bit-Vector of size size
def intList2BV(ls, size):
	if (not(ls)):
		return(buildZERO(size))
	if (any([l >= size for l in ls]) or any([l==None for l in ls])):
		print(warnBV("Bit-Vector will be trimmed", ["l", "size"], [l, size]))
	return(BitVecVal(sum([2**i for i in ls]), size))

default=-1
false = buildZERO(1)
true = buildONE(1)

#____________________________#
# Compare Bit-Vectors        #
#____________________________#

#' Test q[i] == v, where q is a Bit-Vector
#'
#' @param q Z3 Bit-Vector
#' @param i index of the bit to be tested
#'          (integer)
#' @param v value 
#'          (Z3 bit-vector of size 1)
#' @return res condition to add to solver 
def testValue(q, i, v):
	if (v==None):
		print(noBV("NULL BIT-VECTOR", ["q", "i"], [q, i]))
		return(None)
	if (v.size()!=1):
		print(noBV("WRONG", ["|v|", "q", "i"], [v.size(), q, i]))
		return(None)		
	qq = extract1BV(q, i)
	if (qq == None):
		print(regularError("NULL extracted Bit-Vector"))
		return(None)
	return(qq==v)

#' Test q[idx!=i] == q1[idx!=i]
#' where q, q1 are Bit-Vectors
#'
#' @param q Z3 Bit-Vector
#' @param q1 Z3 Bit-Vector
#' @param i index of the bit to be ruled out
#'          (integer)
#' @return res condition to add to solver 
def equalExcept1BV(q, q1, i):
	qq1 = trim1BV(q1, i)
	qq = trim1BV(q, i)
	if (qq1==None or qq==None):
		print(warnBV("One of the vectors to compare is NULL", ["qq", "qq1"], [qq, qq1]))
		return(False)
	return(qq1==qq)

#' Test q[c] == q1[c1], where q, q1 are Bit-Vectors
#'
#' @param q Z3 Bit-Vector
#' @param q1 Z3 Bit-Vector
#' @param c index of the bit to be tested in q
#'          (integer)
#' @param c1 index of the bit to be tested in q1
#'          (integer)
#' @return cond condition to add to solver  
def testValueVec(q, q1, c, c1=None):
	if (c1 == None):
		c1 = c
	qq = extract1BV(q, c)
	qq1 = extract1BV(q1, c1)
	if (qq == None or qq1 == None):
		print(warnBV("One of the vectors to compare is NULL", ["qq", "qq1"], [qq, qq1]))
		return(False)
	return(qq==qq1)

#' Returns condition b1 == b on (coordinate) non-default values in dict
#' else if default value then b1 == 1
#' else b1 == 0
#'
#' @param s Z3 solver
#' @param b Z3 Bit-Vector
#' @param b1 Z3 Bit-Vector
#' @param di dictionary of indices of b1 (keys) and b (values)
#' @return s updated solver 
def evalVec(s, b, b1, di, default):
	n = b1.size()
	if (n<1):
		print(regularError("NULL Bit-Vector"))
		return(None)
	for i in range(n):
		if (i in di.keys()):
			if (di.get(i) != default):
				s.add(testValueVec(b1, b, i, c1=di.get(i)))
			else:
				s.add(testValue(b1, i, true))
		else:
			s.add(testValue(b1, i, false))
	return(s)

#________#
# Tools  #
#________#

## b != e !                             ##
## b before e in x                      ##
getIt = lambda x, b, e : reduce(lambda a, b : a + b, x.split(b)[1].split(e)[0].split(" "))

rev = lambda ls : list(reversed(ls))

## Count number of 1-bits in Bit-Vector ##
def sub(q):
	n = q.size()
	bits = [extract1BV(q, i) for i in range(n)]
	bvs = [Concat(buildZERO(n-1), b) for b in bits]
	return(reduce(lambda a, c: a+c, bvs))

## Full: get list of binary integers    ##
## corresponding to input Bit-Vector r  ##
## Full iff. a size is given            ##
## Not Full: get integer corresponding  ##
## to input Bit-Vector r                ##
def getBinaryDec(r, size=None):
	n = r.size()
	i = int(str(simplify(BV2Int(r))))
        if (size):
		options = '0' + str(int(size)) + 'b'
		return(format(i, options))
	else:
		return(i)

def filterNone(ls):
	return(filter(lambda x : x != None, ls))

def ifthenelse(ifc, thenc, elsec=None):
	return(thenc if (ifc) else elsec)

def ifthenelseFct(ifc, thenc, setc, elsec=None):
	return(filterNone([thenc(i) if (ifc(i)) else elsec for i in range(len(setc))]))

def ifthenelseSet(ifc, thenc, setc, elsec=None):
	return(filterNone([thenc[i] if (ifc(i)) else elsec for i in range(len(setc))]))

def concat(l, sep=""):
	return(reduce(lambda x,y : x + sep + y, l))

def printStates(trajectories, i, C):
	chunksC = [C[j:j + 10] for j in xrange(0, len(C), 10)]
	## Print trajectory number
	print(trajectories[i][0] + "\n")
	path = trajectories[i][1]
	## get step name
	getStep = lambda t : path[t][0][2:].split("_")[0].replace("_", " ")
	gene_pos = 0
	for ck in chunksC:
		## Header for gene names
		print(" "*(len(path[-1][0][2:].split("_")[0])+2) + concat(ck, sep=" "))
		for t in range(len(path)):
			nn = len(path[-1][0])-len(path[t][0])+1
			print(getStep(t) + " " + ifthenelse(nn > 0, " "*nn, "") + concat([path[t][1][gene_pos:(gene_pos+len(ck))][l]+(" "*len(ck[l])) for l in range(len(ck))]))
		gene_pos += len(ck)
