# -*- coding: utf-8 -*-

from global_paths import *

##########################
## Read full REIN file  ##
##########################

#________________#
#   Tools        #
#________________#

#' Get directive markers in RE:IN file
#'
#' @param f Python file object associated with RE:IN file
#' @return res list of directives (with conversion to integer if possible)
def getDirective(f):
	x = f.readline().split(";")[0].split(" ")[2]
	try:
		return(int(x))
	except:
		return(x)

## Convert a well-formatted string into a list of integers                ##
getRegularList = lambda r : map(int, r.split(","))
## Convert a well-formatted string into a list of integers                ##
getRange = lambda l : range(int(l[0]), int(l[1])+1)
## Convert a well-formatted string into string x string x "+" or "-"      ##
getInteraction = lambda x : x[:2] + ["+" if (x[2]=="positive") else "-"]
## Implement custom grep                                                  ##
grep = lambda x, pattern : len(x.split(pattern)) > 1
## Sanitize strings using starting and ending characters                  ##
getIt = lambda x, b, e : reduce(lambda a, b : a + b, x.split(b)[1].split(e)[0].split(" "))
## Append an integer n to each list element of a list x                   ##
updateStep = lambda x, n : [[n]+v[1:] for v in x]
## Sanitize strings using pattern                                         ##
splitNclean = lambda x, pattern : filter(lambda x : x != "", x.split(pattern))

#________________#
#   Model        #
#________________#

#' Parse model from RE:IN file
#'
#' @param f Python file object associated with RE:IN file
#' @param verbose logical for printing messages
#' @return res list that contains the description of the abstract model
#' [C, CRM, length, Idef, Iopt, R, typeT, solmax, uniqueness, limreg, P]
#' - C is the set of genes/nodes
#' - CRM is the list of associated regulated genes of length #nodes (None if the 
#' considered node is not a regulatory module)
#' - length maximum length of the model
#' - Idef set of definite interactions
#' - Iopt set of optional interactions
#' - R list of allowed regulation conditions for each gene/node
#' - typeT transition type in the model ("asynchronous" or "synchronous")
#' - solmax maximum number of solutions returned by the solver
#' - uniqueness string for solution uniqueness condition ("interaction", "full" or "paths")
#' - limreg (not used) limitation on the set of regulation conditions used
#' - P set of known perturbations for each node
def getModel(f, verbose):
	## Parameters      ##
	typeT = getDirective(f) + "hronous"
	length, uniqueness, solmax, limreg = [getDirective(f) for i in range(4)]
	x = [v.split("[") for v in f.readline().split("; ")][:-1]
	## Genes           ##
	C = [r[0] for r in x]
	P = [list(r[1].split("]")[0]) for r in x]
	CRM = [getIt(r[1], "]{", "}") for r in x]
	RR = [getIt(r[1], "}(", ")") for r in x]
	## Regulation      ##
	R = []
	for r in RR:
		l = r.split("..")
		if (len(l) == 1):
			R.append(getRegularList(r))
		else:
			ll = r.split(",")
			if (len(ll) == 1):
				R.append(getRange(l))
			else:
				rls = []
				for rl in ll:
					rll = rl.split("..")
					if (len(rll) == 1):
						rls.append(int(rl))
					else:
						rls += getRange(rll)
				R.append(rls)
	## Interactions    ##
	x = f.readline()
	Idef, Iopt = [], []
	while (x):
		x = x.split(";")[0].split("\t")
		if (len(x) == 4):
			Iopt.append(getInteraction(x))
		elif (len(x) == 3):
			Idef.append(getInteraction(x))
		else:
			print("ERROR PARSING: length(interaction) = " + str(len(x)) + " != 3, 4.")
			return(None)
		x = f.readline()
	## Debugging            ##
	if (verbose):
		print("---- MODEL       ----")
		for k, v in {"C": C, "length": length, "Idef": Idef, "Iopt": Iopt, "R": R, 
				"typeT": typeT, "solmax": solmax, 
				"uniqueness": uniqueness, "limreg": limreg, "P": P}.items():
			print(k + " = " + str(v))
	return([C, CRM, length, Idef, Iopt, R, typeT, solmax, uniqueness, limreg, P])

#________________#
# Experiments    #
#________________#

#' Parse experimental conditions
#' 
#' @param x string associated with RE:IN experiments file
#' @return res list [Fixpoint, summary, condexp]
#' - Fixpoint is a list that contains the step of a fix point, 
#' and the associated experimental condition
#' - summary is a list describing the experimental conditions at each step
#' in each experiment
#' - condexp is a dictionary with keys: condition names, values: list of (node name, value)
#' in the associated condition
def readConditions(x):
	Fixpoint = []
	summary = []
	condexp = dict()
	for line in x:
		line = splitNclean(line, "\r\n")
		m = len(line)
		i = 0
		while (i < m):
			sub = line[i]
			## Experiment definition ##
			if (grep(sub, "|=")):
				sb = sub.split("|=")
				n = int(getIt(sb[0], "[", "]"))
				exp = getIt(sb[0], "#", "[")
				expr = getIt(sb[1], "$", " ")
				summary.append([n, exp, expr])
			## Condition definition  ##
			if (grep(sub, ":=")):
				if (i < m-1):
					sb = reduce(lambda x, y : x + y, line[i:])
				else:
					sb = sub
				li = getIt(sb, "$", ";").split(":=")
				key = li[0]
				value = getIt(li[1], "{", "}").split("and")			
				value = [e.split("=") for e in value]
				value = [[0]+[e[0], int(e[1])] for e in value]
				condexp.setdefault(key, value)
			## Fixpoint definition   ##
			if (grep(sub, "fixpoint")):
				n = int(getIt(sub, "[", "]"))
				exp = getIt(sb[0], "#", "[")
				Fixpoint.append([n, exp])
			## Other: comments       ##		
			i += 1
	return([Fixpoint, summary, condexp])

#' Parse RE:IN experiments file
#' 
#' @param f Python file object associated with RE:IN experiments file
#' @param C set of genes/nodes
#' @param length maximum length of experiment
#' @param verbose logical for printing messages
#' @return res list that describes the experiments
#' [E, KO, FE, Fixpoint]
#' - E is a set describing each "regular" experiment
#' - KO is a set describing KO perturbations
#' - FE is a set describing FE perturbations
#' - Fixpoint is a set describing fix point/steady state constraints
def getExperiments(f, C, length, verbose):
	x = splitNclean(f.read(), ";")
	## Remove last empty line ##
	y = [splitNclean(v, "\r\n") for v in x][:-1]
	E = []
	KO = []
	FE = []
	[Fixpoint, summary, condexp] = readConditions(x)
	cond = None
	stack = []
	e = None
	ee = []
	for [n, exp, expr] in summary:
		## Transition to new exp.##
		if (exp != e):
			## Empty stack          ##
			## If no condition seen ##
			if (stack):
				ee += reduce(lambda x, y : x + y, [ex for ex in stack])
				stack = []
			## Push previous exp. ##
			if (ee and e):
				E.append([e, ee])
			## Initialize new exp.##
			cond = None
			stack = []
			e = exp
			ee = []
		## Knock Down expression ##
		if (grep(expr, "KnockDown")):
			KO.append([exp, condexp.get(expr)])
			continue
		## Over Expression expr. ##
		if (grep(expr, "OverExpression")):
			FE.append([exp, condexp.get(expr)])
			continue
		## Condition expression  ##
		if (grep(expr, "Conditions")):
			cond = condexp.get(expr)
		else:
			if (cond):
				ee += updateStep(condexp.get(expr), n)
				ee += updateStep(cond, n)
				if (stack):
					ee += reduce(lambda x, y : x + y, [ex for ex in stack])
					stack = []
			else:
				stack.append(updateStep(condexp.get(expr), n))
	## Empty stack          ##
	## If no condition seen ##
	if (stack):
		ee += reduce(lambda x, y : x + y, [ex for ex in stack])
		stack = []
	if (ee and e):
		E.append([e, ee])
	## Debugging            ##
	if (verbose):
		print("---- EXPERIMENTS ----")
		print("No. of experiments: " + str(len(E)))
		for k, v in {"E": E, "KO": KO, "FE": FE, "Fixpoint": Fixpoint}.items():
			print(k + " = " + str(v))
	return([E, KO, FE, Fixpoint])

#________________#
#   Re:In file   #
#________________#

#' Parse both model and experiments RE:IN files (located in the same folder)
#'
#' @param model model file name
#' @param experiments experimental file name
#' @param verbose logical for printing messages
#' @return res an instance of the GRN inference problem
def readREINfile(model="model_expanded.net", experiments="observations.spec", verbose=False):
	model, experiments = [path_to_models + e for e in [model, experiments]]
	## Extract model        ##
	with open(model, "r") as f:
		[C, CRM, length, Idef, Iopt, R, typeT, solmax, uniqueness, limreg, P] = getModel(f, verbose=verbose)
	if (verbose):
		print("_________________\n")
	## Extract experiments  ##
	with open(experiments, "r") as f:
		[E, KO, FE, Fixpoint] = getExperiments(f, C, length, verbose=verbose)
	return([C, CRM, length, Idef, Iopt, R, E, typeT, solmax, KO, FE, uniqueness, limreg, P, Fixpoint])

