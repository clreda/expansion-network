# -*- coding: utf-8 -*-
 
from grn_inference import *
from utils import *
from time import time 

#######################
## MODEL PROCESSING  ##
#######################
  
#' Search the value of present regulator bit-vector variable
#'
#' @param C set of genes/nodes
#' @param regulatorsVar list of present regulator bit-vector variables
#' @param M model returned by the solver
#' @return res list of pairs (bit-vector variable, value in model @M)
def getPresentRegulators(C, regulatorsVar, M):
    return([[v, filterNone([ifthenelse(simplify(testValue(M[v], i, true)), C[i]) for i in range(len(C))])] for v in regulatorsVar])
 
#' Give the list of perturbed genes in the experiments
#'
#' @param C set of genes/nodes
#' @param chi dictionary of perturbed genes
#' @param typeP type of considered perturbation (KO or FE)
#' @return res list of known perturbed genes in the experiments
def getPerturbedGenes(C, chi, typeP):
    if (chi):
        return([[typeP + "(" + C[i] + ")", typeP + "(" + C[i] + ")"] for i in chi.keys()])
    return([])
 
#' Generate a model from a solution returned by the solver
#' 
#' @param C set of genes/nodes
#' @param length maximum length of experiment
#' @param Iopt set of optional interactions
#' @param ngenes number of genes/nodes
#' @param t timestamp of computation start
#' @param s solver
#' @param intVar list containing Is, Is is the bit-vector variable
#' associated with the selected optional interactions
#' @param regVar list of #nodes Rs_v, Rs_v is the bit-vector variable
#' associated with the selected regulation condition for node v
#' @param stateVar list of list of state variables at each step for each experiment 
#' @param regulators list of present regulator bit-vector for each gene
#' @param chiDOWN dictionary associated with KO perturbations
#' @param chiUP dictionary associated with FE perturbations
#' @param R list of allowed regulation conditions for each node
#' @param uniqueness string for uniqueness condition to find the next solution:
#' "interaction", "full", or "paths"
#' @param verbose logical for printing messages
#' @param resList current list of model solutions
#' @param solNB number of the currently processed solution
#' @param stopSol logical for not applying the uniqueness condition at the end of the processing
#' @param printSolutions logical for printing message about the processed solution
#' @return res list of updated solution list, new timestamp, updated solver with uniqueness
#' condition, logical indicating if the solver should stop looking for conditions
def getSolution(C, length, Iopt, ngenes, t, s, intVar, regVar, stateVar, regulators, chiDOWN, chiUP, R, uniqueness, verbose, resList=[], solNB=1, stopSol=False, printSolutions=True):
    lenIopt = len(Iopt)
    verboseIt("\nTIME = %.2f" % (time()-t) + " sec\n", printSolutions)
    t = time()
    s.check()
    verboseIt("CHECKING TIME = %.2f" % (time()-t) + " sec\n", printSolutions)
    noMoreModel = False
    M = None
    try:
        M = s.model()
    except:
        if (resList):
            verboseIt("No other model found.\n", printSolutions)
        else:
            verboseIt("No model found.\n", True)
        noMoreModel = True
	return([resList, t, s, noMoreModel])
    if (M):
        dec = lambda lsVar, size : [[v, getBinaryDec(M[v], size)] for v in lsVar]
        intSol = dec(intVar, lenIopt) if (lenIopt) else []
        stateSol = dec(stateVar, ngenes)
        regSol = dec(regVar, None)
        actVar, repVar = [x[0] for x in regulators], [x[1] for x in regulators]
        actSol = getPresentRegulators(C, actVar, M)
        repSol = getPresentRegulators(C, repVar, M)
        koSol = getPerturbedGenes(C, chiDOWN, "KO")
        feSol = getPerturbedGenes(C, chiUP, "FE")
        solution = intSol + regSol + stateSol + actSol + repSol + koSol + feSol
        resList.append(solution)
        verboseIt("Model no. " + str(solNB) + " found:", printSolutions)
        if (lenIopt and printSolutions):
            verboseIt("> Interaction vector: ", True)
            verboseIt(reduce(lambda x,y: x+y, [str(x) for x in rev(Iopt)]), True)
            printPretty(intSol)
	if (printSolutions):
        	verboseIt("> GRFs: ", True)
        	printPretty(regSol)
        if (verbose and printSolutions):
            verboseIt("> States: ", True)
            for n in range(len(stateSol)/(length+1)):
                verboseIt(">> Experiment: ", True)
                verboseIt(reduce(lambda x,y: x+y, rev(C)), True)
                printPretty(stateSol[n*(length+1):(n+1)*(length+1)])
        verboseIt("____________________________________\n", printSolutions)
	if (stopSol):
		return([resList, t, s, True])
        ## Uniqueness of models: interactions                ##
        ## Add condition newIs != Is                         ##
        ## (Is != value of Is in model M)                    ##
        ## "interactions": unique interactions               ##
        ## "full": unique interactions OR unique set of      ##
	## regulation conditions                             ##
        ## "paths": unique interactions OR unique set of     ##
	## regulation conditions OR set of trajectories      ##
	cond1 = diff_interactions(s, intVar, M)
	if (uniqueness == "interactions"):
		s.add(cond1)
		return([resList, t, s, noMoreModel])
	cond2 = different(s, regVar, M)
	if (uniqueness == "full"):
		s.add(Or(cond2, cond1))
		return([resList, t, s, noMoreModel])
	if (uniqueness == "paths"):
		cond3 = different(s, stateVar, M)
		s.add(Or(cond1, Or(cond2, cond3)))
		return([resList, t, s, noMoreModel])
	verboseIt("MSG: Warning! No correct uniqueness condition detected.", True)
	verboseIt("MSG: Current uniqueness condition: \'" + uniqueness + "\'.", True)
	verboseIt("MSG: Uniqueness conditions can only be in [\'interactions\', \'full\', \'paths\'].\n", True)
	noMoreModel = True
	return([resList, t, s, noMoreModel])
 
#######################
## NETWORK INFERENCE ##
#######################

#' Solve an instance of the GRN inference problem
#' 
#' Given the network topology and the regulation
#' function types for each node of the network, 
#' find a subset of optional interactions and 
#' of regulation function types that comply with
#' the provided experiments
#' 
#' @param C node names in the network 
#'          (character string list)
#' @param CRM for each node (eventually) the gene 
#' 	      it directly regulates (if the node is a CRM)
#' @param length maximum experiment length 
#'          (integer)
#' @param Idef set of definite interactions 
#'             (list of lists of character strings)
#' @param Iopt set of optional interactions
#'             (list of lists of character strings)
#' @param R set of set of possible regulation 
#'          conditions for each node
#'          (list of lists of integers)
#' @param E set of experiments: 
#'          list of [experiment name, list of (step x gene x GE value)]
#'          (list of string x (list of (integer x character string x binary integer)))
#' @param typeT either sync (synchronous) or async (asynchronous)
#'              transition
#'              (character string)
#' @param solmax maximal number of solutions to return
#'               (integer)
#' @param KO knock-down perturbations
#'           (list of [experiment name, list of (step x KO(gene) x binary)])
#' @param FE forced expression perturbations
#'           (list of [experiment name, list of (step x FE(gene) x binary)])
#' @param uniqueness character string (see getSolution function)
#' @param limreg character string (to limit regulation functions)
#' @param P perturbation list
#'          (list of size |C| containing lists with 0 element, or "-" or "+" or both)
#' @param Fixpoint fixpoint constraints
#'                 (list of [start step of fix point state, experiment name])
#' @param verbose boolean: prints status messages if set to True
#' @return resList list of models where Is and Rs are 
#'                 the instanciated constrained ABN
#'                 that agree with all the experiments (+ solver)
def grn_solver(C, CRM, length, Idef, Iopt, R, E, typeT, solmax, KO, FE, uniqueness, limreg, P, Fixpoint, verbose=False, printSolutions=True, printmain=True, maximize=False):
    ## Selected interaction number limit                ##
    interaction_limit = 0
    if (not interaction_limit and Iopt):
        interaction_limit = len(Iopt)
    #____________________________________________________#
    #  Initialization of constants and variables         #
    #____________________________________________________#
    if (not solmax):
        solmax = 10
    if (maximize):
	s = Optimize()
    else:
        s = Solver()
    ngenes = len(C)
    zero, one = buildZERO(ngenes), buildONE(ngenes)
    mustHaveActivator = []
    UP, DOWN = [], []
    chiUP, chiDOWN = dict(), dict()
    for i in range(len(P)):
        if ("-" in P[i]):
            DOWN.append(C[i])
            chiDOWN.setdefault(i, len(DOWN)-1)
        if ("+" in P[i]):
            UP.append(C[i])
            chiUP.setdefault(i, len(UP)-1)
        if ("!" in P[i]):
            mustHaveActivator.append(i)
    ## Variable for the subset of optional interactions ##
    Is = BitVec("selected_interactions", len(Iopt)) if (Iopt) else false
    if (maximize):
	from utils import sub
	lIs = BV2Int(sub(Is))
	s.add(lIs > 0)
	s.add(lIs <= len(Iopt))
	objective = s.maximize(lIs)
    intVar = [Is]
    ## Variables for the regulation functions           ##
    Rs = [BitVec("grf_%s" % node, 5) for node in C] 
    regVar = Rs
    ## Variables for regulators                         ##
    regulators = [[BitVec("activators_%s" % gene, ngenes), 
        BitVec("repressors_%s" % gene, ngenes)] for gene in C]
    ## Variables for perturbations                      ##
    ko = [BitVec("ko_%s" % e[0], len(DOWN)) for e in E] if (DOWN) else []
    fe = [BitVec("fe_%s" % e[0], len(UP)) for e in E] if (UP) else []
    regInt = []
    stateVar = []
    exp_names = [e[0] for e in E]
    t = time()
    #____________________________________________________#
    #  Conditions on regulation functions                #
    #____________________________________________________#
    verboseIt("Conditions on regulation functions", verbose)
    s = regulation_condition(s, Rs, R, ngenes)
    #____________________________________________________#
    #  Conditions on perturbations                       #
    #____________________________________________________#
    s = perturbation_condition(s, KO, ko, exp_names, DOWN, "KO") if (ko and KO) else s
    s = perturbation_condition(s, FE, fe, exp_names, UP, "FE") if (fe and FE) else s
    #____________________________________________________#
    #  Conditions on regulators                          #
    #____________________________________________________#
    verboseIt("Computation of interactions", verbose)
    if (any([len(c) > 0 for c in CRM])):
	s = crmInteractions_condition(s, Is, Idef, Iopt, CRM, C)
    if (Iopt):
        s = interaction_condition(s, interaction_limit, Iopt, Is)
    for gene in C:
        idx = C.index(gene)
        ## Build lists of (indices of)    ##
        ## activators and                 ##
        ## repressors for input gene      ##
        ## For interaction edges in Iopt  ##
        ## Keep track of the regulatory   ##
        ## interactions for the gene      ##
        regIntActivators, regIntRepressors = None, None
        if (Iopt):
            [regIntActivators, regIntRepressors] = buildInteractionDict(Iopt, gene, C, opt=True)
        ## For interaction edges in Idef  ##
        if (Idef):
            [regIntActivators, regIntRepressors] = buildInteractionDict(Idef, gene, C, 
            regIntActivators=regIntActivators, regIntRepressors=regIntRepressors)
        regInt.append([regIntActivators, regIntRepressors])
        ## Potential regulator is a       ##
        ## regulator iff. the interaction ##
        ## where it is involved is        ##
        ## selected or definite           ##
        activators, repressors = regulators[idx]
        s = regulators_condition(s, Is, activators, regIntActivators, default)
        s = regulators_condition(s, Is, repressors, regIntRepressors, default)
        if (idx in mustHaveActivator):
            s = mustHaveActivator_condition(s, activators, ngenes)
        regulatorList = lambda regulators, regInt : [regulators+
            " of gene "+gene+": "]+["gene "+C[i]+", interaction "
            + str(regInt[i]) + "; " for i in regInt.keys()]
        verboseIt(strList2Str(regulatorList("ACTIVATORS", regIntActivators)), verbose)
        verboseIt(strList2Str(regulatorList("REPRESSORS", regIntRepressors)), verbose)
    verboseIt("Computation of GRFs", verbose)
    prepreComputation = [prepreCompute(regulators[ci][0], regulators[ci][1]) for ci in range(len(C))]
    #____________________________________________________#
    #  Conditions on experiments                         #
    #____________________________________________________#
    verboseIt("Conditions on experiments", verbose)
    for exp in E:
        verboseIt("--------- EXPERIMENT \'" + exp[0] + "\'", verbose=printmain)
        ## State variables                  ##
        q = [BitVec(getState(n, exp[0]), ngenes) for n in range(length+1)]
        stateVar += q
        ## Adding KO and FE constraints     ## 
        if (KO and FE and ko and fe):
            [ko_e, s] = pert2full(s, ko[exp_names.index(exp[0])], chiDOWN, "ko_" + exp[0] + "_f", ngenes)
            [fe_e, s] = pert2full(s, fe[exp_names.index(exp[0])], chiUP, "fe_" + exp[0] + "_f", ngenes)
            res = lambda x : (x & ~ko_e) | fe_e
        elif (KO and ko):
            [ko_e, s] = pert2full(s, ko[exp_names.index(exp[0])], chiDOWN, "ko_" + exp[0] + "_f", ngenes)
            res = lambda x : x & ~ko_e
        elif (FE and fe):
            [fe_e, s] = pert2full(s, fe[exp_names.index(exp[0])], chiUP, "fe_" + exp[0] + "_f", ngenes)
            res = lambda x : x | fe_e
        else:
            res = lambda x : x
        #____________________________________#
        ## States must define a trajectory  ##
        ## in the search space              ##
        #____________________________________#
	existsFixpoint = filter(lambda x : x[1] == exp[0], Fixpoint)
	## Finds the starting step point    ##
	## for fix point                    ##
	if (existsFixpoint):
		sstep = existsFixpoint[0][0]
	else:
		sstep = None
	## Enforces constraint for all i,   ##
	## 0 <= i < sstep, T(q[i],q[i+1])   ##
	s = to_next_state(s, exp[0], prepreComputation, 0, ifthenelse(sstep == None, length, sstep), q, typeT, length, regulators, ngenes, R, Rs, res, verbose)
        ## Fixpoint constraint              ##
	## for all i, sstep <= i <= length  ##
	## T(q[i], q[i+1]) (synchronous) and##
	## q[i] = q[i+1]                    ##
        s = fixpoint_condition(s, exp[0], prepreComputation, ifthenelse(sstep==None, length+1, sstep), q, typeT, length, regulators, ngenes, R, Rs, res, verbose)
        #____________________________________#
        ## Experiment values should be      ##
        ## satisfied                        ##
        #____________________________________#
        ## For each observation in e        ##
        ## ee = { n, gene, value }          ##
        for [n, gene, value] in exp[1]:
            verboseIt("Experiment=\'" + exp[0] + "\', Step="
                + str(n) + ": grf(" + gene + ")=" + str(value), verbose)
            s = experiment_condition(s, q, n, C, gene, value)
    #____________________________________#
    ## Solution processing              ##
    #____________________________________#
    [resList, t, s, res] = getSolution(C, length, Iopt, ngenes, t, s, intVar, 
        regVar, stateVar, regulators, chiDOWN, chiUP, R, uniqueness, verbose, 
		resList=[], stopSol=(solmax==1), printSolutions=printSolutions)
    if (not len(resList)):
        return([resList, s, regInt])
    sol = 2
    while (sol <= solmax and not res):
        [resList, t, s, res] = getSolution(C, length, Iopt, ngenes, t, s, intVar, 
        	regVar, stateVar, regulators, chiDOWN, chiUP, R, uniqueness, 
			verbose, resList=resList, solNB=sol, stopSol=(solmax == sol), printSolutions=printSolutions)
        if (res):
            return([resList, s, regInt])
        sol += 1
    if (sol == solmax+1):
	verboseIt("Maximum number of solutions reached.\n", printSolutions)
    if (sol > 0):
	verboseIt("There are solutions.\n", printSolutions)
    return([resList, s, regInt])
