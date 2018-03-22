# -*- coding: utf-8 -*-

from shortcuts import *
from utils import default, bv, testValue, equalExcept1BV, verboseIt, evalVec, strList2Str, getIt, testValueVec, ifthenelse, ifthenelseFct, extract1BV
 
######################################
## UTILS for GRN inference problem  ##
######################################

#' Build regulator/interaction dictionaries
#' keys: indices of regulators r of input gene g in the GRN
#' values: indices of associated interactions r->g if optional
#' else a default value is set
#'
#' @param Iset interaction set
#' @param gene gene which regulators shsould be considered
#' @param C set of genes/nodes
#' @param regIntActivators optionally, the partially built 
#' activator/interaction dictionary for gene
#' @param regIntRepressors optionally, the partially built 
#' repressor/interaction dictionary for gene
#' @param opt logical if set to TRUE the index of the optional
#' interaction corresponding to the regulator identifier (key) is set as value in
#' the dictionary, else a default value is set
#' @return res a list containing the activator/interaction 
#' and repressor/interaction dictionaries
def buildInteractionDict(Iset, gene, C, regIntActivators=None, regIntRepressors=None, opt=False):
    ## Initialization of dictionaries                   ##
    if (regIntActivators==None and regIntRepressors==None):
        regIntActivators = dict()
        regIntRepressors = dict()
    for i in range(len(Iset)):
        [regulator, out, sign] = Iset[i]
	## Look for interactions of type r -> gene   ##
        if (out == gene):
            idx = C.index(regulator)
	    ## value is the index of the interaction ##
	    ## if optional                           ##
	    value = i if (opt) else default
            if (sign=="+"):
                regIntActivators.setdefault(idx, value)
            else:
                regIntRepressors.setdefault(idx, value)
    return([regIntActivators, regIntRepressors])

######################################
##  SMT-CONDITIONS                  ##
######################################

## Condition: is Rs[c] == i for a given node c         ##
testRS = lambda rs, i : rs == BitVecVal(i, 5)
## Returns the string associated with a state variable ##
## for step n in experiment exp                        ##
getState = lambda n, exp : (("q_step%s_" % n) + "%s") % exp
## Returns bit-vector selecting value x[idx]           ##
## for x of size ngenes                                ##
getCState = lambda x, idx, ngenes : x & bv(idx, ngenes)

#_______________________________#
#  CONDITIONS ON uniqueness     #
#_______________________________#

#' Build condition for interaction selection uniqueness
#' for model solutions
#'
#' @param s solver
#' @param intVar bit-vector variable associated with 
#' selected optional interactions
#' @param model a solution model found be the solver
#' @return cond condition (BoolRef object) for interaction
#' uniqueness
def diff_interactions(s, intVar, model):
	m = intVar[0] if (model[intVar[0]]==None) else model[intVar[0]]
	return(intVar[0] != m)

#' Build condition of uniqueness for other cases than interaction
#' selection (i.e. difference in regulation conditions or in 
#' experiment paths)
#'
#' @param s solver
#' @param var variable list on which uniqueness should be applied 
#' (i.e. list of selected regulation conditions, list of instancied states)
#' @param model a model solution returned by the solver
#' @return cond uniqueness condition on variables in @var
def different(s, var, model):
	cond = False
	for v in var:
		cond = Or(cond, v != model[v])
	return(cond)

#_______________________________#
#  CONDITIONS ON experiments    #
#_______________________________#

#' Build the constraint on the known values of states
#' according to the experiment file
#'
#' @param s solver
#' @param q state variable list associated with one given experiment
#' @param n step in which a given state value is known
#' @param C set of genes/nodes
#' @param gene gene which state value is known
#' @param value value that appears in the experimental constraints
#' @return s updated solver with condition q[n](gene) == value
def experiment_condition(s, q, n, C, gene, value):
	s.add(testValue(q[n], C.index(gene), BitVecVal(value, 1)))
	return(s)

#____________________________#
#  CONDITIONS ON regulation  #
#____________________________#

#' Build condition on allowed regulation conditions
#'
#' @param s solver
#' @param Rs regulation condition variable list
#' @param R list of allowed regulation condition list for each gene/node
#' @param ngenes number of nodes in the abstract model
#' @return s updated solver with conditions that strictly restrict the values 
#' of a regulation condition variable to ONE of the allowed values for each gene/node
def regulation_condition(s, Rs, R, ngenes):
	for c in range(ngenes):
		cond = False
		for r in R[c]:
			cond = Or(cond, testRS(Rs[c], r))
		s.add(cond)
	return(s)

#_______________________________#
#  CONDITIONS ON perturbations  #
#_______________________________#

#' Build perturbation condition
#'
#' @param s solver
#' @param SETP list of (associated experiment, known perturbations) pairs in the experiment file
#' @param setp perturbation on which the condition should be generated
#' @param exp_names list of all experiment names in the experiment file
#' @param SETEXPR list of (gene, value) pairs in the experiment
#' @param typeP type of perturbation for the gene/node of interest
#' @return s updated solver with conditions that definitely set the known values of all perturbed gene
#' variables in a given experiment
def perturbation_condition(s, SETP, setp, exp_names, SETEXPR, typeP):
    for [e, p] in SETP:
	p_e = setp[exp_names.index(e)]
	for pp in p:
		g = SETEXPR.index(getIt(pp[1], typeP + "(", ")"))
		x = pp[2]
		s.add(testValue(p_e, g, BitVecVal(x, 1)))
    return(s)

#' Creates an equivalent Bit Vector of size #nodes
#' from a bit-vector of size #perturbed genes 
#' 
#' @param s solver
#' @param pert_e bit-vector of size #perturbed genes
#' @param chi dictionary of keys: perturbed gene indices, 
#' values: known value of gene states
#' @param name name to give to the new bit-vector perturbation variable
#' @param ngenes number of genes/nodes in the abstract model
#' @return res pair of bit-vector variable and updated solver
def pert2full(s, pert_e, chi, name, ngenes):
	pert = BitVec(name, ngenes)
	for i in range(ngenes):
		## Non-perturbed genes ##
		if (not (i in chi.keys())):
			s.add(testValue(pert, i, buildZERO(1)))
		## Perturbed genes     ##
		else:
			s.add(testValueVec(pert, pert_e, i, c1=chi.get(i)))
	return([pert, s])

#_______________________________#
#  CONDITIONS ON regulators     #
#_______________________________#

#' Build condition on the regulators in the model
#' depending on the selected optional interactions
#'
#' @param s solver
#' @param Is bit-vector of selected interactions
#' @param regulators list of present regulator variables (activators/repressors) for each node 
#' @param regInt list of regulator/interaction dictionaries for each node
#' @param default default value present in @regInt
#' @return cond condition Is[regInt[i]] == regulators[i] for all i regInt[i]!=default
def regulators_condition(s, Is, regulators, regInt, default):
	return(evalVec(s, Is, regulators, regInt, default))

#_______________________________#
#  CONDITIONS ON transitions    #
#_______________________________#

#' Build transition condition 
#' 
#' @param i regulation condition number
#' @param g gene on which the transition condition should be built
#' @param q1 output state variable of the transition
#' @param res wrapper for old state updated by input rule
#' @param ngenes number of genes
#' @param pcg pre-computation associated with gene g, respectively terms
#' allActivators, allRepressors, noRepressors, notAllActivators, notnoActivators
#' notnoRepressors, threshold regulation condition
#' @return cond transition condition q1[g] == res(rule_i)[g] \in 0, 1
def aux_transition(i, g, q1, res, ngenes, pcg):
	[aA, aR, nA, nR, naA, naR, nnA, nnR, gc] = pcg
	rulei = If(diCompute.get(i)(aA, aR, nA, nR, naA, naR, nnA, nnR, gc), bv(g, ngenes), buildZERO(ngenes))
	return(getCState(q1, g, ngenes) == getCState(res(rulei), g, ngenes))

#' Build synchronous transition condition
#'
#' @param Rs regulation condition variable list
#' @param g gene on which the transition condition should be built
#' @param i regulation condition number
#' @param auxtr function that builds the transition condition part on gene values
#' @return cond condition Rs[g] == i => gene value condition
def transition_condition_sync(Rs, g, i, auxtr):
	return(Implies(testRS(Rs[g], i), auxtr(i)))

#' Build asynchronous transition condition
#'
#' @param Rs regulation condition variable list
#' @param g gene on which the transition condition should be built
#' @param i regulation condition number
#' @param q0 input state variable for transition
#' @param q1 output state variable for transition
#' @param auxtr function that builds the transition condition part on gene values
#' @return cond condition q0[g'] == q1[g'] for g'!=g => (Rs[g] == i => gene value condition)
def transition_condition_async(Rs, g, i, q0, q1, auxtr):
	return(Implies(equalExcept1BV(q0, q1, g), transition_condition_sync(Rs, g, i, auxtr)))

#' Build the full transition condition for all genes at a given step of the experiment
#' 
#' @param s solver
#' @param n step of the experiment where q0 -> q1
#' @param length maximal length of experiment
#' @param q0 input state variable for transition
#' @param q1 output state variable for transition
#' @param preComputation list of precomputed terms for each node
#' @param ngenes number of genes/nodes
#' @param R list of allowed regulation conditions for each node
#' @param Rs list of regulation condition variables for each node
#' @param res wrapper for input state variable according to selected regulation condition
#' @param typeT type of transition, either "asynchronous" or "synchronous"
#' @param regulators list of bit-vectors variables associated with present regulators for each node
#' @return s updated solver with conditions q0 -> q1
def transition_condition(s, n, length, q0, q1, preComputation, ngenes, R, Rs, res, typeT, regulators):
	condAsync = False
	for g in range(ngenes):
		pcg = preComputation[g]
		auxtrr = lambda i : aux_transition(i, g, q1, res, ngenes, pcg)
		[a, r] = regulators[g]
		auxtr17 = lambda i : ifthenelse(i==18, rule18(q0, a, r, g), rule19(q0, a, r, g))
		for i in R[g]:
			## Different auxiliary functions are used to generate transition condition ##
			## when threshold rules are selected                                       ##
			auxtr = ifthenelse(i<18, auxtrr, auxtr17)
			if (typeT == "asynchronous"):
				cond = transition_condition_async(Rs, g, i, q0, q1, auxtr)
			elif (typeT == "synchronous"):
				cond = transition_condition_sync(Rs, g, i, auxtr)
			else:
				print("ERROR: Wrong transition type.")
				return(None)
			s.add(cond)
		if (typeT == "asynchronous"):
			condAsync = Or(condAsync, equalExcept1BV(q0, q1, g))
	if (typeT == "asynchronous"):
		s.add(condAsync)
	return(s)

#' Build the full transition condition for all genes for several steps of the experiment
#' 
#' @param s solver
#' @param expname name of the current experiment (only for printing purposes)
#' @param prepreComputation list of pre-precomputed terms for each node
#' @param startstep starting step of transition
#' @param endstep ending step of transition
#' @param q list of state variables associated with current experiment
#' @param typeT type of transition, either "asynchronous" or "synchronous"
#' @param length maximal length of experiment
#' @param regulators list of bit-vectors variables associated with present regulators for each node
#' @param ngenes number of genes/nodes
#' @param R list of allowed regulation conditions for each node
#' @param Rs list of regulation condition variables for each node
#' @param res wrapper for input state variable according to selected regulation condition
#' @param verbose logical if set to TRUE prints status messages
#' @return s updated solver with conditions q_startstep -> q_(startstep+1) -> ... -> q_endstep
def to_next_state(s, expname, prepreComputation, startstep, endstep, q, typeT, length, regulators, ngenes, R, Rs, res, verbose):
	if (typeT == "fixpoint"):
		for i in range(startstep, endstep):
			verboseIt("+++++++ FIXPOINT STEP #" + str(i) + " TO #" + str(i+1), verbose)
                	verboseIt("@ " + getState(i, expname) + " = " + getState(i+1, expname), verbose)
			s.add(q[i] == q[i+1])
	typeT = ifthenelse(typeT == "fixpoint", "synchronous", typeT)
	for i in range(startstep, endstep):
		verboseIt("******* TRAJECTORY STEP #" + str(i) + " TO #" + str(i+1), verbose)
                verboseIt("@ T(" + getState(i, expname) + ", " + getState(i+1, expname) 
			+ ") <=> " + typeT + " transition", verbose)
		preComputation = [preCompute(q[i], prepreComputation[ci]) for ci in range(ngenes)]
        	s = transition_condition(s, i, length, q[i], q[i+1], preComputation, ngenes, R, Rs, res, typeT, regulators)
	if (startstep == endstep):
		preComputation = [preCompute(q[endstep], prepreComputation[ci]) for ci in range(ngenes)]
        	s = transition_condition(s, endstep, length, q[endstep], q[endstep], 
			preComputation, ngenes, R, Rs, res, typeT, regulators)
	return(s)

#_______________________________#
#  CONDITIONS ON interactions   #
#_______________________________#

#' Build optional condition on the number of selected interactions
#'
#' @param s solver
#' @param interaction_limit integer of the maximum number of selected
#' interactions
#' @param Iopt list of optional interactions
#' @param Is bit-vector of selected interactions
#' @return s updated solver with condition on the maximum number of
#' selected interactions
def interaction_condition(s, interaction_limit, Iopt, Is):
	s.add(UGT(BitVecVal(interaction_limit+1, len(Iopt)), sub(Is)))
	return(s)

#' Implements optional condition "this gene must have at least one activator"                            
#' 
#' @param s solver
#' @param activators bit-vector associated with present activators for a given gene
#' @param ngenes number of genes
#' @return s updated solver with condition 
def mustHaveActivator_condition(s, activators, ngenes):
	s.add(UGT(activators, buildZERO(ngenes)))
	return(s)

#' Implements conditions on regulatory modules: 
#' a TF->RM interaction is selected <-> corresponding RM->gene interaction is selected
#' 
#' @param s solver
#' @param Is bit-vector variable associated with selected optional interactions
#' @param Idef set of definite interactions
#' @param Iopt set of optional interactions
#' @param CRM set of regulatory modules
#' @param C set of genes/nodes
#' @return s updated solver with conditions on regulatory modules
def crmInteractions_condition(s, Is, Idef, Iopt, CRM, C):
	for e in Idef:
		## Interaction TF->RM is definite    ##
		if (len(CRM[C.index(e[1])]) > 0):
			idx = ifthenelseFct(lambda i : Iopt[i][0]==e[1] and Iopt[i][1]==CRM[C.index(e[1])], 
				lambda i : i, Iopt)
			if (len(idx) > 0):
				idx = idx[0]
				## Selects automatically the associated RM->gene  ##
				s.add(extract1BV(Is, idx) == 1)
	for i in range(len(Iopt)):
		e = Iopt[i]
		## Interaction RM->gene is optional  ##
		if (len(CRM[C.index(e[0])]) > 0):
			## Finds all associated TF->RM interactions               ##
			idx = ifthenelseFct(lambda i : Iopt[i][1]==e[1], lambda i : i, Iopt)
			condAllNot = True
			for ii in idx:
				condAllNot = And(condAllNot, extract1BV(Is, ii) == 0)
				s.add(Implies(extract1BV(Is, ii) == 1, extract1BV(Is, i) == 1))
			s.add(Implies(condAllNot, extract1BV(Is, i) == 0))
	return(s)

#_______________________________#
#  CONDITIONS ON fix points     #
#_______________________________#

#' Build conditions on fix points
#'
#' @param s solver
#' @param expname current experiment name (only for printing purposes)
#' @param prepreComputation list of pre-precomputed terms for each node
#' @param sstep step at which starts the fix point condition
#' @param q list of state variables associated with current experiment
#' @param typeT type of transition, either "asynchronous" or "synchronous"
#' @param length maximum length of experiment
#' @param regulators list of present regulator bit-vector for each gene
#' @param ngenes number of genes/nodes
#' @param R list of allowed regulation conditions for each node
#' @param Rs list of selected regulation condition bit-vector for each node
#' @param res wrapper for input state transition
#' @param verbose logical for printing messages
#' @return s updated solver with conditions on fix points from step sstep to 
#' the end of the experiment
def fixpoint_condition(s, expname, prepreComputation, sstep, q, typeT, length, regulators, ngenes, R, Rs, res, verbose):
	if (sstep < length+1):
		s = to_next_state(s, expname, prepreComputation, sstep, length, q, 
			"fixpoint", length, regulators, ngenes, R, Rs, res, verbose)
	return(s)

