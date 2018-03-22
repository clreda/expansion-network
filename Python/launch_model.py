# -*- coding: utf-8 -*-

from grn_solver import grn_solver
from grn_inference import getState
from utils import rev, verboseIt, ifthenelse
from igraph import *
from copy import deepcopy
 
##############################################
## Display the graph visualization of a     ##
## model found by the solver                ##
##############################################
 
#' Display the Boolean Network representation 
#' of the model: nodes are annotated by gene name and
#' selected regulation function, only selected optional
#' and definite interactions are drawn
#' - nodes are coloured in green if they can be FE
#' - nodes are coloured in red if they can be KO
#' - nodes are coloured in yellow if they can be FE or KO
#'
#' @param modelID integer identifier of model in resList
#' @param resList result from the solver containing
#'                the values for each variable
#' @param C the list of nodes
#' @param Idef1 the list of definite interactions: 
#'             regulator x output gene x sign
#' @param Iopt1 the list of optional interactions:
#'             regulator x output gene x sign
#' @param P the list of possible perturbations for every gene
#' @param plotIt boolean for plotting the igraph
#' @param verbose boolean for having comments
#' @return res igraph associated to the model
def model2igraph(modelID, resList, C, Idef1, Iopt1, P, plotIt=False, verbose=False):
    Idef = deepcopy(Idef1)
    Iopt = deepcopy(Iopt1)
    ngenes = len(C)
    resList = resList[modelID]
    Is = resList[0][1]
    I = Idef
    nIdef = len(Idef)
    Iopt = rev(Iopt)
    for i in range(len(Iopt)):
        if (int(Is[i])==1):
            I.append(Iopt[i])
    edges = [[C.index(y) for y in x[:2]] for x in I]
    verboseIt("Edges: " + str([[C[i] for i in e] for e in edges]), verbose)
    elabels = [x[2] for x in I]
    verboseIt("Edge signs: " + str(elabels), verbose)
    color = ["white"]*ngenes
    for i in range(len(P)):
        condKO = "-" in P[i]
        condFE = "+" in P[i]
        if (condKO and condFE):
            color[i] = "yellow"
        elif (condKO):
            color[i] = "red"
        elif (condFE):
            color[i] = "green"
    verboseIt("Node colors: " + str([C[i] + " is " + color[i] for i in range(ngenes)]), verbose)
    grf = [str(x[1]) for x in resList[1:(ngenes+1)]]
    verboseIt("Node GRFs: " + str([C[i] + " uses rule #" + grf[i] for i in range(ngenes)]), verbose)
    vertex_attrs = {"label": [C[i] + " (" + grf[i] + ")" for i in range(ngenes)], 
            "size" : [60+20*(len(c)/2-1) for c in C]*ngenes,
            "color": color}
    edge_attrs = {"label": elabels,
            "color": [ifthenelse(el == "-", "red", "green") for el in elabels],
	    "lty": [0]*nIdef + [0]*(len(I)-nIdef)}
    g = Graph(n = ngenes, 
        edges = edges, 
        directed = True, 
        vertex_attrs = vertex_attrs, 
        edge_attrs = edge_attrs)
    ## Delete 0-degree vertices       ##
    g.vs.select(_degree = 0).delete()
    if (plotIt):
	plot(g, layout = g.layout("kk"), bbox = ((g.vcount()+1)*100, (g.vcount()+1)*100), margin = 100)
    else:
	out = plot(g, layout = g.layout("kk"), bbox = ((g.vcount()+1)*100, (g.vcount()+1)*100), margin = 100)
	out.save("grn_model" + str(modelID) + ".png")
    return(g)
 
##############################################
## From a given initial state and a model,  ##
## perform one or several steps of          ##
## transition and return final (steady)     ##
## state(s)                                 ##
##############################################
 
#' Perform ("unfold") one or several transition
#' steps from a given model and an initial state
#'
#' @param modelID integer identifier of model in resList
#' @param C the list of nodes in the model
#' @param resList result from the solver containing
#'                the values for each variable
#' @param Idef the list of definite interactions: 
#'             regulator x output gene x sign
#' @param Iopt the list of optional interactions:
#'             regulator x output gene x sign
#' @param R list of allowed regulatory functions for each node
#' @param q0 initial state (list of integers 0,1 or booleans, or string character
#'           of 0's and 1's, in the order of nodes in C)
#' @param nstep number of steps to perform
#' @param typeT type of transition: either "synchronous" or "asynchronous"
#' @param KO list of KO experiments
#' @param FE list of FE experiments
#' @param P list of possible perturbations for every gene
#' @param solmax maximum number of trajectories to return
#' @param steadyStates boolean for obtaining the steady states
#' @param verbose boolean for having comments
#' @return trajectories a list of lists of pairs (state x string character) 
#'                     (the string character elements are ordered according to the order in C)
def launch_model(modelID, C, CRM, resList, Idef, Iopt, R, q0, nstep, typeT, KO, FE, P, 
		solmax=10, steadyStates=False, verbose=False):
	resList = resList[modelID]
	## Computing the whole set of interactions I for the given model ##
	Is = resList[0][1]
	I = Idef
	Iopt = rev(Iopt)
	for i in range(len(Iopt)):
		if (int(Is[i])):
			I.append(Iopt[i])
	## Building the constraints for the initial state                ##
	exp_name = "Experiment"
	E = [[exp_name, [[0, C[i], int(q0[i])] for i in range(len(C))]]]
	Fixpoint = ifthenelse(steadyStates, [[nstep, exp_name]], [])
	grfsres = [str(x[0]) for x in resList]
	R = [[resList[grfsres.index("grf_" + c)][1]] for c in C]
	## Get solmax trajectories from the given initial state          ##
	res = grn_solver(C, CRM, nstep, I, [], R, E, typeT, solmax, KO, FE, "paths", "", P, Fixpoint, 
			verbose=verbose, printSolutions=verbose, printmain=verbose)
	try:
		[resList, _, _] = res
	except:
		verboseIt("No model has been returned", verbose=True)
		return(res)
	trajectories = []
	for k in range(len(resList)):
		res = resList[k]
		statesnames = [getState(i, "Experiment") for i in range(nstep+1)]
		resnames = [str(x[0]) for x in res]
		statesIDX = [resnames.index(sn) for sn in statesnames]
		trajectory = [[resnames[i], rev(res[i][1])] for i in statesIDX]
		trajectories.append(["Trajectory #" + str(k+1), trajectory])
	return(trajectories)
