# -*- coding: utf-8 -*-

from shortcuts import *
from subprocess import call, check_output
from utils import filterNone, getIt
from global_paths import *

#' Get index of value of perturbed (of type @typeP) gene @c 
#' in result list names @resnames
#'
#' @param c gene name
#' @param resnames solver result list names
#' @param typeP type of perturbation, either KO or FE
#' @return idx index of the value of perturbed gene variable
#' in result list, None if does not exist
def aux_getPerturbationsGRFs(c, resnames, typeP):
	try:
		return(resnames.index(typeP + "(" + c + ")"))
	except:
		return(None)

#' Build the boolean functions associated with perturbated gene variables
#'
#' @param C set of genes/nodes
#' @param resnames solver result list names
#' @param typeP type of perturbation, either KO or FE
#' @return res list of boolean functions associated with all perturbed genes 
#' (of perturbation @typeP) in the GRN
def getPerturbationsGRFs(C, resnames, typeP):
	tmp = filterNone([aux_getPerturbationsGRFs(c, resnames, typeP) for c in C])
	return([typeP + "_" + getIt(resnames[t], "(", ")") for t in tmp])

## Get GRFs into readable form ##
#' Build the boolean (regulatory) functions associated with all variables
#'
#' @param C set of genes/nodes
#' @param res one solver result list of values for all logical parameter in
#' the associated model solution
#' @param regInt the list that associates to each node the dictionaries of its activators
#' and repressors, and the corresponding optional interaction indices 
#' in selected interaction vector
#' @return grfs list of GRFs for each node of the GRN
def get_grfs(C, res, regInt):
	grfs = []
	sep = "\' = "
	ko_g = lambda g : "KO_" + C[g]
	fe_g = lambda g : "FE_" + C[g]
	## Get regulation template fct no.   ##
	resnames = [str(x[0]) for x in res]
	ko = getPerturbationsGRFs(C, resnames, "KO")
	fe = getPerturbationsGRFs(C, resnames, "FE")
	## Write perturbation functions      ##
	for kog in ko:
		grfs.append(kog + sep + writePerturbedGene(kog))
	for feg in fe:
		grfs.append(feg + sep + writePerturbedGene(feg))
	for idx in range(len(C)):
		## Get activators and repressors for gene C[idx]        ##
		## regInt[idx] = list of two dictionaries               ##
		## keys: index of regulator in C                        ##
		## values: index of corresponding interaction in        ##
		## optional interaction vector                          ##
		[regIntActivators, regIntRepressors] = regInt[idx]
		activators = res[resnames.index('activators_' + C[idx])][1]
		repressors = res[resnames.index('repressors_' + C[idx])][1]
		## Get no. for current gene                             ## 
		r = int(res[resnames.index('grf_' + C[idx])][1])
		## Write corresponding GRF                              ##
		grf = C[idx] + sep + writePerturbation((ko_g(idx) in ko), 
				(fe_g(idx) in fe), 
				diGRF(0).get(r)(activators, repressors, C[idx]),
				ko_g(idx),
				fe_g(idx))
		grfs.append(grf)
	return(grfs)

#' Write GRFs in a text file named title"_model"#.txt in 
#' folder path_to_results
#'
#' @param grfsList list of GRF list (for each model solution)
#' @param title result file name
#' @return None 
def write_grfs(grfsList, title):
	nmodels = len(grfsList)
	for m in range(1, nmodels+1):
		with open(path_to_results + title + "_model" + str(m) + ".txt", "w") as f:
			for grf in grfsList[m-1]:
				f.write(grf + "\n")
	return(None)

#' Simplify GRFs in text files named title located in 
#' folder path_to_results using reducer in folder path_to_reducer
#' Resulting file is named "simplified--"title
#'
#' @param title result file name in which GRFs should be simplified
#' @return None 
def simplify_grfs(title):
	print("MSG: Simplify file " + title)
	check_output("R -e \'suppressMessages(source(\"" + path_to_reducer + "boolean_reducer.R\"));"
		+ " res = boolean_reducer(\""
		+ path_to_results + title 
		+ "\", title=\"" + path_to_results + "simplified--" + title + "\")\'", shell=True)
	return(None)
