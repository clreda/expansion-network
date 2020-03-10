# -*- coding: utf-8 -*-

import sys
from z3 import *
from launch_model import *
from utils import *
from global_paths import *
from grn_inference import getCState, regulation_condition, testRS, aux_transition, transition_condition_sync, transition_condition_async

##############################################
## SOLVE MODEL call                         ##
##############################################
 
def call_run(model=full_toy_model, experiments=full_toy_experiments, simplify=False, visualize=False):
    from models import readREINfile
    from grn_solver import grn_solver
    from get_grfs import get_grfs, write_grfs, simplify_grfs
    if (visualize):
    	from launch_model import model2igraph
    [C, CRM, length, Idef, Iopt, R, E, typeT, solmax, KO, FE, uniqueness, 
            limreg, P, Fixpoint] = readREINfile(model, experiments)
    res = grn_solver(C, CRM, length, Idef, Iopt, R, E, typeT, solmax, 
            KO, FE, uniqueness, limreg, P, Fixpoint, verbose=False)
    try:
        [resList, s, regInt] = res
    except:
        return(res)
    title = "result_" + model.split("/")[0] + "_" + model.split("/")[1]
    grfsList = [get_grfs(C, res, regInt) for res in resList]
    write_grfs(grfsList, title=title)
    print("MSG: Wrote GRFs in " + title + ".txt")
    if (simplify):
        for m in range(1, len(grfsList)+1):
            simplify_grfs(title=title + "_model" + str(m) + ".txt")
	print("MSG: Simplified functions.")
    if (visualize):
        for modelID in range(len(resList)):
            print("MSG: Close plot to resume.")
            model2igraph(modelID, resList, C, Idef, Iopt, P, plotIt=True)
    return(s)
 
def printRunSyntaxError(c):
    if (not c):
        print("MSG: If you wanted to run the solver, then you did not use the correct syntax.")
        print("MSG: Correct syntax is \'run [--simplify] [--visualize] model experiments\'.")
        return(True)
    return(False)
 
if (len(sys.argv) > 1 and sys.argv[1] != "launch"):
	if (len(sys.argv) == 2):
	    cond1 = sys.argv[1] == "run"
	    if (not printRunSyntaxError(cond1)):
		call_run()
	if (len(sys.argv) > 2):
	    cond1 = sys.argv[1] == "run"
	    if (not printRunSyntaxError(cond1)):
		## Tests toy model                   ##
		if (len(sys.argv) == 2):
		    call_run()
		## Tests toy model w/ 1 option       ##
		if (len(sys.argv) == 3):
		    cond2 = sys.argv[2] == "--simplify"
		    cond3 = sys.argv[2] == "--visualize"
		    if (not printRunSyntaxError(cond2 or cond3)):
		        call_run(simplify=cond2, visualize=cond3)
		## Tests other models w/o options    ##
		if (len(sys.argv) == 4):
		    cond4 = sys.argv[2] != "--simplify"
		    cond5 = sys.argv[2] != "--visualize"
		    cond6 = sys.argv[3] != "--simplify"
		    cond7 = sys.argv[3] != "--visualize"
		    if (not printRunSyntaxError((cond4 and cond5 and cond6 and cond7) or ((not cond4 and not cond7) or (not cond6 and not cond5)))):
		        if (cond4 and cond5 and cond6 and cond7):
		            call_run(model=sys.argv[2], experiments=sys.argv[3])
		        ## Tests toy model w/ 2 options      ##
		        else:
		            call_run(simplify=True, visualize=True)
		## Tests other models w/ 1 option    ##
		if (len(sys.argv) == 5):
		    cond12 = sys.argv[2] == "--simplify"
		    cond13 = sys.argv[2] == "--visualize"
		    cond14 = sys.argv[3] != "--simplify"
		    cond15 = sys.argv[3] != "--visualize"
		    cond16 = sys.argv[4] != "--simplify"
		    cond17 = sys.argv[4] != "--visualize"
		    if (not printRunSyntaxError((cond12 or cond13) and cond14 and cond15 and cond16 and cond17)):
		        if (cond12):
		            call_run(model=sys.argv[3], experiments=sys.argv[4], simplify=True)
		        else:
		            call_run(model=sys.argv[3], experiments=sys.argv[4], visualize=True)
		## Tests other models w/ 2 options   ##
		if (len(sys.argv) == 6):
		    cond18 = sys.argv[2] == "--simplify"
		    cond19 = sys.argv[3] == "--visualize"
		    cond20 = sys.argv[3] == "--simplify"
		    cond21 = sys.argv[2] == "--visualize"
		    cond22 = sys.argv[4] != "--simplify"
		    cond23 = sys.argv[4] != "--visualize"
		    cond24 = sys.argv[5] != "--simplify"
		    cond25 = sys.argv[5] != "--visualize"
		    if (not printRunSyntaxError(((cond18 and cond19) or (cond20 and cond21)) and cond22 and cond23 and cond24 and cond25)):
		        call_run(model=sys.argv[4], experiments=sys.argv[5], simplify=True, visualize=True)
		if (len(sys.argv) > 6):
		    printRunSyntaxError(False)

##############################################
## LAUNCH MODEL call                        ##
##############################################

def getArgument(x, argv, default):
	xx = "--" + x
	if (any([arg == xx for arg in argv])):
		i = argv.index(xx)
		if (len(argv) < i+2):
			return(default)
		return(argv[i+1])
	return(default)

def getTrailingArguments(x, argv, default):
	xx = "--" + x
	if (any([arg == xx for arg in argv])):
		i = argv.index(xx)
		if (len(argv) < i+2):
			return(default)
		return(argv[i+1:])
	return(default)

from models import *

def getConditionsExp(experiments):
	with open(path_to_models + experiments, "r") as f:
		[_, summary, condexp] = readConditions(splitNclean(f.read(), ";"))
	return(condexp)

if (len(sys.argv) > 1 and sys.argv[1] == "launch"):
	if (len(sys.argv) == 2):
		print("MSG: If you wanted to launch the solver, then you probably forgot the model name.")
		print("MSG: Correct syntax is \'launch model_name [options]\'.")
	else:
		emodel = getArgument("model", sys.argv, "model_expanded")
		model = sys.argv[2] + "/" + emodel + ".net"
		experiments = sys.argv[2] + "/" + getArgument("experiments", sys.argv, "observations") + ".spec"
		## To modify for plotting graph resulting from one of the models in resList ##
		## For examples, see file resList_test.py in "examples/"                    ##
		sys.path.insert(0, '../examples')
		from resList_test import *
		if (sys.argv[2]=="pluripotency"):
			resList = resList_Dunn_regular if (emodel=="model_expanded") else resList_Dunn_expanded
		elif (sys.argv[2]=="collombet"):
			resList = resList_Collombet_regular if (emodel=="model_expanded") else resList_Collombet_expanded
		else:
			resList = []
		[C, CRM, length, Idef, Iopt, R, E, typeT, solmax, KO, FE, uniqueness, 
			limreg, P, Fixpoint] = readREINfile(model, experiments)
		if (len(sys.argv) > 3 and "igraph" in sys.argv):
			## Avoids adding colours to the nodes according to their perturbations      ##
			P = [""]*len(C)
			if (not resList):
				addGRF = False
				R = [["?"]]*len(C)
				resList = [[['Is', [1]*len(Iopt)]] + [[C[i], R[i][0]] for i in range(len(C))]]
			else:
				addGRF = True
			## Boolean addGRF is set to True iff. indices corresponding to GRF  ##
			## templates should be printed                                      ## 
			model2igraph(0, resList, C, Idef, Iopt, P, model=sys.argv[2] + "_" + emodel, plotIt=True, addGRF=addGRF)
		else:
			print("-- START")
			if (len(Iopt) > 0):
				print("Solving abstract model...")
				[resList, _, _] = grn_solver(C, CRM, length, Idef, Iopt, R, E, typeT, 
					solmax, KO, FE, uniqueness, limreg, P, Fixpoint, printSolutions=False)
				print("... done!")
			else:
				resList = [[['Is', []]] + [["grf_"+C[i], R[i][0]] for i in range(len(C))]]
			condexp = getConditionsExp(experiments)
			idm = int(getArgument("modelID", sys.argv, 0))
			modelID = idm if (not idm) else idm-1
			q0 = getArgument("q0", sys.argv, "1"*len(C))
			nstep = getArgument("nstep", sys.argv, length)
			solmax = getArgument("solmax", sys.argv, 10)
			steadyStates = getArgument("steadyStates", sys.argv, 0)
			expnames = getTrailingArguments("expnames", sys.argv, [])
			if (q0 in condexp.keys()):
				q0 = condexp.get(q0)
			## In order to add long-lasting perturbations
			ko = getArgument("KO", sys.argv, "")
			fe = getArgument("FE", sys.argv, "")
			if (ko in condexp.keys()):
				KO = [[ko, condexp.get(ko)]]
			else:
				KO = []
			if (fe in condexp.keys()):
				FE = [[fe, condexp.get(fe)]]
			else:
				FE = []
			trajectories = launch_model(modelID, C, CRM, resList, Idef, Iopt, R, 
					q0, int(nstep), typeT, KO, FE, P, int(solmax), 
					steadyStates=bool(steadyStates))
			print("----------------------------------------------------------------")
			print("modelID = " + str(modelID) + "; q0 = {" + (reduce(lambda x,y : x+", "+y, list(map(lambda x: x[1]+"="+str(x[2]), q0))) if (str(type(q0)) == "<type 'list'>") else q0) + "} ; nstep = " + str(nstep))
			if (KO):
				print("KO perturbations = { "+reduce(lambda x,y: x+", "+y, list(filter(lambda x :x, [x[1][3:-1] if (x[2] > 0) else None for x in KO[0][1]])))+" }")
			if (FE):
				print("FE perturbations = { "+reduce(lambda x,y: x+", "+y, list(filter(lambda x :x, [x[1][3:-1] if (x[2] > 0) else None for x in FE[0][1]])))+" }")
			print("")
			npaths = len(trajectories)
			print("#trajectories = " + str(npaths))
			print("----------------------------------------------------------------")
			chunksC = [C[i:i + 10] for i in xrange(0, len(C), 10)]
			for i in range(npaths):
				printStates(trajectories, i, C)
				path = trajectories[i][1]
				print("\n")
				for exp in expnames:
					if (exp in condexp.keys()):
						summary = condexp.get(exp)
						nodesIDX = [C.index(x[1]) for x in summary]
						values = concat([str(x[2]) for x in summary])
						appear = False
						for j in range(len(path)):
							valuesPath = concat([path[j][1][idx] for idx in nodesIDX])
							if (values == valuesPath):
								print("Condition \'" + exp 
									+ "\' appears at step " + str(j) + ".")
								appear = True
						if (not appear):
							print("Condition \'" + exp + "\' does not appear in trajectory.")
					else:
						print("The condition \'" + exp + "\' does not exist.")
				print("\n_______________________________________________\n\n")
			print("--END")
