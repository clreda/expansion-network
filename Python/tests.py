# -*- coding: utf-8 -*-
 
import sys
from z3 import *
from random import randint
from utils import *
from grn_inference import getCState, regulation_condition, testRS, aux_transition, transition_condition_sync, transition_condition_async
 
###################
## TOOLS         ##
###################
 
is_test = lambda fct : len(sys.argv) == 1 and sys.argv[1] == (fct)
 
generate_bv = lambda size : filter(lambda x : x >= 0, [i if (randint(0, 1)==1) else -1 for i in range(size)])
 
def generate_state(size, prelist=[]):
    if (prelist == []):
        prelist = generate_bv(size)
    if (not len(prelist)):
        return(buildZERO(size))
    return(intList2BV(prelist, size))
 
def generate_regulators(idx, size, noRegulator=False, selfRegulator=False, prelist=[], prohibitedIdx=[]):
    if (prelist == []):
        prelist = generate_bv(size)
    if (not len(prelist)):
        return(buildZERO(size))
    regulators = buildZERO(size)
    if (noRegulator):
        return(regulators)
    regulators = intList2BV(filter(lambda i : i!=idx and not (i in prohibitedIdx), prelist), size)
    if (selfRegulator):
        regulators = regulators | bv(idx, size)
    return(regulators)
 
def custom_int2regulators(lsi, size):
    lsi = [[ifthenelse(int(rev(i)[k]), k, None) for k in range(len(i))] for i in lsi]
    lsi = [filterNone(ls) for ls in lsi]
    return([intList2BV(ls, size) for ls in lsi])
 
## These operations are NOT bit-wise ##
printRule = {
         0: "allActivators & noRepressors",
         1: "~noActivators & noRepressors",
         2: "allActivators & ~allRepressors",
         3: "(noRepressors & ~noActivators) | (~allRepressors & allActivators)",
         4: "allActivators",
     5: "allActivators | (noRepressors & ~noActivators)",
     6: "~noActivators & ~allRepressors",
     7: "(~noActivators & ~allRepressors) | allActivators",
     8: "~noActivators",
     9: "noRepressors",
     10: "noRepressors | (~allRepressors & allActivators)",
     11: "noRepressors | (~noActivators & ~allRepressors)",
     12: "~allRepressors",
     13: "noRepressors | allActivators",
     14: "(noRepressors | allActivators) | (~allRepressors & ~noActivators)",
     15: "~allRepressors | allActivators",
     16: "noRepressors | ~noActivators",
     17: "~allRepressors | ~noActivators",
     18: "#activators > #repressors | #activators = #repressors & Q(gene)",
         19: "#activators > #repressors"
}
 
get_variables = lambda ls : ["S" + str(i) for i in ls]
get_valuation = lambda ls, q : [str(simplify(extract1BV(q, k))) for k in ls]
get_value = lambda ls, vs : [ls[i] + "=" + vs[i] for i in range(len(vs))]
 
## Returns two states q0, q1 of size size such as q0[idx] != q1[idx]
## and for all i != idx, q0[i] = q1[i]
def ensure_only_1_different(size, idx):
    q0 = generate_state(size)
    keep = idx
    q1 = simplify((q0 & intList2BV(filter(lambda x : x != keep, range(size)), 
        size)) | intList2BV([keep] if (simplify(extract1BV(q0, keep) == 0)) else [], size)) 
    print("for all i != idx, q0[i] = q1[i] ?")
    print(str(simplify(equalExcept1BV(q0, q1, idx))) + " should be True")
    print("q0[idx] != q1[idx] ?")
    print(str(not simplify(testValueVec(q0, q1, idx))) + " should be True")
    return([q0, q1])
 
#######################
## SHORTCUTS.py test ##
#######################
 
def test_shortcuts():
    from shortcuts import notInducible, neg, notRepressible, castBVinCond, nxor, wGRF, noActivatorsGRF, allActivatorsGRF, noRepressorsGRF, allRepressorsGRF, listActivators, listRepressors, prepreCompute, preCompute, diGRF, diCompute, rule18, rule19
    print("--------- START TEST")
    idx = 0
    size = 4
    C = get_variables(range(size))
    print(C)
    prelistQ = generate_bv(size)
    prelistQ = [3]  
    q = generate_state(size, prelist=prelistQ)
    prelistA = filter(lambda x : x!=idx, generate_bv(size))
    prelistA = [3]
    prelistR = filter(lambda x : not (x in prelistA + [idx]), generate_bv(size))
    activators = generate_regulators(idx, size, prelist=prelistA)
    if (not prelistR):
        repressors = buildZERO(size)
    else:
        repressors = generate_regulators(idx, size, prelist=prelistR)
    print("* State list: ")
    print(get_variables(prelistQ))
    print(">> State q: " + getBinaryDec(q, size))
    print("* Activators list: ")
    print(get_variables(prelistA))
    print(">> Activators: " + getBinaryDec(activators, size))
    print("* Repressors list: ")
    print(get_variables(prelistR))
    print(">> Repressors: " + getBinaryDec(repressors, size))
    [aA, aR, nA, nR, naA, naR, nnA, nnR, gc] = preCompute(q, prepreCompute(activators, repressors), test=True)
    inducibleRegulationGQ = neg(neg(notInducible(activators)) & notRepressible(repressors)) | nnA
    repressibleRegulationGQ = neg(notRepressible(repressors)) & notInducible(activators) & nR
    print("-----------------------")
    print("For gene " + C[idx] + ":")
    print("$ notInducible")
    print(simplify(notInducible(activators) == buildONE(size)))
    print(getBinaryDec(notInducible(activators), size))
    print("$ notRepressible")
    print(simplify(notRepressible(repressors) == buildONE(size)))
    print(getBinaryDec(notRepressible(repressors), size))
    print("$ repressibleRegulation")
    print(simplify(repressibleRegulationGQ == buildONE(size)))
    print(getBinaryDec(repressibleRegulationGQ, size))
    print("$ neg(repressibleRegulation)")
    print(simplify(neg(repressibleRegulationGQ) == buildONE(size)))
    print(getBinaryDec(neg(repressibleRegulationGQ), size))
    print("$ inducibleRegulation")
    print(simplify(inducibleRegulationGQ == buildONE(size)))
    print(getBinaryDec(inducibleRegulationGQ, size))
    print("$ neg(inducibleRegulation)")
    print(simplify(neg(inducibleRegulationGQ) == buildONE(size)))
    print(getBinaryDec(neg(inducibleRegulationGQ), size))
    print("$ noActivators")
    print(simplify(nA == buildONE(size)))
    print(getBinaryDec(nA, size))
    print("$ neg(noActivators)")
    print(simplify(nnA == buildONE(size)))
    print(getBinaryDec(nnA, size))
    print("$ neg(notInducible)")
    print(simplify(neg(notInducible(activators)) == buildONE(size)))
    print(getBinaryDec(neg(notInducible(activators)), size))
    print("$ allActivators")
    print(simplify(aA == buildONE(size)))
    print(getBinaryDec(aA, size))
    print("$ neg(allActivators)")
    print(simplify(naA == buildONE(size)))
    print(getBinaryDec(naA, size))
    print("$ noRepressors")
    print(simplify(nR == buildONE(size)))
    print(getBinaryDec(nR, size))
    print("$ neg(noRepressors)")
    print(simplify(nnR == buildONE(size)))
    print(getBinaryDec(nnR, size))
    print("$ neg(notRepressible)")
    print(simplify(neg(notRepressible(repressors)) == buildONE(size)))
    print(getBinaryDec(neg(notRepressible(repressors)), size))
    print("$ allRepressors")
    print(simplify(aR == buildONE(size)))
    print(getBinaryDec(aR, size))
    print("$ neg(allRepressors)")
    print(simplify(naR == buildONE(size)))
    print(getBinaryDec(naR, size))
    vA = get_valuation(prelistA, q)
    vR = get_valuation(prelistR, q)
    vC = get_valuation([idx], q)
    vvA = get_value(get_variables(prelistA), vA)
    vvR = get_value(get_variables(prelistR), vR)
    vvC = get_value(get_variables([idx]), vC)
    for i in range(20):
        print("_________________________________\n")
        print(">>> RULE #I = %d" % i)
        print("RULE: " + printRule.get(i))
        print(diGRF(0).get(i)([C[k] for k in prelistA], [C[k] for k in prelistR], C[idx]))
        if (i < 18):
            print(diGRF(0).get(i)(vA, vR, vC))
            d = diCompute.get(i)(aA, aR, nA, nR, naA, naR, nnA, nnR, gc)
            print("Result: " + str(simplify(d)))
            print("Binary form: " + getBinaryDec(d, size))
            print("== True (a.k.a. " + str(simplify(buildONE(size))) + " = " + getBinaryDec(buildONE(size), size) + ")?")
            print(simplify(d == buildONE(size)))
        else:
            print([vvA, vvR, vvC])
            print(simplify(ifthenelse(i==18, rule18(q, activators, repressors, idx), 
                rule19(q, activators, repressors, idx))))
    print("--------- END TEST")
 
#######################
## UTILS.py test     ##
#######################
 
def test_utils():
    print("------ START TEST")
    print(">>> Test getBinaryDec:")
    print(getBinaryDec(BitVecVal(3,3), 3) + " == 011")
    print(getBinaryDec(BitVecVal(3,3), 5) + " == 00011")
    print(str(getBinaryDec(BitVecVal(3,3))) + " == 3")
    print(">>> Test bv:")
    print(getBinaryDec(bv(0, 3), 3) + " == 001")
    print(getBinaryDec(bv(1, 3), 3) + " == 010")
    print(getBinaryDec(bv(2, 3), 3) + " == 100")
    print(getBinaryDec(bv(3, 3), 3) + " == 000")
    print(">>> Test buildZERO:")
    print(getBinaryDec(buildZERO(3), 3) + " == 000")
    print(getBinaryDec(buildZERO(6), 6) + " == 000000")
    print(">>> Test buildONE:")
    print(getBinaryDec(buildONE(3), 3) + " == 111")
    print(getBinaryDec(buildONE(6), 6) + " == 111111")
    print(">>> Test trim1BV:")
    print(">> First try:")
    print(getBinaryDec(BitVecVal(15, 5), 5) + " == 01111")
    print(getBinaryDec(trim1BV(BitVecVal(15, 5), 0), 4) + " == 0111")
    print(getBinaryDec(trim1BV(BitVecVal(15, 5), 1), 4) + " == 0111")
    print(getBinaryDec(trim1BV(BitVecVal(15, 5), 2), 4) + " == 0111")
    print(getBinaryDec(trim1BV(BitVecVal(15, 5), 3), 4) + " == 0111")
    print(getBinaryDec(trim1BV(BitVecVal(15, 5), 4), 4) + " == 1111")
    print(">> Second try:")
    print(getBinaryDec(BitVecVal(12, 5), 5) + " == 01100")
    print(getBinaryDec(trim1BV(BitVecVal(12, 5), 0), 4) + " == 0110")
    print(getBinaryDec(trim1BV(BitVecVal(12, 5), 1), 4) + " == 0110")
    print(getBinaryDec(trim1BV(BitVecVal(12, 5), 2), 4) + " == 0100")
    print(getBinaryDec(trim1BV(BitVecVal(12, 5), 3), 4) + " == 0100")
    print(getBinaryDec(trim1BV(BitVecVal(12, 5), 4), 4) + " == 1100")
    print(">>> Test extract1BV:")
    print(str(simplify(extract1BV(BitVecVal(3, 3), 0))) + " == 1")
    print(str(simplify(extract1BV(BitVecVal(3, 3), 1))) + " == 1")
    print(str(simplify(extract1BV(BitVecVal(3, 3), 2))) + " == 0")
    print(">>> Test intList2BV:")
    print(getBinaryDec(intList2BV([2, 3], 4), 4) + " == 1100")
    print(getBinaryDec(intList2BV([0], 4), 4) + " == 0001")
    print(getBinaryDec(intList2BV([], 4), 4) + " == 0000")
    print(">>> Test testValue:")
    print(str(simplify(testValue(BitVecVal(3, 3), 0, true))) + " == True")
    print(str(simplify(testValue(BitVecVal(3, 3), 0, false))) + " == False")
    print(">>> Test equalExcept1BV:")
    print(getBinaryDec(BitVecVal(3, 3), 3))
    print(getBinaryDec(BitVecVal(1, 3), 3))
    print("Coordinate 0: ")
    print(str(simplify(equalExcept1BV(BitVecVal(3, 3), BitVecVal(1, 3), 0))) + " == False")
    print("Coordinate 1: ")
    print(str(simplify(equalExcept1BV(BitVecVal(3, 3), BitVecVal(1, 3), 1))) + " == True")
    print("Coordinate 2: ")
    print(str(simplify(equalExcept1BV(BitVecVal(3, 3), BitVecVal(1, 3), 2))) + " == False")
    print(">>> Test testValueVec:")
    print("(1): " + getBinaryDec(BitVecVal(3, 3), 3))
    print("(2): " + getBinaryDec(BitVecVal(1, 3), 3))
    print("Coordinate 0: ")
    print(str(simplify(testValueVec(BitVecVal(3, 3), BitVecVal(1, 3), 0))) + " == True")
    print("Coordinate 1: ")
    print(str(simplify(testValueVec(BitVecVal(3, 3), BitVecVal(1, 3), 1))) + " == False")
    print("Coordinate 2: ")
    print(str(simplify(testValueVec(BitVecVal(3, 3), BitVecVal(1, 3), 2))) + " == True")
    print("Coordinate 0 (1) and 2 (2): ")
    print(str(simplify(testValueVec(BitVecVal(3, 3), BitVecVal(1, 3), 0, c1=2))) + " == False")
    print("Coordinate 1 (1) and 0 (2): ")
    print(str(simplify(testValueVec(BitVecVal(3, 3), BitVecVal(1, 3), 1, c1=0))) + " == True")
    print(">>> Test evalVec:")
    print("Dictionary: " + str({0:1, 1:2, 2:0}))
    print("Two bit-vector variables x > 0 and y > 0 of size 3...")
    x = BitVec("x", 3)
    y = BitVec("y", 3)
    s = Solver()
    s.add(x > buildZERO(3))
    s.add(y > buildZERO(3))
    print("Conditions x[1] == y[0], x[2] == y[1], x[0] == y[2]")
    s = evalVec(s, x, y, {0:1, 1:2, 2:0}, default)
    ## To see the conditions in native Z3 format       ##
    #print(s.sexpr())
    print(s.check())
    print(" == sat")
    M = s.model()
    print("x = " + getBinaryDec(M[x], 3))
    print("y = " + getBinaryDec(M[y], 3))
    print(">>> Test sub:")
    print(str(simplify(sub(BitVecVal(7,3)))) + " == 3")
    print(str(simplify(sub(BitVecVal(0,3)))) + " == 0")
    print(">>> Test rev:")
    print(str(rev(range(3))) + " == [2, 1, 0]")
    print(">>> Test getIt:")
    print(getIt("bl#ablal*ba", "#", "*") + " == ablal")
    print("------ END TEST")
 
###########################
## GRN_INFERENCE.py test ##
###########################
 
def testTransition(size, Rs, res, a, r, q0, q1, typeT):
    from grn_inference import preCompute, prepreCompute
    from shortcuts import rule18, rule19
    typeT += "hronous"
    verbose = False
    s = Solver()
    ## Checks if at least one index satisfies equalExcept1BV(q0, q1, g) ##
    condAsync = False
    s = regulation_condition(s, Rs, [range(20)]*size, size)
    for g in range(size):
        pcg = preCompute(q0, prepreCompute(a[g], r[g]))
        auxtrr = lambda i : aux_transition(i, g, q1, res, size, pcg)
        auxtr17 = lambda i : ifthenelse(i==18, rule18(q0, a[g], r[g], g), rule19(q0, a[g], r[g], g))
        for i in range(20):
            auxtr = ifthenelse(i<18, auxtrr, auxtr17)
            conclusion = str(simplify(auxtr(i)))
            if (typeT == "asynchronous"):
                hypothesis = str(simplify(testRS(Rs[g], i)))
                hypothesis2 = str(simplify(equalExcept1BV(q0, q1, g)))
                verboseIt(hypothesis2 + " => " + "(" + hypothesis 
                        + " => " + conclusion + ")", verbose)
                if (hypothesis2 == "True" and conclusion == "False"):
                    verboseIt(">> Add this condition to solver:", verbose)
                    verboseIt("NOT(" + hypothesis + ")", verbose)
                cond = transition_condition_async(Rs, g, i, q0, q1, auxtr)
            elif (typeT == "synchronous"):
                hypothesis = str(simplify(testRS(Rs[g], i)))
                verboseIt(hypothesis + " => " + conclusion, verbose)
                cond = transition_condition_sync(Rs, g, i, auxtr)
            else:
                print("ERROR: Wrong transition type.")
                return(None)
            verboseIt("Using rule #" + str(i) + ": " + str(simplify(cond)), verbose)
            s.add(cond)
        if (typeT == "asynchronous"):
            condAsync = Or(condAsync, equalExcept1BV(q0, q1, g))
    if (typeT == "asynchronous"):
        s.add(condAsync)
    print("* Initial state:")
    print(">> " + getBinaryDec(q0, size))
    print("* Final state:")
    print(">> " + getBinaryDec(q1, size))
    print("* Regulators:")
    print(">> Activators: " + str([getBinaryDec(aa, size) for aa in a]))
    print(">> Repressors: " + str([getBinaryDec(rr, size) for rr in r]))
    print(s.check())
    try:
        M = s.model()
        print(M)
    except:
        print("No model")
 
def test_grn_inference():
    from grn_inference import buildInteractionDict, preCompute, prepreCompute
    from shortcuts import diCompute, rule18, rule19
    print("------ START TEST")
    size = 4
    idx = 1
    print("Test buildInteractionDict:")
    C = get_variables(range(size))
    gene = C[idx]
    Iopt = [["S1", "S2", "+"], ["S0", "S1", "+"], ["S2", "S3", "-"], ["S2", "S1", "+"], ["S1", "S1", "-"]]
    regIntActivators, regIntRepressors = buildInteractionDict(Iopt, gene, C, opt=True)
    print("---")
    print("Optional interaction set:")
    print(Iopt)
    print(range(len(Iopt)))
    print("Gene list:")
    print(C)
    print(range(size))
    print("Activators: ")
    print(regIntActivators)
    print("Repressors: ")
    print(regIntRepressors)
    Idef = [["S1", "S0", "-"], ["S3", "S1", "+"]]
    regIntActivators, regIntRepressors = buildInteractionDict(Idef, gene, C, regIntActivators=regIntActivators, regIntRepressors=regIntRepressors)
    print("---")
    print("Definite interaction set:")
    print(Idef)
    print(range(len(Idef)))
    print("Gene list:")
    print(C)
    print(range(size))
    print("Activators: ")
    print(regIntActivators)
    print("Repressors: ")
    print(regIntRepressors)
    print("---")
    print("Test testRS:")
    print(str(simplify(testRS(BitVecVal(3, 5), 3))) + " == True")
    print(str(simplify(testRS(BitVecVal(0, 5), 3))) + " == False")
    print("Test getCState:")
    print(str(simplify(getCState(BitVecVal(3, 20), 1, 20))) + " == 2")
    print(str(simplify(getCState(BitVecVal(3, 20), 2, 20))) + " == 0")
    print("Test aux_transition:")
    q = generate_state(size)
    print("* State:")
    print(">> " + getBinaryDec(q, size))
    prelistA = generate_bv(size)
    a = generate_regulators(idx, size, prelist=prelistA)
    print("* Activators:")
    print(">> " + getBinaryDec(a, size))
    r = generate_regulators(idx, size, prohibitedIdx=prelistA)
    print("* Repressors:")
    print(">> " + getBinaryDec(r, size))
    [aA, aR, nA, nR, naA, naR, nnA, nnR, gc] = preCompute(q, prepreCompute(a, r))
    res = lambda x : x
    for i in range(20):
        if (i < 18):
            rss = int(bool(simplify(diCompute.get(i)(aA, aR, nA, nR, naA, naR, nnA, nnR, gc))))
            print("Rule #" + str(i) + ": " + str(rss) + " == "
                + str(simplify(extract1BV(If(diCompute.get(i)(aA, aR, nA, nR, naA, naR, nnA, nnR, gc), 
                bv(idx, size), buildZERO(size)), idx))))
        if (i == 18):
            rss = int(bool(simplify(rule18(q, a, r, idx))))
            print("Rule #18: " + str(rss) + " == "
                + str(simplify(extract1BV(If(rule18(q, a, r, idx), 
                bv(idx, size), buildZERO(size)), idx))))
        else:
            rss = int(bool(simplify(rule19(q, a, r, idx))))
            print("Rule #19: " + str(rss) + " == "
                + str(simplify(extract1BV(If(rule19(q, a, r, idx), 
                bv(idx, size), buildZERO(size)), idx))))
    print("Test conditions:")
    Rs = [BitVec("Rs_%d" % i, 5) for i in range(size)]
    print("_______ transition_condition_sync:")
    ## To generate randomly states and regulators:     ##
    #q0 = generate_state(size)
    #q1 = generate_state(size)
    #prelists = [generate_bv(size) for i in range(size)]
    #a = [generate_regulators(idx, size, prelist=prelists[i]) for i in range(size)]
    #r = [generate_regulators(idx, size, prohibitedIdx=prelists[i]) for i in range(size)]
    print(">>>> First try:")
    q0 = intList2BV([2], size)
    q1 = intList2BV([0, 2, 3], size)
    a = custom_int2regulators(['0000', '0000', '1100', '0001'], size)
    r = custom_int2regulators(['0100', '0001', '0000', '1100'], size)   
    testTransition(size, Rs, res, a, r, q0, q1, typeT="sync")
    print(" ==> unsat (fails because no regulation function can satisfy the transition condition)")
    print(">>>> Second try:")
    q0 = intList2BV(range(1, 4), size)
    q1 = intList2BV([0, 2, 3], size)
    a = [intList2BV([0, 2], size), intList2BV([0, 2], size), intList2BV([0, 3], size), intList2BV([0, 2], size)]
    r = [intList2BV([], size), intList2BV([3], size), intList2BV([2], size), intList2BV([3], size)]
    testTransition(size, Rs, res, a, r, q0, q1, typeT="sync")
    print(" ==> sat")
    print("_______ transition_condition_async:")
    print(">>>> First try:")
    q0 = intList2BV(range(4), size)
    q1 = intList2BV([2, 3], size)
    a = [intList2BV([0, 2], size), intList2BV([0, 2], size), intList2BV([0, 3], size), intList2BV([0, 2], size)]
    r = [intList2BV([], size), intList2BV([3], size), intList2BV([2], size), intList2BV([3], size)]
    testTransition(size, Rs, res, a, r, q0, q1, typeT="async")
    print(" ==> unsat (fails because q1 cannot be reached in one async. step from q0)")
    print(">>>> Second try:")
    q0 = bv(2, size)
    q1 = intList2BV([0, 2], size)
    print("q0 = " + getBinaryDec(q0, size) + ", q1 = " + getBinaryDec(q1, size))
    print("(q1 can maybe be reached in one async. step from q0)") 
    a = [bv(2, size), intList2BV([0, 3], size), intList2BV([0, 2, 3], size), intList2BV([0, 3], size)]
    r = [bv(3, size)] + [buildZERO(size)]*3
    testTransition(size, Rs, res, a, r, q0, q1, typeT="async")
    ## To see the conditions on regulation functions   ##
    #pcg = preCompute(q0, prepreCompute(a[idx], r[idx]))
    #auxtrr = lambda i : aux_transition(i, idx, q1, res, size, pcg)
    #auxtr17 = lambda i : ifthenelse(i==18, rule18(q0, a[idx], r[idx], idx), rule19(q0, a[idx], r[idx], idx))
    #for i in range(20):
    #   print("Rule #" + str(i) + ": q1[0] == r_" + str(i) + "(q0) ?")
    #   auxtr = ifthenelse(i<18, auxtrr, auxtr17)
    #   print(str(simplify(auxtr(i))))
    print(" ==> sat (and actually, in this case, any regulation function except for the threshold functions for gene S0 can be used)")
    print(">>>> Third try:")
    [q0, q1] = [bv(0, size), intList2BV([0, 1], size)]
    a = [buildZERO(size), intList2BV([0, 3], size), bv(3, size), intList2BV([0, 3], size)]
    r = [bv(3, size), buildZERO(size), bv(2, size), bv(2, size)]
    testTransition(size, Rs, res, a, r, q0, q1, typeT="async")
    ## To see the conditions on regulation functions   ##
    #pcg = preCompute(q0, prepreCompute(a[idx], r[idx]))
    #auxtrr = lambda i : aux_transition(i, idx, q1, res, size, pcg)
    #auxtr17 = lambda i : ifthenelse(i==18, rule18(q0, a[idx], r[idx], idx), rule19(q0, a[idx], r[idx], idx))
    #for i in range(20):
    #   print("Rule #" + str(i) + ": q1[1] == r_" + str(i) + "(q1) ?")
    #   auxtr = ifthenelse(i<18, auxtrr, auxtr17)
    #   print(str(simplify(auxtr(i))))
    print(" ==> sat (and any regulation function EXCEPT FOR rules 0, 2 and 4 can be used for gene S1)")
    return(None)
    print("------ END TEST")
 
##########################
## LAUNCH_MODEL.py test ##
##########################
 
def test_launch_model():
    from models import readREINfile
    from launch_model import model2igraph, launch_model
    from grn_solver import grn_solver
    from grn_inference import getState
    print("------- START TEST")
    print("Test model2igraph:")
    print(">>> With non-coloured nodes:")
    C = ["Gene_" + str(i) for i in range(3)]
    Idef = [[C[1], C[2], "+"], [C[0], C[2], "-"]]
    Iopt = [[C[0], C[0], "+"], [C[1], C[0], "+"], [C[2], C[1], "-"]]
    resList = [[["Is", "011"]] + [["R", 10]]*3]
    P = [[]*3]
    model2igraph(0, resList, C, Idef, Iopt, P, plotIt=True, verbose=True)
    print(">>> With coloured nodes:")
    P = [["-", "+"], ["-"], ["+"]]
    model2igraph(0, resList, C, Idef, Iopt, P, plotIt=True, verbose=True)
    print("Test launch_model:")
    print("On toy model:")
    [C, CRM, length, Idef, Iopt, R, E, typeT, solmax, KO, FE, uniqueness, limreg, P, Fixpoint] = readREINfile(model="toy/model_expanded.net", experiments="toy/observations.spec")
    res = grn_solver(C, CRM, length, Idef, Iopt, R, E, typeT, 1, KO, FE, uniqueness, limreg, P, Fixpoint)
    try:
        [resList, _, _] = res
    except:
        resList = None
    if (resList):
        ## Random generation of state                      ##
        #q0 = [randint(0, 1) for i in range(len(C))]
        q0 = ["0", "1", "1", "0", "1"]
        print("Test on first model found:")
        modelID = 0
	printStates([[">>> Initial state :", [["q0", q0]]]], 0, C)
        sstep = 6
        print(">>> Number of steps: " + str(sstep))
        print(">>> Not steady state:")
        states = launch_model(modelID, C, CRM, resList, Idef, Iopt, R, q0, 
		sstep, typeT, KO, FE, P)
	for i in range(len(states)):
		printStates(states, i, C)
        print("\n")
        print(">>> Steady state:")
        states = launch_model(modelID, C, CRM, resList, Idef, Iopt, R, q0, 
		sstep, typeT, KO, FE, 
                P, steadyStates=True)
	for i in range(len(states)):
		printStates(states, i, C)
    print("------- END TEST")
 
##########################
## CALL                 ##
##########################
 
tests = [test_shortcuts, test_utils, test_grn_inference, test_launch_model]
names = ["shortcuts", "utils", "grn_inference", "launch_model"]
i = 0
lenargvC = len(sys.argv) == 2
 
def printRunSyntaxError(c):
    if (not c):
        print("MSG: If you wanted to run the solver, then you did not use the correct syntax.")
        print("MSG: Correct syntax is \'run [--simplify] [--visualize] model experiments\'.")
        return(True)
    return(False)
 
done = False
while (i < 4 and lenargvC):
	if (names[i] == sys.argv[1]):
		tests[i]()
		done = True
		break
	i += 1
if (i == 4 and lenargvC and not done):
	print("MSG: If you wanted to run a test on a file, maybe you mistyped the filename.")
	print("MSG: You may want to test one of these: " + str(names))
