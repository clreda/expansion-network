#################
## Libraries   ##
#################

library(QCA)
options(warn=-1)

#################
## Utils       ##
#################

maxnvar = 8

tail <- function(ls) return(ls[2:length(ls)])

rm <- function(ls, patt) return(ls[ls != patt])

#' Change logical operators in string
#'
#' @param y character string = boolean formulae
#' @return y character string with replacements 
#'
#' @export
changeOperators <- function(y) {
	AND="*"
	OR="+"
	NOT="!"
	return(gsub("AND", AND, gsub("NOT", NOT, gsub("OR", OR, y))))
}

#' Get causal condition names for logical formulae
#'
#' @param y character string = boolean formulae
#' @return cond character string vector for conditions
#'
#' @export
getConditions <- function(y) {	
	AND="*"
	OR="+"
	NOT="!"
	f <- function(op) return(paste0("[", op, "]"))
	AND <- f(AND)
	OR <- f(OR)
	NOT <- f(NOT)
	y <- paste0(tail(strsplit(y, "= ")[[1]]), collapse="= ")
	y <- gsub("[)]", "", gsub("[(]", "", gsub(AND, "", 
			gsub(NOT, "", gsub(OR, "", y)))))
	y <- gsub("[ ]+", " ", y)
	cond <- unique(strsplit(y, "[ ]")[[1]])
	return(cond)
}

#' Get boolean function corresponding to logical formulae
#'
#' @param y character string = boolean formulae
#' @return f associated boolean function: evaluation of formulae
#'               respect to given variable assignment
#'
#' @export
getFormulae <- function(y) {
	y <- strsplit(gsub("[ ]", "", y), "=")[[1]][2]
	f <- function(args) {
		argsLength <- unlist(sapply(names(args), function(e) length(strsplit(e, "")[[1]])))
		argsOrder <- sort(argsLength, decreasing=TRUE, index.return=T)$ix
		args <- args[argsOrder]
		for (i in 1:length(args)) y <- gsub(names(args)[i], 
			as.character(args[i]), y)
		y <- gsub("0", "FALSE", y)
		y <- gsub("1", "TRUE", y)
		y <- gsub("True","TRUE", gsub("False","FALSE", y))
		y <- gsub("![(]TRUE[)]","(!TRUE)", gsub("![(]FALSE[)]","(!FALSE)", y))
		res <- as.numeric(as.logical(eval(parse(text=y))))
		return(res)
	}
	return(f)
}

#' Get outcome names for logical formulae
#'
#' @param y character string = boolean formulae
#' @return out outcome for this logical formulae
#'
#' @export
getOutcomes <- function(y) return(strsplit(y, "' = ")[[1]][1])

#' Return named object
#'
#' @param obj nameable R object
#' @param cc character string vector of names of same length
#' @return obj named with cc
#'
#' @export
nameIt <- function(obj, cc, make.list=T) {
	obj <- if (make.list) sapply(obj, list) else obj
	names(obj) <- cc
	return(obj)	
}

#' Return named object by itself
#'
#' @param obj character string vector
#' @return obj named with itself: names(obj) = obj
#'
#' @export
nameItself <- function(obj) return(nameIt(obj, obj))

#################
## BFL object  ##
#################

#' Return object Boolean Formulae List
#'
#' @param conditions named (by outcome) list of character string vectors
#' @param outcome list of outcome variable names
#' @param formulas named (by outcome) list of boolean functions
#' @return generics a Boolean Formulae List object
#'
#' @export
createBooleanFormulaeList <- function(conditions, outcomes, formula, x) {
	formulas <- formula
	idx <- which(sapply(conditions, function(e) length(e) > maxnvar))
	if (length(idx) > 0) {
		formulas[idx] <- sapply(x[idx], function(e) gsub("!", "~", strsplit(e, " = ")[[1]][2]))
	}
	generics <- list(conditions=conditions,
			outcomes=outcomes, 
			formula=formulas)
	return(generics)
}

#' Get Boolean Formulae from object Boolean Formulae List
#'
#' @param generics Boolean Formulae List object
#' @param out outcome associated to output Boolean Formulae
#' @return generic conditions, outcome name, and function associated
#'                     with out
#'
#' @export
getBooleanFormulaeFromBFL <- function(generics, out) {
	return(list(conditions=generics$conditions[[out]], 
			outcome=out, 
			formula=generics$formula[[out]]))
}

###################
## Input/Output  ##
###################

#' Read RE:IN formatted resulting model
#' 
#' @param filename name of file that contains the model returned by RE:IN
#' @return generics a Boolean Formulae List object: a Boolean formulae is a list 
#'                  with conditions, outcome, formulae (Boolean function) 
#'
#' @export
readModel <- function(filename) {
	## Scan the lines of file                            ##
	x <- scan(filename, what="", sep="\n")
	## Get lines into list                               ##
	x <- sapply(x, list)
	## Make outcome variable names be names of each line ##
	names(x) <- sapply(x, getOutcomes)
	## Make the logical operators evaluable by R         ##
	x <- sapply(x, changeOperators)
	## Get outcome variable names                        ##
	outcomes <- names(x)
	## Get condition variable names                      ##
	conditions <- sapply(x, getConditions)
	conditionsIDX <- sapply(conditions, function(cc) !(cc %in% c("", "True", "False")))
	conditions <- sapply(outcomes, function(g) conditions[[g]][conditionsIDX[[g]]])
	## Get Boolean functions associated with the Boolean ##
	## formulas in the model                             ##
	formula <- sapply(x, getFormulae)
	## Create the associated Boolean Formula List object ##
	generics <- createBooleanFormulaeList(conditions, outcomes, formula, x)
	return(generics)
}

formatSolution_aux <- function(y, node) {
	return(paste0(sapply(y[[node]]$solution,
				function(solution) 
					paste0(
					"GRF(", node, ") = ",
						paste0(solution, 
						collapse=" + "))
				),
			collapse="\n"))
}

formatSolution <- function(y, node) {
	if (!is.character(y[[node]])) return(formatSolution_aux(y, node))
	else return(paste0("GRF(", node, ") = ", sub(" <=> OUTCOME", "", y[[node]])))
}

#' Write the simplified model
#' 
#' @param y result of Boolean simplification
#' @param title name of the file to write
#'
#' @export
writeModel <- function(y, title) {
	str=""
	nodes <- names(y)
	for (node in nodes) {
		str <- paste0(str, "\n", formatSolution(y, node))
	}
	write(sub("\n", "", str), file=title)
	return(NULL)
}

###########################
## Truth Table creation  ##
###########################

#' Enumerate all possibilities for n Boolean variables
#' 
#' @param n number of Boolean variables
#' @return res matrix of size nrow x n, 
#'                nrow = number of possibilities (2^n) and ncol = n
#'                such as each line is a unique vector of values for each variable 
#'
#' @export
enumerateGrayCode <- function(n) return(do.call(expand.grid, rep(list(0:1), n)))

## not == function(x) return(as.numeric(!as.logical(x))) ##
not <- function(x) return(as.numeric(as.logical(x-1)))

#' Get the truth table from a Boolean Formulae object
#' 
#' @param generic Boolean Formulae object
#' @param neg.out should the outcome be negated?
#' @return tt the truth table associated to the Boolean Formulae object 
#'
#' @export
generic2truthTable <- function(generic, neg.out) {
	## Get attributes of generic       ##
	conditions <- generic$conditions
	outcome <- generic$outcome
	f <- generic$formula
	## Get the matrix of possibilities ##
	gcode <- enumerateGrayCode(length(conditions))
	## Evaluate the Boolean function   ##
	## in each combination             ##
	outcome <- apply(gcode, 1, function(x) {
			names(x) <- conditions;
			res <- if (neg.out) not(f(x)) else f(x);
			res
			})
	## Get the truth table             ##
	res <- as.data.frame(cbind(gcode, outcome))
	colnames(res) <- c(conditions, "OUTCOME")
	## QCA wrapper                     ##
	return(truthTable(res, outcome="OUTCOME", conditions=conditions))
}

###########################################
## Quine-McCluskey boolean minimization  ##
###########################################

retOneCond <- function(b) return(paste0(as.character(b)," <=> OUTCOME"))

#' Apply QCA enchanced Quine-McCluskey minimization to
#'               truth table
#' 
#' @param tt truth table
#' @return y QCA eqmcc result object
#'
#' @export
truthTable2simplified <- function(tt) {
	outcomeSum <- sum(as.numeric(tt$tt$OUT))
	if (outcomeSum == 0) return(retOneCond(FALSE))
	if (any(colnames(tt$tt) == "FALSE")) {
		## Remove cases where sth*FALSE <=> OUT ##
		zeroIdx <- tt$tt["FALSE"] == 0
		zeroOut <- tt$tt$OUT[zeroIdx]
		zeroFalse <- tt$tt["FALSE"][zeroIdx]
		if (all(zeroOut == zeroFalse)) {
			return(retOneCond(FALSE))
		}
	}
	idxCond <- 1:(ncol(tt$tt)-5)
	if (outcomeSum == length(tt$tt$OUT)) return(retOneCond(TRUE))
	return(eqmcc(tt, outcome="OUT", relation="sufnec", use.tilde=T, 
			conditions=colnames(tt$tt)[idxCond]))
}

###############
## Aux fcts  ##
###############

#' Get simplified boolean formulae
#' 
#' @param generics Boolean Formulae List object
#' @param out outcome variable name
#' @param neg.out should the outcome be negated?
#' @return res QCA result object, or result string
#'
#' @export
boolean_reducer_aux_aux <- function(generics, out, neg.out) {
	print(out)
	generic <- getBooleanFormulaeFromBFL(generics, out)
	conditions <- generic$conditions
	if (length(conditions) == 0) return(retOneCond("False"))
	if (length(conditions) == 1) {
		isActive <- nameIt(1, conditions)
		f <- generic$formula
		if (f(isActive) == 1) return(retOneCond(conditions))
		if (f(isActive) == 0) return(retOneCond(paste0("~", conditions)))		
		return(paste0("Error for ", out))
	}
	## Cases of threshold GRFs                ##
	if (any(c(">", "=") %in% conditions)) {
		trim <- function(b, e) return(rm(rm(conditions[b:e], "="), out))
		idx <- which(conditions == ">")
		act <- paste0(trim(1, idx-1), collapse=" ")
		rr <- paste0(trim(idx+1, length(conditions)), collapse=" ")
		rep <- ifelse(idx == length(conditions) || rr == "" || is.na(rr), act, rr)
		if ("=" %in% conditions) {
			return(retOneCond(paste0("(", act, " > ", rep, ") OR ((", act, " = ", rep, ") AND ", out, ")")))
		}
		else return(retOneCond(paste0("(", act, " > ", rep, ")")))
	}
	if (length(conditions) > maxnvar) return(retOneCond(generic$formula))
	res <- truthTable2simplified(generic2truthTable(generic, neg.out=neg.out))
	## Remove cases where sth + FALSE <=> OUT ##
	if (!is.character(res)) res$solution <- lapply(res$solution, 
			function(solution) solution[solution!="FALSE"])
	return(res)
}

#' Get simplified boolean formulae for each outcome variable
#' 
#' @param generics Boolean Formulae List object
#' @param neg.out should the outcome be negated?
#' @return res QCA result object, or result string list
#'
#' @export
boolean_reducer_aux <- function(generics, neg.out) {
        y <- sapply(generics$outcomes, 
		function(out) boolean_reducer_aux_aux(generics, out, neg.out))
	return(y)
}

###############
## Oneliner  ##
###############

#' Get simplified boolean formulae for each outcome variable in result RE:IN model
#' 
#' @param filename file where result model is stored
#' @param neg.out should the outcome be negated?
#' @param title name of the resulting file
#' @return res QCA result object, or result string list
#'
#' @export
boolean_reducer <- function(filename, neg.out=F, 
			title=paste0("simplified--", filename)) {
	generics <- readModel(filename)
	y <- boolean_reducer_aux(generics, neg.out)
	writeModel(y, title)
	return(y)
}
