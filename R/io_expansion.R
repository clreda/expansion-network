source("utils_expansion.R")

#############################
## GENERIC MODEL STRUCTURE ##
#############################

#' Create generic model structure (see specifications file)
#'
#' @param directives -named- list of RE:IN  directives in the following order: 
#'                  updates, length, uniqueness, limit, regulation
#' @param nodes vector of character strings i.e. gene names
#' @param conditions vector of size |nodes| of character strings (of type a,b or a..b)
#' @param perturbations vector of size |nodes| of character strings (-+ (KO/FE), + (FE), - (KO))
#' @param optionals -named- list of optional interactions i.e. 3-sized vectors: 
#'                   regulator in nodes, regulated gene in nodes, sign (positive, negative)
#' @param definites -named- list of definite interactions i.e. 3-sized vectors: 
#'                   regulator in nodes, regulated gene in nodes, sign (positive, negative)
#' @return generic model
#' 
#' @export
createModel <- function(directives, nodes, conditions, perturbations, optionals, definites) {
	return(list(directives=directives, 
		nodes=nodes, 
		conditions=conditions, 
		perturbations=perturbations, 
		definites=definites, 
		optionals=optionals))
}

#############################
## INPUT for MODELS        ##
#############################

#' Read RE:IN formatted model from file such as model.net
#' 
#' @param filename name of file, ending with .net, that contains the model
#' @return model generic model structure encoding for the model
#'
#' @export
readModel <- function(filename) {
	## Scan the lines of file                ##
	x <- scan(filename, what="", sep="\n")
	## Process directive lines               ##
	## (the first five lines)                ##
	x[1:5] <- sapply(x[1:5], function(y) gsub(";", "", sub("[[:space:]]", "_", y)))
	## Process node line                     ##
	## (the sixth line)                      ##
	x[6] <- paste("nodes", gsub(",[[:space:]]", ",", gsub(";", "", x[6])))
	## Process edge lines                    ##
	## (from the 7th line to the end of file)##
	x[7:length(x)] <- sapply(x[7:length(x)], 
		function(y)
			paste(
				paste(ifelse(grepl("optional", y), "optional", "definite"), 
					which(x == y), 
					sep=""), 
				sub("optional", "", gsub(";", "", y))
			)
	)
	## Get a named list of lines             ##
	## (1 line = 1 edge)                     ##
	y <- strsplit(x, "[[:space:]]+")
	names(y) <- sapply(y, `[[`, 1)
	y <- lapply(y, `[`, -1)
	## Get RE:IN directives                  ##
	directives <- grepl("directive", names(y))
	## Get optional interactions             ##
	optionals <- grepl("optional", names(y))
	## Get definite interactions             ##
	definites <- grepl("definite", names(y))
	## Get nodes                             ##
	nodes <- nameItself(sapply(y$nodes, function(x) strsplit(x, "[[]")[[1]][1]))
	## Get conditions for each node          ##
	conditions <- nameIt(sapply(y$nodes, function(x) strsplit(strsplit(x, "[(]")[[1]][2], "[)]")[[1]][1]),
			 nodes)
	## Get perturbations for each node       ##
	perturbations <- nameIt(sapply(y$nodes, function(x) strsplit(strsplit(x, "[[]")[[1]][2], "[]]")[[1]][1]),
			 nodes)
	## Create model                          ##
	model <- createModel(directives=y[directives], 
				nodes=nodes, 
				conditions=conditions, 
				perturbations=perturbations, 
				definites=y[definites], 
				optionals=y[optionals])
	return(model)
}

#############################
## OUTPUT for MODELS       ##
#############################

#' Write generic model into RE:IN formatted file
#' 
#' @param model generic model object
#' @param title name for the resulting model file
#' @param format "RE:IN"-like or "expanded"
#'
#' @export
writeModel <- function(model, title="newmodel.net", format="RE:IN") {
	str <- ""
	## Writes RE:IN directives                   ##
	dir <- c("updates", "length", "uniqueness", "limit", "regulation")
	dirstr <- paste0(
			sapply(1:length(dir), 
				function(i) {
					idx <- getIdxOfColumn(dir[i], model$dir)
					paste("directive", 
						dir[i], 
						paste0(model$dir[[idx]][[1]], ";"), sep=" ") 
					}
				), 
			collapse="\n")
	str <- paste0(str, dirstr)
	## Writes nodes                              ##
	## (Last white space after the list of nodes ##
	## is TREMENDOUSLY IMPORTANT)                ##
	nodestr <- paste0(paste0(sapply(model$nodes, function(nd) {
				## Associated index in perturbation list ##
				idxP <- getIdxOfData(nd, model$pert)
				## Associated index in condition list    ##
				idxC <- getIdxOfData(nd, model$cond)
				ascrm <- strsplit(nd, "_")[[1]]
				paste0(nd, 
					"[", model$pert[idxP], "]", 
					ifelse(format=="expanded", paste0("{", ifelse(is.na(ascrm[2]), "", ascrm[1]), "}"), ""),
					"(", gsub(",", ", ", model$cond[idxC]), ");"
					)
				}), 
				collapse=" "), 
			" ")
	str <- paste0(str, "\n", nodestr)
	## Writes optional interactions              ##
	## (if there exists one)                     ##
	if (length(model$opt) > 0) {
		optstr <- paste0(sapply(model$opt, 
				function(e) paste0(c(e, "optional;"), collapse="\t")), 
				collapse="\n")
		str <- paste0(str, "\n", optstr)
	}
	## Writes definite interactions              ##
	## (if there exists one)                     ##
	if (length(model$def) > 0) {
		defstr <- paste0(sapply(model$def, function(e) paste0(paste0(e, 
									collapse="\t"), ";")), 
								collapse="\n")
		str <- paste0(str, "\n", defstr)
	}
	## Writes into text file                     ##
	write(str, file=title)
}

