require(XML)
require(igraph)
source("io_expansion.R")

################
## Utils      ##
################

## Removes white spaces     ##
ws <- function(e) return(ifelse(length(gsub(" ", "", as.vector(e)))==0, "0", gsub(" ", "", as.vector(e))))
## Returns a named R object ##
nameIt <- function(e, nms) {
	names(e) <- nms
	return(e)
}
## Returns a named set, with##
## names: "header1", 	    ##
## "header2", etc.          ##
nameInt <- function(set, header) {
	if (!is.null(set)) {
		return(nameIt(set, sapply(1:length(set), function(i) paste0(header, i))))
	}
	return(set)
}
#' Returns the "tail" of a list 
#'
#' @param ls R list
#' @param from start index of the tail
#' @return t tail of the list. If from > length of the list, 
#'         returns a NULL element.
#'
#' @export
get_list_from <- function(ls, from) if (length(ls) < from) return(NULL) else return(ls[from:length(ls)])

####################
## SBML parsing   ##
####################

get_node <- function(e) return(e['id'])
get_regulators <- function(e) return(as.vector(sapply(e$listOfInputs, function(ee) as.vector(ee['qualitativeSpecies']))))
get_regulated <- function(e) return(as.vector(e$listOfOutputs$output['qualitativeSpecies']))
get_sign <- function(e) return(sapply(e$listOfInputs, function(ee) as.vector(ee['sign'])))
get_default_value <- function(e) return(as.vector(e$listOfFunctionTerms$defaultTerm$.attrs['resultLevel']))
## Recursively gets the values of each regulator in the current function term ##
get_1_grf <- function(eftm) {
	## Get operator         ##
	operator <- names(eftm)[1]
	## Either "eq" or "geq" ##
	if (grepl("eq", operator)) {
		op <- ifelse(operator=="eq", "=", ifelse(operator=="geq", ">=", "<="))
		return(paste0(ws(eftm['ci']), op, ws(eftm$cn['text'])))
	}
	## Either "and" or "or" ##
	else {
		return(paste0("(", 
			paste0(sapply(get_list_from(eftm, 2), get_1_grf), collapse=paste0(") ", operator, " (")),
			")")
			)
	}
}
## Summary of the default value and all regulatory functions for a given term ##
## associated with a gene                                                     ##
get_grf <- function(e) {
	default = paste0("default-value:", get_default_value(e))
	fcts = sapply(get_list_from(e$listOfFunctionTerms, 2), function(eft) {
				paste0("function: ", 
				get_1_grf(eft$math$apply), 
				" -> ",
				ws(eft$.attrs['resultLevel'])
				)
			})
	if (length(fcts) == 0) fcts <- ""
	else fcts <- paste0(c("", fcts), collapse=", ")
	return(paste0(default, fcts))
}
## Builds an igraph object                                                    ##
get_graph <- function(edges, signs, grfs) {
	g <- make_graph(edges=edges, directed=TRUE) %>% set_edge_attr("sign", value=signs)
	nodes <- get.vertex.attribute(g)$name
	## Accounts for non-output nodes ##
	idx <- which(nodes %in% names(grfs))
	g <- g %>% set_vertex_attr("grf", index = idx, grfs[nodes[idx]])
	return(g)
}
## Get node-pairwise interactions                                             ##
get_interactions <- function(edges, signs) {
	int <- lapply(seq(1, length(edges), 2), function(i) c(edges[i], edges[i+1], ""))
	left <- NULL
	opt <- NULL
	idxopt <- NULL
	for (i in 1:length(int)) {
		if (signs[i] == "dual") {
			left <- c(left, list(c(int[[i]][1], int[[i]][2], "negative")))
			int[[i]][3] <- "positive"
		}
		else {
			if (signs[i] == "unknown") {
				opt <- c(opt, list(c(int[[i]][1], int[[i]][2], "positive")))
				opt <- c(opt, list(c(int[[i]][1], int[[i]][2], "negative")))
				idxopt <- c(idxopt, i)
			}
			else int[[i]][3] <- signs[i]
		}
	}
	int <- c(int[!sapply(1:length(int), function(i) i %in% idxopt)], left)
	return(list(def=nameInt(int, "definite"), opt=nameInt(opt, "optional")))
}

##########################
## Calling function     ##
##########################
## Current parameters are set for conversion of the Collombet model ##

#' Converts a SBML model (as a XML file) into a RE:IN-like model
#' or a igraph object
#'
#' @param filename string character of the XML file (with the extension) present in the current folder
#' @param updates string character for updates: c("async", "sync")
#' @param length string character of the maximum length of the experiments in associated experiment file
#' @param uniqueness string character of the uniqueness constraint for solution models: 
#' 			c("interactions", "full", "paths")
#' @param limit string character of the maximum number of solution models to find
#'              (maximum number of interactions to select in RE:IN!)
#' @param regulation string character of the regulatory function templates to use:
#'		c("threshold", "default", "legacy")
#' @param conditions_default string character of the default regulatory function templates
#' 	associated with every gene node created: e.g. ("0..17", "17,18", etc.)
#' @param perturbations_default string character of the default perturbations
#' 	associated with every gene node created: c("''", "'-'", "'+'", "'-+'")
#' @param igraph boolean parameter: if true, the function returns the igraph object, else the generic model
#' @return res if igraph=true then res is an igraph object associated with the model
#' 	if igraph=false then res is a generic model such as defined in the expansion code
#' 	also when igraph=false this function writes a RE:IN-like formatted model
#'
#' @export
sbml2rein <- function(filename="MODEL1610240000.xml", updates="async", length=20, 
			uniqueness="interactions", limit=1, regulation="default", 
			conditions_default="0..17", perturbations_default="", igraph=FALSE) {
	data <- xmlParse(filename)
	xml_data <- xmlToList(data)
	model <- xml_data$model
	edges <- as.vector(unlist(lapply(xml_data$model$listOfTransitions, function(e) {
			regulators <- get_regulators(e)
			regulated <- get_regulated(e)
			lapply(regulators, function(r) c(r, regulated))
		})))
	signs <- as.vector(unlist(lapply(xml_data$model$listOfTransitions, get_sign)))
	if (igraph) {
		nms <- sapply(xml_data$model$listOfTransitions, 
				function(e) as.vector(e$listOfOutputs$output['qualitativeSpecies']))
		grfs <- nameIt(lapply(xml_data$model$listOfTransitions, get_grf), nms)
		return(get_graph(edges, signs, grfs))
	}
	nodes <- as.vector(sapply(xml_data$model$listOfQualitativeSpecies, get_node))
	interactions <- get_interactions(edges, signs)
	model <- list(directives=list(directive_updates=updates,
					directive_length=length,
					directive_uniqueness=uniqueness,
					directive_limit=limit,
					directive_regulation=regulation
				), 
	 		nodes=nameIt(nodes, nodes), 
			conditions=nameIt(rep(conditions_default, length(nodes)), nodes),
			perturbations=nameIt(rep(perturbations_default, length(nodes)), nodes), 
			definites=interactions$def, 
			optionals=interactions$opt
		)
	writeModel(model, title="model.net", format="RE:IN")
	return(model)
}
