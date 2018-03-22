source("io_expansion.R")

#########################
## PARTIAL EXPANSION   ##
#########################

#' Trim the model to get a partial expansion
#'
#' @param crms Regulatory Module Set object
#' @return crms "trimmed" Regulatory Module Set object
#'
#' @export
partialExp <- function(crms) {
	## Only selects nonempty CRMs from genes having strictly more than  1   ##
	## nonempty CRM                                                         ##
	## An empty CRM is a CRM with NO known TF binding to it                 ##
	## Gets the indices of those "invalid" CRMs                             ##
	idx <- sapply(crms$tfs, function(lsls) length(removeNULLInList(lsls)) < 2)
	## Removes the associated columns from the model                        ##
	crms$crms[idx] <- NULL
	crms$tfs[idx] <- NULL
	## Remove empty CRMs among the CRMs of genes                                  ##
	idx <- lapply(crms$tfs, function(lsls) which(!sapply(lsls, is.null)))
	crms$tfs <- sapply(names(crms$crms), function(n) crms$tfs[[n]][idx[[n]]])
	crms$crms <- sapply(names(crms$crms), function(n) crms$crms[[n]][idx[[n]]])
	return(crms)
}

#################
## EXPANSION   ##
#################

#' Expand model
#'
#' @param model generic model obtained with function readModel
#' @param method = c("positive", "negative", c("positive", "negative")) 
#'         type of interaction signs expected for CRM->gene links
#' @param type = c("partial", "full") type of expansion
#' @param inferTF boolean: should the solver try to infer direct TF binding to CRMs?
#' @param crms Regulatory Module Set object
#' @return model expanded model as generic model (see specifications file)
#'
#' @export
expansion_aux <- function(model, method, type, inferTF, crms) {
	if (type=="partial") crms <- partialExp(crms)
	nodesPertCond <- getNewNodes(model, crms)
	## Add CRM-gene links in model       ##
	crmLinks <- getCRMLinks(model, crms, method)
	## Add known TF-CRM (optional) links ##
	tfLinks <- getTFCRMLinks(model$optionals, crms)
	## Add known TF-CRM (definite) links ##
	tfLinksDef <- getTFCRMLinks(model$definites, crms)
	## Add unknown TF-CRM links in model ##
	## /!\ Only if the TF bindings data  ##
	## /!\ is really trustworthy         ##
	#allTFCRMLinks <- getAllTFCRMLinks(crms)
	#unknownTFLinks <- removeInt(c(tfLinksDef, tfLinks), 
	#		allTFCRMLinks, f=posInListList3)
	## Comment the following line in     ##
	## case of the precedings one be     ##
	## decommented                       ##
	unknownTFLinks <- NULL
	## Gene-gene interactions            ##
	## (no assessed TF binding)          ##
	geneLinks <- removeInt(tfLinks, model$opt)
	geneLinksDef <- removeInt(tfLinksDef, model$def)
	## Add an optional interaction       ##
	## between g1 and every CRM of g2    ##
	## for every interaction g1 -> g2    ##
	## where g1 is a TF                  ##
	if (inferTF) {
		## Get all TFs               ##
		idx <- names(sapply(crms$isTF, function(e) e))
		tfsNodes <- crms$nodes[unlist(sapply(idx, function(id) which(crms$nodes == id)))]
		## Screen every gene-gene    ##
		## interaction g1 -> g2      ##
		inferTFLinks <- unwrap(sapply(c(geneLinks, geneLinksDef), function(e) {
				## If g1 is a TF     ##
				ls <- names(crms$tfs[e[2]][[1]])
				if (e[1] %in% tfsNodes && length(ls) > 0) {
					lapply(names(crms$tfs[e[2]][[1]]), 
						function(crm) 
						c(e[1], paste0(e[2], "_", crm), e[3]))
				}
			}))
		## Remove the associated gene-gene   ##
		## interactions from the sets        ##
		geneLinks <- removeInt(inferTFLinks, geneLinks, f=posInListList3)
		geneLinksDef <- removeInt(inferTFLinks, 
					geneLinksDef, 
					f=posInListList3)		
		unknownTFLinks <- c(unknownTFLinks, inferTFLinks)		
	}
	else {
		inferTFLinks <- NULL
		geneLinks <- trimUnknownTFLinks(geneLinks, crms)
		geneLinksDef <- trimUnknownTFLinks(geneLinksDef, crms)
	}
	model$nodes <- nodesPertCond[[1]]
	model$perturbations <- nodesPertCond[[2]]
	model$conditions <- nodesPertCond[[3]]
	model$optionals <- c(crmLinks, tfLinks, unknownTFLinks, geneLinks)
	## Removes NULL elements ##
	model$optionals <- removeNULLInList(model$optionals)
	model$definites <- c(tfLinksDef, geneLinksDef)
	## Removes NULL elements ##
	model$definites <- removeNULLInList(model$definites)
	## Prints stats          ##
	statsExpandedModel(unknownTFLinks, tfLinks, tfLinksDef, 
			geneLinks, geneLinksDef, crmLinks, 
			model, inferTFLinks)
	return(model)
}

#################
## ONE-LINER   ##
#################

#' Get and write expanded model
#'
#' @param crms Regulatory Module Set object
#' @param filename RE:IN .net file that contains the RE:IN formatted model
#' @param title the name of resulting model file
#' @param method = c("positive", "negative", c("positive", "negative")) 
#'         type of interaction signs expected for CRM->gene links
#' @param type = c("partial", "full") type of expansion
#' @param inferTF should the solver try to infer direct TF binding to CRMs?
#' @param format "RE:IN"-like or "expanded"
#' @return model expanded model as a generic model object
#'
#' @export
expansion <- function(crms, filename, title, 
		method=c("positive", "negative"), 
		type="partial", 
		inferTF=c(T, F), format="RE:IN") {
	model <- readModel(filename)
	modelCRM <- expansion_aux(model=model, 
				method=method, 
				type=type, 
				inferTF=inferTF, 
				crms=crms)
	writeModel(modelCRM, title=title, format=format)
	return(modelCRM)
}

#' Get and write reduced model from expanded model
#' i.e. merge RM nodes with their corresponding target gene nodes
#' WARNING: Interactions TF->gene present in the original model and removed in
#' the expanded model because the TF did not bind to any of the RMs of the gene
#' cannot be retrieved
#'
#' @param filename RE:IN .net file that contains the RE:IN formatted expanded model
#' @param title the name of resulting model file
#' @param format "RE:IN"-like or "expanded"
#' @return model expanded model as a generic model object
#'
#' @export
reduction <- function(filename, crms, title=NULL, format="RE:IN") {
	title <- if (is.null(title)) paste0(strsplit(filename, ".net")[[1]][1], "_reduced.net") else title
	model <- readModel(filename)
	crmsList <- lapply(1:length(crms$crms), 
		function(i) lapply(crms$crms[[i]], 
					function(crm) paste0(names(crms$crms)[i], "_", crm)
			))
	names(crmsList) <- names(crms$crms)
	directives <- model$directives
	nodes <- model$nodes[sapply(model$nodes, function(e) !(e %in% unlist(crmsList)))]
	conditions <- model$conditions[nodes]
	perturbations <- model$perturbations[nodes]
	optionals <- reduceInteractions(model$optionals, crmsList)
	definites <- reduceInteractions(model$definites, crmsList)
	modelReduced <- createModel(directives, nodes, conditions, perturbations, optionals, definites)
	writeModel(modelReduced, title=title, format=format)
	return(modelReduced)
}
