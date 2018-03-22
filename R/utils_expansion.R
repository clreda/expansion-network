###########
## TOOLS ##
###########

getIdxOfColumn <- function(coln, ls) return(which(grepl(coln, names(ls))))
getIdxOfData <- function(col, ls) return(which(col == names(ls)))

#' Return the indices of NULL elements in input list list
#'
#' @param lsls input -named- list list
#' @return idx vector of integer indices of NULL elements in input list list
#'
#' @export
idxNULLInList <- function(lsls) {
	if (is.null(lsls)) return(NULL)
	return(which(sapply(lsls, function(ll) is.null(ll[[1]]))))
}

#' Remove NULL elements in input list list
#'
#' @param lsls input -named- list list
#' @return lsls input list list without any NULL element
#'
#' @export
removeNULLInList <- function(lsls) {
	idx <- idxNULLInList(lsls)
	if (length(idx) > 0) return(lsls[setdiff(1:length(lsls), idx)])
	return(lsls)
}

#' Concatenate strings in list
#'
#' @param v list or vector of strings
#' @return str concatenation of all strings in v
#'
#' @export
paste1 <- function(v) return(paste0(v, collapse=""))
paste2 <- function(ls) return(sapply(ls, paste1))
paste3 <- function(ls) return(sapply(ls, 
			function(e) paste1(c(e[1], strsplit(e[2], "_")[[1]][1], e[3]))))

#' Return indices of occurrences of element in list list
#'
#' @param e element (that can be casted into a list of character strings)
#' @param lsls list of list of e type elements
#' @return i index of first occurrence of e in lsls
#'
#' @export
posInListList <- function(e, lsls) return(which(grepl(paste1(e), paste2(lsls))))
posInListList2 <- function(e, lsls) return(which(grepl(paste1(e), paste3(lsls))))
posInListList3 <- function(e, lsls) return(which(grepl(paste1(e[1:2]), paste2(lsls))))

## Remove interactions from set b1 present in set b2       ##
removeInt <- function(b1, b2, f=posInListList2) {
	if (length(b2) == 0 || length(b1) == 0) return(b2)
	idx <- unique(unlist(sapply(1:length(b2), 
			function(i) if (length(f(b2[[i]], b1))>0) i)))
	return(b2[setdiff(1:length(b2), idx)])
}

nameIt <- function(cc, nm) {
	names(cc) <- nm
	return(cc)
}
nameItself <- function(cc) return(nameIt(cc, cc))

unwrap <- function(e, times=1) {
	while (times > 0) {
		e <- unlist(e, recursive=F)
		times <- times-1
	}
	return(e)
}

#################################
## EXPANSION                   ##
#################################

##___________________________________________________##
## Operations on nodes, conditions and perturbations ##
##___________________________________________________##

getNewNodes_aux <- function(crmNodes, crms, modelcol) {
	return(nameIt(sapply(crmNodes, function(x) { 
			y <- strsplit(x, "_")[[1]][1]
			if (!is.null(crms$crms[y][[1]])) { 
				modelcol[getIdxOfData(y, modelcol)]
			}
		}), crmNodes))
}

rep_nodes <- function(crmNodes, value) {
	return(nameIt(rep(value, length(crmNodes)), crmNodes))
}

getNewNodes <- function(model, crms) {
	## Add CRMs as nodes               ##
	crmNodes <- nameItself(unlist(sapply(model$nodes, function(x) 
			if (!is.null(crms$crms[x][[1]])) {
				paste0(x, "_", crms$crms[x][[1]])
			}
		)))
	## with no perturbations           ##
	crmPNodes <- rep_nodes(crmNodes, "")
	## with all regulation conditions  ##
	## except for threshold ones       ##
	crmCNodes <- rep_nodes(crmNodes, "0..17")
	return(list(c(model$nodes, crmNodes), 
		c(model$perturbations, crmPNodes), 
		c(model$conditions, crmCNodes)
		))
}

##___________________________________________________##
## Operations on interactions CRM -> gene            ##
##___________________________________________________##

getCRMLinks <- function(model, crms, method) {
	addCRMGene <- function(crm, g) {
		return(lapply(method, function(s) c(paste0(g, "_", crm), g, s)))
	}
	crmLinks <- unwrap(sapply(crms$nodes, 
		function(n) sapply(unlist(crms$crms[[n]]), 
				function(crm) addCRMGene(crm, n)
				)
		))
	return(crmLinks)
}

##___________________________________________________##
## Operations on interactions TF -> CRM              ##
##___________________________________________________##

getTFCRMLinks <- function(modelInt, crms) {
	addTFCRM <- function(tf, crm, g, s) {
		return(c(tf, paste0(g, "_", crm), s))
	}
	tfLinks <- sapply(modelInt, function(e) {
			## Looks into TFs of regulated gene     ##
			## in considered optional interaction e ##
			lapply(1:length(crms$tfs[e[2]][[1]]), function(i) {
				## Looks for a CRM of regulated gene ##
				## having regulator as TF            ##
				ls <- unlist(crms$tfs[e[2]][[1]][[i]])
				crm <- names(crms$tfs[e[2]][[1]])[i]
				## If there are such CRMs, create    ##
				## link between regulator and them   ##
				if (length(ls)>0 && e[1] %in% ls) {
					addTFCRM(e[1], crm, e[2], e[3])
				}
			})		
		})
	if (length(tfLinks) == 0) tfLinks <- NULL
	tfLinks <- unwrap(tfLinks)[!sapply(unwrap(tfLinks), function(e) is.null(e[[1]]))]
	return(tfLinks)
}

##___________________________________________________##
## Gets all possible TF -> CRM interactions from     ##
## the CRM object                                    ##
##___________________________________________________##

getAllTFCRMLinks <- function(crms) {
	addTFCRM <- function(tf, crm, g, ss) {
		return(lapply(ss, function(s) c(tf, paste0(g, "_", crm), s)))
	}
	tfLinks <- sapply(crms$nodes, function(n) {
			lapply(1:length(crms$tfs[n][[1]]), function(i) {
				ls <- unlist(crms$tfs[n][[1]][[i]])
				crm <- names(crms$tfs[n][[1]])[i]
				if (length(ls) > 0) {
					lss <- lapply(ls, 
					function(tf) {
					addTFCRM(tf, crm, n, c("positive", "negative"))
					})
				}
				else lss <- NULL
			})		
		})
	if (length(tfLinks) == 0) tfLinks <- NULL
	tfLinks <- unwrap(tfLinks, 3)[!sapply(unwrap(tfLinks, 3), 
				function(e) is.null(e[[1]]))]
	return(tfLinks)	
}

##___________________________________________________##
## Removes interactions g1 -> g2 where g1 is a TF and##
## g2 has CRMs (assuming known TF bindings have      ##
## already been trimmed)                             ##
##___________________________________________________##

trimUnknownTFLinks <- function(geneLinks, crms) {
	## TF->gene, where gene has at most 1 CRM, MUST be kept ##
	return(geneLinks[!sapply(geneLinks, function(e) !is.null(crms$isTF[[e[1]]]) && crms$isTF[[e[1]]] && length(crms$crms[[e[2]]]) > 1)])
}

#################################
## REDUCTION                   ##
#################################

reduceInteractions <- function(interactions, crmsList) {
	trimmedInt <- lapply(interactions, function(e) {
			if (e[1] %in% unlist(crmsList)) NULL
			else {
				if (all(!(e %in% unlist(crmsList)))) e
				else {
					idx <- sapply(crmsList, function(ls) e[2] %in% ls)
					c(e[1], names(crmsList)[idx], e[3])
				}
		}
	})
	trimmedInt <- trimmedInt[!sapply(trimmedInt, is.null)]
	trimmedInt <- lapply(unique(sapply(trimmedInt, function(e) paste0(e, collapse="#"))), function(e) strsplit(e, "#")[[1]])
	return(trimmedInt)
}

#################################
## PRINTING                    ##
#################################

printPretty <- function(args) return(print(paste0(args, collapse=" ")))

statsModel <- function(model) {
	printPretty(c("STATS-----------------------"))
	printPretty(c("Nodes:", length(model$nodes)))
	printPretty(c("CRM nodes:", length(model$nodes[grepl("_CM", names(model$nodes))])))
	printPretty(c("Optional edges:", length(model$opt)))
	printPretty(c("Definite edges:", length(model$def)))
}

statsExpandedModel <- function(unknownTFLinks, tfLinks, tfLinksDef, 
				geneLinks, geneLinksDef, crmLinks, 
				model, inferTFLinks) {
	printPretty(c("EXPANSION-------------------"))
	printPretty(c("Unsupported (opt.) TF-gene edges:", length(unknownTFLinks)))
	printPretty(c("Supported opt. TF-gene edges:", length(tfLinks)))
	printPretty(c("Supported def. TF-gene edges:", length(tfLinksDef)))
	printPretty(c("Gene-gene opt. edges:", length(geneLinks)))
	printPretty(c("Gene-gene def. edges:", length(geneLinksDef)))
	printPretty(c("CRM-gene edges:", length(crmLinks)))
	if (!is.null(inferTFLinks)) {
		printPretty(c("TF-CRM edges to infer:", length(inferTFLinks)))
	}
	printPretty(c("Nodes:", length(model$nodes)))
	printPretty(c("CRM nodes:", length(model$nodes[grepl("_CM", names(model$nodes))])))
	printPretty(c("Optional edges:", length(model$opt)))
	printPretty(c("Definite edges:", length(model$def)))
	printPretty(c("----------------------------"))
}
