###########################################
## Compute Overlap Matrix                ##
###########################################

library(gplots)

## From outer function ##
## only for vectors    ##
## of same length      ##
outerR <- function(x, y, f) {
    dx <- length(x)
    x <- rep(x, times=dx)
    y <- rep(y, rep.int(dx, dx))
    robj <- t(matrix(sapply(1:(dx*dx), function(i) f(x[i], y[i])), ncol=dx))
    return(robj)
}

## From outer function ##
## only for matrices   ##
## of same dimension   ##
## Operation column by ##
## column              ##
outerRM <- function(x, y, f) {
    dx <- nrow(x)
    x <- matrix(rep(x, times=dx), ncol=dx*dx)
    y <- matrix(apply(y, 2, function(r) rep(r, times=dx)), ncol=dx*dx)
    robj <- matrix(sapply(1:(dx*dx), function(i) f(x[,i], y[,i])), ncol=dx)
    return(t(robj))
}

#' Compute overlap matrix
#'
#' @param crms matrix that describes the CRMs in the CRN
#' @return overlap matrix M where Mij = |CRM_i \inter CRM_j|/min(|CRM_i|, |CRM_j|)
#'
#' @export
compute_overlap_matrix <- function(crms) {
	degree <- apply(crms, 2, sum)
	n <- ncol(crms)
	degree_min <- outerR(degree, degree, min)
	diag(degree_min) <- degree
	## (Col i + Col k)_j = 0 if gene j does not belong to CRM i or k
	## (Col i + Col k)_j = 1 if gene j belongs to i XOR k
	## (Col i + Col k)_j = 2 if gene j belongs to both i or k
	jn <- apply(crms, 2, function(x) apply(crms, 2, function(y) sum(as.logical(x + y == 2))))
	dmatrix <- jn/degree_min
	## 0-1 scale
	#dmatrix <- (dmatrix-min(dmatrix))/max(dmatrix)
	dmatrix <- dmatrix/max(dmatrix)	
	## Solve case i = k
	diag(dmatrix) <- 1
	colnames(dmatrix) <- colnames(crms)
	rownames(dmatrix) <- colnames(crms)
	return(as.matrix(dmatrix))
}

###########################################
## Regulatory Module Analysis and        ##
## Visualization                         ##
###########################################

colors = c("pink", "red", "blue", "yellow", "orange", "blue4", "gray", "purple", "cyan", "magenta", "violet", "indigo", "lightblue", "antiquewhite", "aliceblue", "aquamarine", "bisque3", "blueviolet", "burlywood", "chartreuse", "coral3", "darkorchid", "darkseagreen", "forestgreen", "gold", "deepskyblue", "lavender")

## Build a CRN graph using DOT language with CRM file ##
build_graph <- function(name) {
	source(paste0(name, "_crms_cisview.R"))
	nodes <- crms$nodes
	crmscrms <- crms$crms
	crmstfs <- crms$tfs
	colorscrms <- colors[1:length(nodes)]
	names(colorscrms) <- nodes
	## Node definition        ##
	str <- paste0("digraph g1 {\n", 
		paste0(sapply(nodes, function(n) {
			paste0(c(paste0("node [style=filled color=", colorscrms[n], "]\n", n), 
				sapply(crmscrms[[n]], function(crm) paste0(n, "_", crm))
				), collapse=";\n")
			}), collapse=";\n"), ";\n")
	## CRM-Gene definition    ##
	str <- paste0(str,
		paste0(sapply(nodes, function(n) {
			paste0(sapply(crmscrms[[n]], function(crm) paste0(n, "_", crm, " -> ", n, "[fillcolor=green];\n")), collapse="")
			}),
		collapse="\n"))
	## TF binding definition  ##
	str <- paste0(str,
		paste0(sapply(nodes, function(n) {
			paste0(sapply(1:length(crmstfs[[n]]), function(i) {
					paste0(sapply(crmstfs[[n]][[i]], function(tfs) {
						paste0(sapply(tfs, function(tf) paste0(tf, " -> ", n, "_", names(crmstfs[[n]])[i], "[fillcolor=", colorscrms[n], "];\n")), collapse="")
					}), collapse="")
				}), collapse="")
			}), collapse="")		
		)
	str <- paste0(str, "\n}")
	strr <- strsplit(str, "\n")[[1]]
	str <- paste0(strr[strr != ""], collapse="\n")
	write(str, file=paste0("graph_cisview_", name, ".dot"))
}

## Build matrix for computation of overlap/heatmaps ##
## M_{i,j} = 1 if TF_j in CRM_i, otherwise 0        ##
build_matrix <- function(name, save=F) {
	source(paste0(name, "_crms_cisview.R"))
	nodes <- crms$nodes
	crmscrms <- crms$crms
	crmstfs <- crms$tfs
	rowNames <- unlist(sapply(nodes, function(n) {
				ls <- unlist(crmscrms[[n]])
				names(ls) <- rep(n, length(ls))
				ls
			}))
	names(rowNames) <- sapply(names(rowNames), function(s) strsplit(s, "[.]")[[1]][1])
	crms_modules <- matrix(0, ncol=length(nodes), nrow=length(rowNames))
	colnames(crms_modules) <- nodes
	rownames(crms_modules) <- sapply(1:length(rowNames), function(i) paste0(names(rowNames)[i], "_", rowNames[i]))
	for (i in 1:length(rowNames)) for (j in 1:length(nodes)) {
		crms_modules[i, j] <- ifelse(nodes[j] %in% crmstfs[[names(rowNames)[i]]][[rowNames[i]]], 1, 0)
	}
	if (save) save(crms_modules, file=paste0(name, "crms_cisview.rData"))
	return(crms_modules)
}

## Heatmap functions                                ##
## For visualization of CRMs                        ##
heatmap_c <- function(mm, title="") {
	return(heatmap.2(mm, col=RColorBrewer::brewer.pal(9,"Blues"),
		Colv=NA, Rowv=NA, dendrogram="none",
		revC=F, symm=T, trace="none",
		main=title, cexRow=1, cexCol=1, margins=c(12, 12), srtCol=45))
}

## For visualization of matrix classified by gene   ##
heatmap_c2 <- function(mm, title="") {
	return(heatmap.2(mm, col=RColorBrewer::brewer.pal(9,"Blues"),
		revC=T, symm=T, trace="none",
		main=title, cexRow=1, cexCol=1, margins=c(12, 12), srtCol=45))
}

margins=c(7, 4, 4, 7)+0.1
size=c(800, 750)

heatmapIt <- function(name, byGene=FALSE) {
	crms_modules <- build_matrix(name)
	tt <- t(crms_modules)[, apply(t(crms_modules), 2, function(x) var(x) != 0)]
	cortt <- 1-cor(tt, method="spearman")
	par(mar=margins) 
	png(filename=ifelse(byGene, 'heatmap_by_gene.png', 'heatmap.png'), width=size[1], height=size[2])
	if (byGene) heatmap_c2(cortt) else heatmap_c(cortt)
	graphics.off()
}

overlapMatrix <- function(name) {
	crms_modules <- build_matrix(name)
	tt <- t(crms_modules)[, apply(t(crms_modules), 2, function(x) var(x) != 0)]
	M <- compute_overlap_matrix(tt)
	par(mar=margins) 
	png(filename='overlap.png', width=size[1], height=size[2])
	heatmap_c(1-M)
	graphics.off()
}


