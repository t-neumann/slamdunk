#!/usr/bin/env Rscript

# Script to plot pairwise correlations and PCA
# 
# Author: Tobias Neumann, Zuber group, Institute for Molecular Pathology
# Email: tobias.neumann@imp.ac.at
###############################################################################

# Helper

my_panel_cor <- function(x, y, digits=2, prefix="", cex.cor, ...)
{
	usr <- par("usr"); on.exit(par(usr))
	par(usr = c(0, 1, 0, 1))
	
	
	toUse = which(is.finite(x) & is.finite(y) & (x|y>0))
	r <- abs(cor(x[toUse], y[toUse]))
	
	
	txt <- format(c(r, 0.123456789), digits=digits)[1]
	txt <- paste(prefix, txt, sep="")
	if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
	text(0.5, 0.5, txt, cex = cex.cor * r)
}

my_panel_smooth <- function(x, y,lcol="red")

{
	smoothScatter(x,y,add=T)
	abline(0,1,col=lcol)
}

# Load packages only from local Rslamdunk library 
libLoc = .libPaths()[grep("Rslamdunk",.libPaths())]

# Check if libraries are available, install otherwise
source(paste(libLoc,'/../checkLibraries.R',sep=""))

checkLib(libLoc)

library(getopt, lib.loc = libLoc)

spec = matrix(c(
				'help'      , 'h', 0, "logical","print the usage of the command",
				'sampleTab', "i", 2,"character","csv table of sample counts",
				'outputPrefix', "o", 2,"character","output file name prefix"
		),ncol = 5,byrow=T)

opt = getopt(spec)

if ( !is.null(opt$help) || length(opt)==1 ) {
	#get the script name
	cmd = commandArgs(FALSE)
	self = strsplit(cmd[grep("--file",cmd)],"=")[[1]][2]
	cat(basename(self),": Compute sample comparison statistics from sample counts.\n\n")
	#print a friendly message and exit with a non-zero error code
	cat(getopt(spec,command = self,usage=T))
	q(status=1);
}


if ( is.null(opt$sampleTab) ) stop("arg sampleTab must be specified")
if ( is.null(opt$outputPrefix) ) { opt$outputPrefix = "sampleCorrelation" }

rates = read.table(opt$sampleTab,header=TRUE,sep=";", comment.char = "")

if (ncol(rates) < 6) {
	print("No need for calculating pairwise statistics for single sample")
	quit(status=0)
}

library(RColorBrewer, lib.loc = libLoc)
library(lattice, lib.loc = libLoc)
library(matrixStats, lib.loc = libLoc)

values = data.matrix(rates[,c(5:ncol(rates))])

##################################################
# PCA
##################################################

rowVariances = rowVars(data.matrix(values))

select = order(rowVariances, decreasing = TRUE)[seq_len(min(500, length(rowVariances)))]

pca = prcomp(t(values[select, ]))

if (ncol(values) == 2) {
	col = brewer.pal(3, "Paired")[1:2]	
} else if (ncol(values) > 12) {
	getPalette = colorRampPalette(brewer.pal(9, "Set1"))
	col = getPalette(ncol(values))
} else {
	col = brewer.pal(ncol(values), "Paired")
}

# Get amount of explained variance (see summary(pca))
varianceProportion = pca$sdev ^ 2 / sum(pca$sdev ^ 2)

pdf(paste(opt$outputPrefix,"_PCA.pdf",sep=""))

if (ncol(values) > 12) {

	xyplot(PC2 ~ PC1, groups = colnames(values), data = as.data.frame(pca$x),
		pch = 20, cex = 2, aspect = "iso", col = col, xlab = paste("PC1 (", round(varianceProportion[1],digits=2), " variance)",sep=""),
		ylab = paste("PC2 (", round(varianceProportion[2],digits=2), " variance)",sep=""),
	)

} else {
	
	xyplot(PC2 ~ PC1, groups = colnames(values), data = as.data.frame(pca$x),
			pch = 20, cex = 2, aspect = "iso", col = col, main = draw.key(key = list(rect = list(col = col),
							text = list(colnames(values)), rep = FALSE)), xlab = paste("PC1 (", round(varianceProportion[1],digits=2), " variance)",sep=""),
			ylab = paste("PC2 (", round(varianceProportion[2],digits=2), " variance)",sep=""),
	)
	
}

dev.off()

##################################################
# Pairwise correlations
##################################################

if (ncol(values) <= 12) {

	pdf(paste(opt$outputPrefix,"_pairwiseCorrelation.pdf",sep=""))

	pairs(values,upper.panel=my_panel_smooth,lower.panel=my_panel_cor)

	dev.off()
}

#signal success and exit.
q(status=0)		
