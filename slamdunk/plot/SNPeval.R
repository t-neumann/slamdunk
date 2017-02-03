#!/usr/bin/env Rscript

# Script to look at SNP distributions along UTRs ranked by # T>C SNPs
# 
# Author: Tobias Neumann, Zuber group, Institute for Molecular Pathology
# Email: tobias.neumann@imp.ac.at
###############################################################################

# Load packages only from local Rslamdunk library 
libLoc = .libPaths()[grep("Rslamdunk",.libPaths())]

# Check if libraries are available, install otherwise
source(paste(libLoc,'/../checkLibraries.R',sep=""))

checkLib(libLoc)

library(getopt, lib.loc = libLoc)

spec = matrix(c(
				'help'      , 'h', 0, "logical","print the usage of the command",
				'inputFile', "i", 2,"character","tsv table of snp vs tc count files",
				'coverageCutoff', "c", 2,"numeric","coverage cutoff for calling variants",
				'variantFraction', "v", 2,"numeric","variant fraction cutoff for calling variants",
				'outputFile', "o", 2,"character","output pdf file name"
		),ncol = 5,byrow=T)

opt = getopt(spec)

if ( !is.null(opt$help) || length(opt)==1 ) {
	#get the script name
	cmd = commandArgs(FALSE)
	self = strsplit(cmd[grep("--file",cmd)],"=")[[1]][2]
	cat(basename(self),": Create mismatch plots from rate tabs.\n\n")
	#print a friendly message and exit with a non-zero error code
	cat(getopt(spec,command = self,usage=T))
	q(status=1);
}


if ( is.null(opt$inputFile) ) stop("arg rateTab must be specified")
if ( is.null(opt$outputFile) ) { opt$outputFile = "out.pdf" }
if ( is.null(opt$coverageCutoff) ) { opt$coverageCutoff = 0 }
if ( is.null(opt$variantFraction) ) { opt$variantFraction = 0 }

tricubeMovingAverage <-function (x, span = 0.5, full.length = TRUE) {
	
	n <- length(x)
	width <- span * n
	hwidth <- as.integer(width%/%2L)
	if (hwidth <= 0L) 
		return(x)
	width <- 2L * hwidth + 1L
	u <- seq(from = -1, to = 1, length = width) * width/(width + 
				1)
	tricube.weights <- (1 - abs(u)^3)^3
	tricube.weights <- tricube.weights/sum(tricube.weights)
	if (!full.length) 
		return(as.vector(filter(x, tricube.weights), mode = "numeric")[(hwidth + 1):(n - hwidth)])
	z <- numeric(hwidth)
	x <- as.vector(filter(c(z, x, z), tricube.weights), mode = "numeric")[(hwidth + 1):(n + hwidth)]
	cw <- cumsum(tricube.weights)
	x[1:hwidth] <- x[1:hwidth]/cw[(width - hwidth):(width - 1)]
	x[(n - hwidth + 1):n] <- x[(n - hwidth + 1):n]/cw[(width - 1):(width - hwidth)]
	x
}

rescale <- function(x, new, old = range(x)) {
	new[1] + (x - old[1])/(old[2] - old[1]) * 
			(new[2] - new[1])
}

GSEAplot <- function(counts, snps, ...) {
	
	num <- length(counts)
	
	sel = rep.int(FALSE, num)
	
	sel[snps] <- TRUE
	
	countOrder <- order(counts, na.last = TRUE, decreasing = TRUE)
	counts <- counts[countOrder]
	sel <- sel[countOrder]
	
	snpLoc <- which(sel)
	
	col.bars <- "black"
	
	ylim <- c(-1, 1.5)
	
	plot(1:num, xlim = c(0, num), ylim = c(0, 2.1), type = "n", 
			axes = FALSE, ylab = "", ...)
	
	lwd <- 50/length(snpLoc)
	lwd <- min(1.9, lwd)
	lwd <- max(0.2, lwd)
	
	barlim <- ylim[2] - c(1.5, 0.5)
	rect.yb <- 0
	rect.yt <- 0.5
	rect(0.5, 0, num + 0.5, 0.5, col = "pink", border = NA)
	
	segments(snpLoc, barlim[1], snpLoc, barlim[2]/2, lwd = lwd, col = "black")
	segments(snpLoc, barlim[2]/2, snpLoc, barlim[2]/2 * 2, lwd = lwd, col = "black")
	
	axis(side = 2, at = 0.5, padj = 3.8, cex.axis = 0.85, 
			labels = "High # T>C reads", tick = FALSE)
	axis(side = 4, at = 0.5, padj = -3.8, cex.axis = 0.85, 
			labels = "Low # T>C reads", tick = FALSE)
	prob <- (10:0)/10
	axis(at = seq(1, num, len = 11), side = 1, cex.axis = 0.7, 
			las = 2, labels = format(quantile(counts, p = prob), 
					digits = 1))
	
	ave.enrich1 <- length(snpLoc)/num
	worm1 <- tricubeMovingAverage(sel, span = 0.45)/ave.enrich1
	
	r.worm1 <- c(0, max(worm1))
	worm1.scale <- rescale(worm1, new = c(1.1 , 2.1 ), old = r.worm1)
	
	lines(x = 1:num, y = worm1.scale, col = "black", lwd = 2)
	abline(h = rescale(1, new = c(1.1 , 2.1), old = r.worm1), lty = 2)
	axis(side = 2, at = c(1.1 , 2.1 ), cex.axis = 0.8, labels = c(0, format(max(worm1), digits = 2)))
	axis(side = 2, labels = "Enrichment", at = 1.6 , padj = -0.6, tick = FALSE, cex.axis = 0.8)
}

pdf(opt$outputFile)

minCounts = round(opt$coverageCutoff * opt$variantFraction)

table = read.delim(opt$inputFile,header=FALSE,col.names = c("name","count","unmasked","masked","snp"))

table = table[table$unmasked >= minCounts,]

table = table[table$count >= quantile(table$count, 0.75),]

# testTab = table
# testTab$unmasked = rank(testTab$unmasked)
# testTab = testTab[order(testTab$unmasked, decreasing=TRUE),]
# testTab$masked = rank(testTab$masked)

par(mfrow=c(2,1))

blindTest = wilcox.test(table$unmasked ~ table$snp == "1", alternative = "less")
maskedTest = wilcox.test(table$masked ~ table$snp == "1", alternative = "less")

blindPvalue = blindTest$p.value

if (blindPvalue < 0.01) {
	blindPvalue = "< 0.01"
} else {
	blindPvalue =  paste("= ",round(blindPvalue,digits=2),sep="")
}

maskedPvalue = maskedTest$p.value

if (maskedPvalue < 0.01) {
	maskedPvalue = "< 0.01"
} else {
	maskedPvalue =  paste("= ",round(maskedPvalue,digits=2),sep="")
}

GSEAplot(table$unmasked, which(table$snp == 1), main="Blind", xlab = paste("Mann-Whitney-U:  p-value ",blindPvalue,sep=""))
GSEAplot(table$masked, which(table$snp == 1), main="SNP-masked", xlab = paste("Mann-Whitney-U:  p-value ",maskedPvalue,sep=""))

# wilcox.test(testTab$masked ~ testTab$snp == "1")
# wilcox.test(testTab$unmasked ~ testTab$snp == "1")
#
# ks.test(testTab$masked, testTab$unmasked)
# ks.test(testTab$unmasked[testTab$snp == "0"],testTab$unmasked[testTab$snp == "1"])

dev.off()