#!/usr/bin/env Rscript

# Load packages only from local Rslamdunk library 
libLoc = .libPaths()[grep("Rslamdunk",.libPaths())]

# Check if libraries are available, install otherwise
source(paste(libLoc,'/../checkLibraries.R',sep=""))

checkLib(libLoc)

library(getopt, lib.loc = libLoc)

spec = matrix(c(
				'help'      , 'h', 0, "logical","print the usage of the command",
				'utr' 		, 'u', 0, "logical","utr plotting",
				'inputFile', "i", 2,"character","tsv table of mutations per position",
				'outputFile', "o", 2,"character","output pdf file name"
		),ncol = 5,byrow=T)

opt = getopt(spec)

if ( !is.null(opt$help) || length(opt)==1 ) {
	#get the script name
	cmd = commandArgs(FALSE)
	self = strsplit(cmd[grep("--file",cmd)],"=")[[1]][2]
	cat(basename(self),": Create mismatches per read/UTR position plots.\n\n")
	#print a friendly message and exit with a non-zero error code
	cat(getopt(spec,command = self,usage=T))
	q(status=1);
}

positionLabel = "Position on read"
mutationLabel = "% of reads with mutation"

if( !is.null(opt$utr)) {
	positionLabel = "Position at 3' UTR end (200 bp upstream)"
	mutationLabel = "% of UTRs with mutation"
}

if ( is.null(opt$inputFile) ) stop("arg input must be specified")
if ( is.null(opt$outputFile) ) { opt$outputFile = paste(opt$inputFile, ".pdf", sep="") }


mut = read.table(opt$inputFile, comment.char = "#")

if (is.null(mut$V6)) {
	mut$V6 = mut$V5
}

#mut = read.table("test_mut_bowtie.csv")

#totalFwd = mut[1,1]
#totalRev = mut[1,2]
#tcFwd = mut[1,3]
#tcRev = mut[1,4]

#mut = mut[-1,]

counts = rbind(c(mut$V1)/c(mut$V5) * 100, c(mut$V2)/c(mut$V6) * 100)
countsTC = rbind(c(mut$V3)/c(mut$V5) * 100, c(mut$V4)/c(mut$V6) * 100)

##################################################################
# Workaround for 0 counts (need to work out what's going on there

counts[is.nan(counts)] = 0
countsTC[is.nan(countsTC)] = 0

##################################################################
pdf(opt$outputFile, width=10, height=10)
par(mfrow=c(2,1))

# Scale to next 10
barplot(counts, beside=T, names.arg=1:nrow(mut), main="All mutations", ylim=c(0,max(10,ceiling(counts / 10) * 10)), xlab=positionLabel, ylab=mutationLabel, legend=c("forward", "reverse"))
#barplot(counts, beside=T, names.arg=1:nrow(mut), main="All mutations", ylim=c(0,10), xlab=positionLabel, ylab=mutationLabel, legend=c("forward", "reverse"))
# Scale to next 1
barplot(countsTC, beside=T, names.arg=1:nrow(mut), main="T->C on fwd, A->G on rev", ylim=c(0,max(1,ceiling(countsTC))), xlab=positionLabel, ylab=mutationLabel, legend=c("forward", "reverse"))
#barplot(countsTC, beside=T, names.arg=1:nrow(mut), main="T->C on fwd, A->G on rev", ylim=c(0,1), xlab=positionLabel, ylab=mutationLabel, legend=c("forward", "reverse"))

dev.off()