#!/usr/bin/env Rscript

# Script to plot TC context rates of reads
###############################################################################
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
  'rateTab', "f", 2,"character","tsv table of rate files",
  'outputFile', "O", 2,"character","output pdf file name"
),ncol = 5,byrow=T)

opt = getopt(spec)

if ( !is.null(opt$help) || length(opt)==1 ) {
  #get the script name
  cmd = commandArgs(FALSE)
  self = strsplit(cmd[grep("--file",cmd)],"=")[[1]][2]
  cat(basename(self),": Create mismatch plots from rate tabs.\n\n")
  #print a friendly message and exit with a non-zero error code
  cat(getopt(spec,command = self,usage=T))
  q(status=1)
}


if ( is.null(opt$rateTab) ) stop("arg rateTab must be specified")
if ( is.null(opt$outputFile) ) { opt$outputFile = "out.pdf" }

library(ggplot2, lib.loc = libLoc)
library(gridExtra, lib.loc = libLoc)
library(dplyr, lib.loc = libLoc)

rates = read.table(opt$rateTab,stringsAsFactors=FALSE,col.names = c("sample","file"),comment.char = "")

pdf(opt$outputFile)
plotList = list()

for (i in 1:nrow(rates)) {
  curTab = read.table(rates$file[i],stringsAsFactors=FALSE,header=TRUE)
  
  subFront = curTab[1:2,]
  subBack = curTab[4:5,]
  names(subBack) = curTab[3,]
  
  #subFront = read.table(rates$file[i],stringsAsFactors=FALSE,header=TRUE, nrow=1)
  #subBack = read.table(rates$file[i],stringsAsFactors=FALSE,header=TRUE, nrow=1,skip=2)
  
  printTabFront = data.frame(contexts=rep(names(subFront),each=2),strand = factor(rep(c("+","-"),ncol(subFront)),levels=c("+","-")),
                             rate_percent = as.numeric(unlist(subFront)))
  printTabBack = data.frame(contexts=rep(names(subBack),each=2),strand = factor(rep(c("+","-"),ncol(subBack)),levels=c("+","-")),
                            rate_percent = as.numeric(unlist(subBack)))
  
  printTabFront$rate_percent = printTabFront$rate_percent / sum(printTabFront$rate_percent)
  printTabBack$rate_percent = printTabBack$rate_percent / sum(printTabBack$rate_percent)
  
  # Ignore N contexts for now
  printTabFront = printTabFront[-grep("NT",printTabFront$contexts),]
  printTabBack = printTabBack[-grep("TN",printTabBack$contexts),]
  
  curPlot = qplot(x=contexts, y=rate_percent, fill=strand,data=printTabFront) + geom_bar(stat="identity") + geom_text(aes(label = round(rate_percent,digits=2)), size = 3, hjust = 0.5, vjust = 1.5, position = "stack") + ylab("TC context percent %") + xlab(rates$sample[i]) +
    theme(text = element_text(size=6),axis.text.x = element_text(size=6), plot.title = element_text(size=10))
  plotList[[length(plotList)+1]] <- curPlot + ylim(0.0,1.0) + ggtitle("5' T->C context")
  curPlot = qplot(x=contexts, y=rate_percent, fill=strand,data=printTabBack) + geom_bar(stat="identity") + geom_text(aes(label = round(rate_percent,digits=2)), size = 3, hjust = 0.5, vjust = 1.5, position = "stack") + ylab("TC context percent %") + xlab(rates$sample[i]) +
    theme(text = element_text(size=6),axis.text.x = element_text(size=6),plot.title = element_text(size=10))
  plotList[[length(plotList)+1]] <- curPlot + ylim(0.0,1.0) + ggtitle("3' T->C context")
}

do.call(grid.arrange,  plotList)

dev.off()

#signal success and exit.
q(status=0)		