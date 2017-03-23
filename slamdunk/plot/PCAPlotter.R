#!/usr/bin/env Rscript

# Plot PCA based on readcounts in UTRs

# Copyright (c) 2015 Tobias Neumann, Philipp Rescheneder.
#
# This file is part of Slamdunk.
# 
# Slamdunk is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as
# published by the Free Software Foundation, either version 3 of the
# License, or (at your option) any later version.
# 
# Slamdunk is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Affero General Public License for more details.
# 
# You should have received a copy of the GNU Affero General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

# Load packages only from local Rslamdunk library 
libLoc = .libPaths()[grep("Rslamdunk",.libPaths())]

# Check if libraries are available, install otherwise
source(paste(libLoc,'/../checkLibraries.R',sep=""))

checkLib(libLoc)

library(getopt, lib.loc = libLoc)

spec = matrix(c(
  'help'      , 'h', 0, "logical","print the usage of the command",
  'fileTab', "f", 2,"character","tsv table of rate files",
  'outputPDF', "O", 2,"character","output pdf file name",
  'outputPCA', "P", 2,"character","output PCA transformations file name"
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


if ( is.null(opt$fileTab) ) stop("arg fileTab must be specified")
if ( is.null(opt$outputPDF) ) { opt$outputFile = "out.pdf" }

library(ggplot2, lib.loc = libLoc)

samples = read.table(opt$fileTab,stringsAsFactors=FALSE,col.names = c("sample","file"), comment.char = "")

if (nrow(samples) <= 1) {
	cat('# slamdunk PCA\n',  file=opt$outputPCA)
	cat(paste(samples$sample,0,"0\n",sep="\t"),append=TRUE,file=opt$outputPCA)
	#signal success and exit.
	q(status=0)
}

countsList = list()

for (i in 1:nrow(samples)) {
  curTab = read.delim(samples$file[i],stringsAsFactors=FALSE, comment.char="#")
  
  countsList[[samples$sample[i]]] = curTab$TcReadCount
  
}

countMatrix = do.call(cbind, countsList)

variances = apply(countMatrix, 1, var)

sel = order(variances, decreasing=TRUE)[seq_len(min(500, length(variances)))]

pca = prcomp(t(countMatrix[sel,]))

PoV = pca$sdev ^ 2 / sum(pca$sdev ^ 2)

plotTab = data.frame(sample = row.names(pca$x), PC1 = pca$x[,1], PC2 = pca$x[,2])

pdf(opt$outputPDF)

ggplot(plotTab, aes(x=PC1, y=PC2, color = sample)) + geom_point(size = 3) +
  xlab(paste("PC1 (", round(PoV[1],digits=2), " % variance)",sep="")) +
  ylab(paste("PC2 (", round(PoV[2],digits=2), " % variance)",sep="")) +
  theme(legend.position="bottom", legend.title=element_blank()) + ggtitle("Slamdunk PCA")

dev.off()

cat('# slamdunk PCA\n',  file=opt$outputPCA)
write.table(plotTab,file=opt$outputPCA,append=TRUE,quote=FALSE,sep="\t",row.names=FALSE,col.names=FALSE)

#signal success and exit.
q(status=0)		
