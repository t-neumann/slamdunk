#!/usr/bin/env Rscript

# Plot overall conversion rates

###############################################################################
# Author: Tobias Neumann
# Email: tobias.neumann@imp.ac.at
###############################################################################

library(getopt)

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
  q(status=1);
}


if ( is.null(opt$rateTab) ) stop("arg rateTab must be specified")
if ( is.null(opt$outputFile) ) { opt$outputFile = "out.pdf" }

require(ggplot2)
require(gridExtra)
require(tidyr)

rates = read.table(opt$rateTab,stringsAsFactors=FALSE,col.names = c("sample","file"), comment.char = "")

pdf(opt$outputFile)

plotList = list()

for (i in 1:nrow(rates)) {
  curTab = read.delim(rates$file[i],stringsAsFactors=FALSE)
  
  plotTab = gather(curTab,"class","values",A_A,A_C,A_G,A_T,A_N,C_A,C_C,C_G,C_T,C_N,G_A,G_C,G_G,G_T,G_N,T_A,T_C,T_G,T_T,T_N,N_A,N_C,N_G,N_T,N_N)
  plotTab$values = plotTab$values / sum(plotTab$values) * 100
  plotTab = plotTab[!plotTab$class %in% c("A_A","T_T","G_G","C_C","N_N","N_A","N_C", "N_T","N_G", "A_N", "C_N", "T_N", "G_N"),]
  plotTab$highlight = "no"
  plotTab$highlight[plotTab$class == "T_C"] = "yes"
  plotTab$class = sub("_", ">", plotTab$class)
  
  curPlot = ggplot(plotTab, aes(x=class,y=values,fill=highlight,col=highlight)) + geom_boxplot() + xlab("") + ylab("Mutation rate per UTR[%]") +
  scale_fill_manual(values=c("black","red")) + scale_color_manual(values=c("black", "red")) + theme_classic() +  theme(axis.ticks.x = element_blank(), axis.line.y = element_line(color="black"), legend.position = "none")
  
  plotList[[length(plotList)+1]] <- curPlot + ylim(0.0,3.5)
}

do.call(grid.arrange,  plotList)

dev.off()

#signal success and exit.
q(status=0)		
