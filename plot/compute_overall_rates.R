#!/software/R/R-3.1.2/bin/Rscript
#!/usr/bin/Rscript
#!/software/R/R-3.1.2/bin/Rscript
#!/biosw/debian7-x86_64/R/3.0.2/bin/Rscript

#Script to overlap public database file 
# 
# Author: Tobias Neumann, Zuber group, Institute for Molecular Pathology
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

rates = read.table(opt$rateTab,stringsAsFactors=FALSE,col.names = c("sample","file"))

pdf(opt$outputFile)

plotList = list()

for (i in 1:nrow(rates)) {
	curTab = read.table(rates$file[i],stringsAsFactors=FALSE)
	
	#curTab = read.table(file,stringsAsFactors=FALSE)
	
	
	printTab = data.frame(rates=c(rep("AT",2),rep("AC",2),rep("AG",2),
					rep("TA",2),rep("TC",2),rep("TG",2),
					rep("CA",2),rep("CT",2),rep("CG",2),
					rep("GA",2),rep("GT",2),rep("GC",2)), strand = rep(c("+","-"),12),
			rate_percent = c(curTab["A","T"],curTab["A","t"],curTab["A","C"],curTab["A","c"],curTab["A","G"],curTab["A","g"],
					curTab["T","A"],curTab["T","a"],curTab["T","C"],curTab["T","c"],curTab["T","G"],curTab["T","g"],
					curTab["C","A"],curTab["C","a"],curTab["C","T"],curTab["C","t"],curTab["C","G"],curTab["C","g"],
					curTab["G","A"],curTab["G","a"],curTab["G","T"],curTab["G","t"],curTab["G","C"],curTab["G","c"])
	)
	

	fwdATot = max(1, sum(curTab["A",c("A", "C", "G", "T", "N")]))
	fwdCTot = max(1, sum(curTab["C",c("A", "C", "G", "T", "N")]))
	fwdGTot = max(1, sum(curTab["G",c("A", "C", "G", "T", "N")]))
	fwdTTot = max(1, sum(curTab["T",c("A", "C", "G", "T", "N")]))

	revATot = max(1, sum(curTab["A",c("a", "c", "g", "t", "n")]))
	revCTot = max(1, sum(curTab["C",c("a", "c", "g", "t", "n")]))
	revGTot = max(1, sum(curTab["G",c("a", "c", "g", "t", "n")]))
	revTTot = max(1, sum(curTab["T",c("a", "c", "g", "t", "n")]))

	total = c(rep(c(fwdATot, revATot), 3), rep(c(fwdTTot, revTTot), 3), rep(c(fwdCTot, revCTot), 3), rep(c(fwdGTot, revGTot), 3) )

	printTab$rate_percent = printTab$rate_percent / total * 100
	#printTab$rate_percent = printTab$rate_percent / sum(printTab$rate_percent) * 100
	
	curPlot = qplot(x=rates, y=rate_percent, fill=strand,data=printTab, geom="bar", stat="identity") + geom_text(aes(label = round(rate_percent,digits=2)), size = 1, hjust = 0.5, vjust = 1.5, position = "stack") + ylab("Rate percent %") + xlab(rates$sample[i]) +
			theme(text = element_text(size=6),axis.text.x = element_text(size=6)) 
	curPlot + xlim(0,35)	
	plotList[[length(plotList)+1]] <- curPlot + ylim(0.0,3.5)
}

do.call(grid.arrange,  plotList)

dev.off()

#signal success and exit.
q(status=0)		
