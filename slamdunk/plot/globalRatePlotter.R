#!/usr/bin/env Rscript

# Plot overall conversion rates per UTR

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

library(ggplot2)
library(gridExtra)

rates = read.table(opt$rateTab,stringsAsFactors=FALSE,col.names = c("sample","file"), comment.char = "")

pdf(opt$outputFile)

plotList = list()

for (i in 1:nrow(rates)) {
	curTab = read.delim(rates$file[i],stringsAsFactors=FALSE,comment.char='#')

	plusTab = curTab[curTab$Strand == "+",]
	minusTab = curTab[curTab$Strand == "-",]

	# "Name"      "Chr"       "Start"     "End"       "Strand"    "ReadCount"
	# "A_A"       "A_C"       "A_G"       "A_T"       "A_N"       "C_A"
	# "C_C"       "C_G"       "C_T"       "C_N"       "G_A"       "G_C"
	# "G_G"       "G_T"       "G_N"       "T_A"       "T_C"       "T_G"
	# "T_T"       "T_N"       "N_A"       "N_C"       "N_G"       "N_T"
	# "N_N"

	names(minusTab) = c("Name", "Chr", "Start", "End", "Strand", "ReadCount",
			"T_T", "T_G", "T_C", "T_A", "NNN", "G_T",
			"G_G", "G_C", "G_A", "NNN", "C_T", "C_G",
			"C_C", "C_A", "NNN", "A_T", "A_G", "A_C",
			"A_A", "NNN", "NNN", "NNN", "NNN", "NNN",
			"NNN")

	plusTab = plusTab[,c(1,grep("N",names(plusTab),invert=TRUE))]
	minusTab = minusTab[,grep("NNN",names(minusTab),invert=TRUE)]

	plusTab = plusTab[,c(-1,-2,-3,-4,-5,-6)]
	plusTab = plusTab[rowSums(plusTab) > 0,]

	plusTab$Asum = plusTab$A_A + plusTab$A_C + plusTab$A_G + plusTab$A_T
	plusTab$Csum = plusTab$C_A + plusTab$C_C + plusTab$C_G + plusTab$C_T
	plusTab$Gsum = plusTab$G_A + plusTab$G_C + plusTab$G_G + plusTab$G_T
	plusTab$Tsum = plusTab$T_A + plusTab$T_C + plusTab$T_G + plusTab$T_T

	plusTab$A_A = plusTab$A_A / plusTab$Asum
	plusTab$A_C = plusTab$A_C / plusTab$Asum
	plusTab$A_G = plusTab$A_G / plusTab$Asum
	plusTab$A_T = plusTab$A_T / plusTab$Asum

	plusTab$C_A = plusTab$C_A / plusTab$Csum
	plusTab$C_C = plusTab$C_C / plusTab$Csum
	plusTab$C_G = plusTab$C_G / plusTab$Csum
	plusTab$C_T = plusTab$C_T / plusTab$Csum

	plusTab$G_A = plusTab$G_A / plusTab$Gsum
	plusTab$G_C = plusTab$G_C / plusTab$Gsum
	plusTab$G_G = plusTab$G_G / plusTab$Gsum
	plusTab$G_T = plusTab$G_T / plusTab$Gsum

	plusTab$T_A = plusTab$T_A / plusTab$Tsum
	plusTab$T_C = plusTab$T_C / plusTab$Tsum
	plusTab$T_G = plusTab$T_G / plusTab$Tsum
	plusTab$T_T = plusTab$T_T / plusTab$Tsum

	plusTab = plusTab[,grep("sum",names(plusTab),invert=TRUE)]

	plusTab = plusTab * 100

	minusTab = minusTab[,c(-1,-2,-3,-4,-5,-6)]
	minusTab = minusTab[rowSums(minusTab) > 0,]

	minusTab$Asum = minusTab$A_A + minusTab$A_C + minusTab$A_G + minusTab$A_T
	minusTab$Csum = minusTab$C_A + minusTab$C_C + minusTab$C_G + minusTab$C_T
	minusTab$Gsum = minusTab$G_A + minusTab$G_C + minusTab$G_G + minusTab$G_T
	minusTab$Tsum = minusTab$T_A + minusTab$T_C + minusTab$T_G + minusTab$T_T

	minusTab$A_A = minusTab$A_A / minusTab$Asum
	minusTab$A_C = minusTab$A_C / minusTab$Asum
	minusTab$A_G = minusTab$A_G / minusTab$Asum
	minusTab$A_T = minusTab$A_T / minusTab$Asum

	minusTab$C_A = minusTab$C_A / minusTab$Csum
	minusTab$C_C = minusTab$C_C / minusTab$Csum
	minusTab$C_G = minusTab$C_G / minusTab$Csum
	minusTab$C_T = minusTab$C_T / minusTab$Csum

	minusTab$G_A = minusTab$G_A / minusTab$Gsum
	minusTab$G_C = minusTab$G_C / minusTab$Gsum
	minusTab$G_G = minusTab$G_G / minusTab$Gsum
	minusTab$G_T = minusTab$G_T / minusTab$Gsum

	minusTab$T_A = minusTab$T_A / minusTab$Tsum
	minusTab$T_C = minusTab$T_C / minusTab$Tsum
	minusTab$T_G = minusTab$T_G / minusTab$Tsum
	minusTab$T_T = minusTab$T_T / minusTab$Tsum

	minusTab = minusTab[,grep("sum",names(minusTab),invert=TRUE)]

	minusTab = minusTab * 100

	plotTab = rbind(plusTab, minusTab)

	plotTab = plotTab[,c("A_C","A_G","A_T","C_A","C_G","C_T","G_A","G_C","G_T","T_A","T_C","T_G")]
	quantiles = lapply(plotTab, function(x) {
				return(quantile(x, na.rm=TRUE, p=0.75) + 1.5 * IQR(x, na.rm=TRUE))
			})

	ymax = ceiling(max(unlist(quantiles)))

	plotTab = rbind(
			data.frame(class = "A_C", values = plotTab$A_C),
			data.frame(class = "A_G", values = plotTab$A_G),
			data.frame(class = "A_T", values = plotTab$A_T),
			data.frame(class = "C_A", values = plotTab$C_A),
			data.frame(class = "C_G", values = plotTab$C_G),
			data.frame(class = "C_T", values = plotTab$C_T),
			data.frame(class = "G_A", values = plotTab$G_A),
			data.frame(class = "G_C", values = plotTab$G_C),
			data.frame(class = "G_T", values = plotTab$G_T),
			data.frame(class = "T_A", values = plotTab$T_A),
			data.frame(class = "T_C", values = plotTab$T_C),
			data.frame(class = "T_G", values = plotTab$T_G)
	)

	plotTab$highlight = "no"
	plotTab$highlight[plotTab$class == "T_C"] = "yes"
	plotTab$class = sub("_", ">", plotTab$class)
	plotTab$group = "A"
	plotTab$group[plotTab$class %in% c("C>A","C>G","C>T")] = "C"
	plotTab$group[plotTab$class %in% c("G>A","G>C","G>T")] = "G"
	plotTab$group[plotTab$class %in% c("T>A","T>C","T>G")] = "T"

	plotTab = plotTab[!is.na(plotTab$values),]

	curPlot = ggplot(plotTab, aes(x=class,y=values,fill=highlight,col=highlight)) + stat_boxplot(geom ='errorbar') + geom_boxplot(outlier.shape = NA,lwd=0.8,fatten=2) + facet_grid(~group, scales="free", space="free") + xlab("") + ylab("Mutation rate per UTR base [%]") +
			scale_fill_manual(values=c("white","white")) + scale_color_manual(values=c("black", "red")) + theme(axis.ticks.x = element_blank(), legend.position = "none")

	plotList[[length(plotList)+1]] <- curPlot + ylim(0,ymax) + ggtitle(rates$sample[i])

	anovaTest = aov(values ~ class, data = plotTab)
	print(paste("Sample: ",rates$sample[i],sep=""))
	print(TukeyHSD(x = anovaTest, 'class', conf.level = 0.95)$class)

}

do.call(grid.arrange,  plotList)

dev.off()

#signal success and exit.
q(status=0)
