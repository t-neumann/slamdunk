#!/usr/bin/env Rscript

# Plot overall conversion rates per UTR

###############################################################################
# Author: Tobias Neumann
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
	q(status=1);
}


if ( is.null(opt$rateTab) ) stop("arg rateTab must be specified")
if ( is.null(opt$outputFile) ) { opt$outputFile = "out.pdf" }

library(ggplot2, lib.loc = libLoc)
library(gridExtra, lib.loc = libLoc)
library(dplyr, lib.loc = libLoc)
library(tidyr, lib.loc = libLoc)

rates = read.table(opt$rateTab,stringsAsFactors=FALSE,col.names = c("sample","file"), comment.char = "")

pdf(opt$outputFile)

plotList = list()

for (i in 1:nrow(rates)) {
	curTab = read.delim(rates$file[i],stringsAsFactors=FALSE,comment.char='#')
	
	plusTab = curTab %>% dplyr::filter(Strand == "+")
	minusTab = curTab %>% dplyr::filter(Strand == "-") %>%
			dplyr::select(A_A = T_T, G_G = C_C, C_C = G_G, T_T = A_A, A_C = T_G, A_G = T_C, A_T = T_A, C_A = G_T, C_G = G_C, C_T = G_A, G_A = C_T, G_C = C_G, G_T = C_A, T_A = A_T, T_C = A_G, T_G = A_C)
	
	#plusTab = plusTab %>% dplyr::select(-contains("N")) %>% dplyr::select(-one_of("Chr","Start")) %>% mutate(rowsum = rowSums(.)) %>% filter(rowsum > 0) %>% mutate_each(funs(. / rowsum)) %>% dplyr::select(-one_of(c("rowsum","A_A","C_C","G_G","T_T")))
	plusTab = plusTab %>% dplyr::select(-contains("N")) %>%
			dplyr::select(-one_of("Chr","Start")) %>%
			filter(rowSums(.) > 0) %>%
			mutate(Asum = A_A + A_C + A_G + A_T) %>%
			mutate(Csum = C_A + C_C + C_G + C_T) %>% 
			mutate(Gsum = G_A + G_C + G_G + G_T) %>%
			mutate(Tsum = T_A + T_C + T_G + T_T) %>%
			mutate_each(funs(. / Asum),matches("A_")) %>%
			mutate_each(funs(. / Csum),matches("C_")) %>%
			mutate_each(funs(. / Gsum),matches("G_")) %>%
			mutate_each(funs(. / Tsum),matches("T_")) %>%
			dplyr::select(-one_of(c("Asum","Csum","Gsum","Tsum"))) %>%
			mutate_each(funs(. * 100))
	
	minusTab = minusTab %>% filter(rowSums(.) > 0) %>%
			mutate(Asum = A_A + A_C + A_G + A_T) %>%
			mutate(Csum = C_A + C_C + C_G + C_T) %>% 
			mutate(Gsum = G_A + G_C + G_G + G_T) %>%
			mutate(Tsum = T_A + T_C + T_G + T_T) %>%
			mutate_each(funs(. / Asum),matches("A_")) %>%
			mutate_each(funs(. / Csum),matches("C_")) %>%
			mutate_each(funs(. / Gsum),matches("G_")) %>%
			mutate_each(funs(. / Tsum),matches("T_")) %>%
			dplyr::select(-one_of(c("Asum","Csum","Gsum","Tsum"))) %>%
			mutate_each(funs(. * 100))
	
	#minusTab = minusTab %>% mutate(rowsum = rowSums(.)) %>% filter(rowsum > 0) %>% mutate_each(funs(. / rowsum)) %>% dplyr::select(-rowsum)
	
	#plotTab = gather(curTab,"class","values",A_A,A_C,A_G,A_T,C_A,C_C,C_G,C_T,G_A,G_C,G_G,G_T,T_A,T_C,T_G,T_T)
	#plotTab = rbind(gather(plusTab,"class","values",A_C,A_G,A_T,C_A,C_G,C_T,G_A,G_C,G_T,T_A,T_C,T_G),
	#                gather(minusTab,"class","values",A_C,A_G,A_T,C_A,C_G,C_T,G_A,G_C,G_T,T_A,T_C,T_G)
	#)
	plotTab = rbind(plusTab, minusTab)
	plotTab = plotTab %>% dplyr::select(A_C,A_G,A_T,C_A,C_G,C_T,G_A,G_C,G_T,T_A,T_C,T_G)
	quantiles = lapply(plotTab, function(x) {
				return(quantile(x, na.rm=TRUE, p=0.75) + 1.5 * IQR(x, na.rm=TRUE))
			})
	
	#ymax = 0
	#for (base in names(quantiles)) {
	#  whiskerTop = quantiles[[base]][4] + 1.5 * (quantiles[[base]][4] - quantiles[[base]][2])
	#  print(whiskerTop)
	#  ymax = max(ymax, whiskerTop)
	#}
	#ymax = ceiling(ymax)
	ymax = ceiling(max(unlist(quantiles)))
	#plotTab$values = plotTab$values / sum(plotTab$values) * 100
	plotTab = gather(plotTab,"class","values",A_C,A_G,A_T,C_A,C_G,C_T,G_A,G_C,G_T,T_A,T_C,T_G)
	plotTab$highlight = "no"
	plotTab$highlight[plotTab$class == "T_C"] = "yes"
	plotTab$class = sub("_", ">", plotTab$class)
	plotTab$group = "A"
	plotTab$group[plotTab$class %in% c("C>A","C>G","C>T")] = "C"
	plotTab$group[plotTab$class %in% c("G>A","G>C","G>T")] = "G"
	plotTab$group[plotTab$class %in% c("T>A","T>C","T>G")] = "T"
	
	plotTab = plotTab[!is.na(plotTab$values),]
	
	#curPlot = ggplot(plotTab, aes(x=class,y=values,fill=highlight,col=highlight)) + stat_boxplot(geom ='errorbar') + geom_boxplot(outlier.shape = NA,lwd=0.8,fatten=2) + facet_grid(~group, scales="free", space="free") + xlab("") + ylab("Mutation rate per UTR base [%]") +
	#scale_fill_manual(values=c("white","white")) + scale_color_manual(values=c("black", "red")) + theme_classic() +  theme(axis.ticks.x = element_blank(), axis.line.y = element_line(color="black"), legend.position = "none")
	
	#plotList[[length(plotList)+1]] <- curPlot + ylim(0,ymax)
	
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
