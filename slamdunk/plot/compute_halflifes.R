#!/usr/bin/env Rscript
#
# Script to compute half-lifes from SlamSeq data 
# 
# Author: Bhat Pooja, Philipp Rescheneder
# Email: pooja.bhat@imba.oeaw.ac.at, philipp.rescheneder@univie.ac.at 
###############################################################################

# Load packages only from local Rslamdunk library 
libLoc = .libPaths()[grep("Rslamdunk",.libPaths())]

# Check if libraries are available, install otherwise
source(paste(libLoc,'/../checkLibraries.R',sep=""))

checkLib(libLoc)

library(getopt, lib.loc = libLoc)

spec = matrix(c(
  'help'      , 'h', 0, "logical","print the usage of the command",
  'slamdunk', "f", 2,"character","Comma seperated list of SlamDunk results",
  'timepoints', "t", 2,"character","Comma seperated list of time points",
  'output', "o", 2,"character","Output tsv"
),ncol = 5,byrow=T)

opt = getopt(spec)

if ( !is.null(opt$help) || length(opt)==3 ) {
  #get the script name
  cmd = commandArgs(FALSE)
  self = strsplit(cmd[grep("--file",cmd)],"=")[[1]][2]
  cat(basename(self),": Create mismatch plots from rate tabs.\n\n")
  #print a friendly message and exit with a non-zero error code
  cat(getopt(spec,command = self,usage=T))
  q(status=1);
}

if ( is.null(opt$slamdunk) ) stop("arg slamdunk must be specified")
if ( is.null(opt$output) ) stop("arg output must be specified")
if ( is.null(opt$timepoints) ) stop("arg timepoints must be specified")

slamDunkFiles = opt$slamdunk
#slamDunkFiles = "/project/ngs/philipp/slamdunk-analysis/simulation/simulation_1/slamdunk/count/pooja_UTR_annotation_examples_1_0min_reads_slamdunk_mapped_filtered_tcount.csv,/project/ngs/philipp/slamdunk-analysis/simulation/simulation_1/slamdunk/count/pooja_UTR_annotation_examples_2_15min_reads_slamdunk_mapped_filtered_tcount.csv,/project/ngs/philipp/slamdunk-analysis/simulation/simulation_1/slamdunk/count/pooja_UTR_annotation_examples_3_30min_reads_slamdunk_mapped_filtered_tcount.csv,/project/ngs/philipp/slamdunk-analysis/simulation/simulation_1/slamdunk/count/pooja_UTR_annotation_examples_4_60min_reads_slamdunk_mapped_filtered_tcount.csv,/project/ngs/philipp/slamdunk-analysis/simulation/simulation_1/slamdunk/count/pooja_UTR_annotation_examples_5_180min_reads_slamdunk_mapped_filtered_tcount.csv,/project/ngs/philipp/slamdunk-analysis/simulation/simulation_1/slamdunk/count/pooja_UTR_annotation_examples_6_360min_reads_slamdunk_mapped_filtered_tcount.csv,/project/ngs/philipp/slamdunk-analysis/simulation/simulation_1/slamdunk/count/pooja_UTR_annotation_examples_7_720min_reads_slamdunk_mapped_filtered_tcount.csv,/project/ngs/philipp/slamdunk-analysis/simulation/simulation_1/slamdunk/count/pooja_UTR_annotation_examples_8_1440min_reads_slamdunk_mapped_filtered_tcount.csv"
filesSlamDunk = as.character(ordered(strsplit(slamDunkFiles, ",")[[1]]))
outputFile = opt$output
#outputFile = "/project/ngs/philipp/slamdunk-analysis/simulation/simulation_1/eval/halflife_per_gene_eval_plots.tsv"
#timesParameter = "0,15,30,60,180,360,720,1440"
timesParameter = opt$timepoints
times = as.numeric(strsplit(timesParameter, ",")[[1]])
times = times / 60


mergeRates <- function(times, files, perRead) {
  mergedRates = data.frame()
  for(i in 1:length(times)) {
    time = times[i]
    print(time)
    simDataFile = files[i]
    simulation = read.table(simDataFile)
    colnames(simulation) = c("chr", "start", "stop", "name", "strand", "conversionRate", "readsCPM", "tCount", "tcCount", "readCount", "convertedReads", "multiMapCount")
    if(nrow(mergedRates) == 0) {
      mergedRates = simulation[, c("chr", "start", "stop", "name", "strand")]
      mergedRates$avgReadsCPM = simulation$readsCPM
      mergedRates$avgMultimapper = simulation$multiMapCount
      if(perRead == TRUE) {
        mergedRates$conversionRate = simulation$convertedReads / simulation$readCount
      } else {
        mergedRates$conversionRate = simulation$conversionRate
      }
    } else {
      mergedRates$avgReadsCPM = mergedRates$avgReadsCPM + simulation$readsCPM
      mergedRates$avgMultimapper = mergedRates$avgMultimapper + simulation$multiMapCount
      if(perRead == TRUE) {
        mergedRates = cbind(mergedRates, simulation$convertedReads / simulation$readCount)
      } else {
        mergedRates = cbind(mergedRates, simulation$conversionRate)
      }
    }
  }
  colnames(mergedRates) = c("chr", "start", "stop", "name", "strand", "readsCPM", "multiMapCount", times)
  mergedRates$readsCPM = mergedRates$readsCPM / length(times)
  mergedRates$multiMapCount = mergedRates$multiMapCount / length(times)
  mergedRates
}

computeHalfLife <- function(rates, timepoints) {
  # Infere half life from data
  a_start<-max(rates) #param a is the y value when x=0
  k_start = log(2, base = exp(1))/5
  
  halfLifePred = NA
  C = NA
  k = NA
  
  tryCatch( {
    fit = nls(rates ~ a*(1-exp(-k*(timepoints))), start=list(a=a_start,k=k_start))
    halfLifePred = log(2, base = exp(1))/coef(fit)[2] * 60
    C = coef(fit)[1]
    k = coef(fit)[2]
  }, error=function(e){})
  summary(fit)
  
  RSS.p <- sum(residuals(fit)^2)
  TSS <- sum((rates - mean(rates))^2)
  rsquared = 1 - (RSS.p/TSS)
  
  c(halfLifePred, C, k, rsquared)
}

perRead = F
slamDunkMergedRates = mergeRates(times, filesSlamDunk, perRead)

halfLifeTable = data.frame()

for(utr in 1:nrow(slamDunkMergedRates)) {
  #utr = 8
  slamDunkMergedRates[utr,]
  pulseSlamDunk = data.frame(y = as.numeric(t(slamDunkMergedRates[utr, 8:(7 + length(times))])[,1]), x = times)
  
  result = computeHalfLife(pulseSlamDunk$y, pulseSlamDunk$x)
  #rates = pulseSlamDunk$y
  #timepoints = pulseSlamDunk$x
  halfLifeTable = rbind(halfLifeTable, cbind(slamDunkMergedRates[utr, c("chr", "start", "stop", "name", "strand", "readsCPM", "multiMapCount")], result[1]))
}  

colnames(halfLifeTable) = c("#chr", "start", "stop", "name", "strand", "readsCPM", "multiMapCount", "score")

write.table(halfLifeTable, outputFile, sep = "\t", quote = F, row.names = F, col.names = T)