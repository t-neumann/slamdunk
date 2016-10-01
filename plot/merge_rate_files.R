#!/usr/bin/env Rscript
#
# Script to merge SlamDunk count files
# 
# Author: Philipp Rescheneder
# Email: philipp.rescheneder@univie.ac.at 
###############################################################################

library(getopt)

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
#times = as.numeric(strsplit(timesParameter, ",")[[1]])
#times = times / 60
times = strsplit(timesParameter, ",")[[1]]

mergeRates <- function(times, files, perRead) {
  mergedRates = data.frame()
  for(i in 1:length(times)) {
    time = times[i]
    print(time)
    simDataFile = files[i]
    simulation = read.table(simDataFile)
    colnames(simulation) = c("chr", "start", "stop", "name", "length", "strand", "conversionRate", "readsCPM", "tContent", "tCount", "tcCount", "readCount", "convertedReads", "multiMapCount")
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

perRead = F
slamDunkMergedRates = mergeRates(times, filesSlamDunk, perRead)

write.table(slamDunkMergedRates, outputFile, sep = "\t", quote = F, row.names = F, col.names = T)