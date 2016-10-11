#!/usr/bin/env Rscript

# Load packages only from local Rslamdunk library 
libLoc = .libPaths()[grep("Rslamdunk",.libPaths())]

# Check if libraries are available, install otherwise
source(paste(libLoc,'/../checkLibraries.R',sep=""))

checkLib(libLoc)

library(getopt, lib.loc = libLoc)

spec = matrix(c(
  'help'      , 'h', 0, "logical","print the usage of the command",
  'simulated', "s", 2,"character","Comma seperated list of simulated files",
  'slamdunk', "f", 2,"character","Comma seperated lost of SlamDunk results",
  'output', "o", 2,"character","Output pdf",
  'conversionrate', "c", 2,"character","Simulated conversion rate"
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


if ( is.null(opt$simulated) ) stop("arg simulated must be specified")
if ( is.null(opt$slamdunk) ) stop("arg slamdunk must be specified")
if ( is.null(opt$output) ) stop("arg output must be specified")
if ( is.null(opt$conversionrate) ) { opt$conversionrate = 0.03 }

#simulatedFiles = "/project/ngs/philipp/slamdunk-analysis/simulation/simulation_1/pooja_UTR_annotation_examples_1_0min_utrsummary.csv,/project/ngs/philipp/slamdunk-analysis/simulation/simulation_1/pooja_UTR_annotation_examples_2_15min_utrsummary.csv,/project/ngs/philipp/slamdunk-analysis/simulation/simulation_1/pooja_UTR_annotation_examples_3_30min_utrsummary.csv,/project/ngs/philipp/slamdunk-analysis/simulation/simulation_1/pooja_UTR_annotation_examples_4_60min_utrsummary.csv,/project/ngs/philipp/slamdunk-analysis/simulation/simulation_1/pooja_UTR_annotation_examples_5_180min_utrsummary.csv,/project/ngs/philipp/slamdunk-analysis/simulation/simulation_1/pooja_UTR_annotation_examples_6_360min_utrsummary.csv,/project/ngs/philipp/slamdunk-analysis/simulation/simulation_1/pooja_UTR_annotation_examples_7_720min_utrsummary.csv,/project/ngs/philipp/slamdunk-analysis/simulation/simulation_1/pooja_UTR_annotation_examples_8_1440min_utrsummary.csv"
#slamDunkFiles = "/project/ngs/philipp/slamdunk-analysis/simulation/simulation_1/slamdunk/count/pooja_UTR_annotation_examples_1_0min_reads_slamdunk_mapped_filtered_tcount.csv,/project/ngs/philipp/slamdunk-analysis/simulation/simulation_1/slamdunk/count/pooja_UTR_annotation_examples_2_15min_reads_slamdunk_mapped_filtered_tcount.csv,/project/ngs/philipp/slamdunk-analysis/simulation/simulation_1/slamdunk/count/pooja_UTR_annotation_examples_3_30min_reads_slamdunk_mapped_filtered_tcount.csv,/project/ngs/philipp/slamdunk-analysis/simulation/simulation_1/slamdunk/count/pooja_UTR_annotation_examples_4_60min_reads_slamdunk_mapped_filtered_tcount.csv,/project/ngs/philipp/slamdunk-analysis/simulation/simulation_1/slamdunk/count/pooja_UTR_annotation_examples_5_180min_reads_slamdunk_mapped_filtered_tcount.csv,/project/ngs/philipp/slamdunk-analysis/simulation/simulation_1/slamdunk/count/pooja_UTR_annotation_examples_6_360min_reads_slamdunk_mapped_filtered_tcount.csv,/project/ngs/philipp/slamdunk-analysis/simulation/simulation_1/slamdunk/count/pooja_UTR_annotation_examples_7_720min_reads_slamdunk_mapped_filtered_tcount.csv,/project/ngs/philipp/slamdunk-analysis/simulation/simulation_1/slamdunk/count/pooja_UTR_annotation_examples_8_1440min_reads_slamdunk_mapped_filtered_tcount.csv"
#outputFile = "/project/ngs/philipp/slamdunk-analysis/simulation/simulation_1/eval/conversion_rate_per_gene_eval_plots.pdf"
#timesParameter = "0,15,30,60,180,360,720,1440"
#conversionRate = 0.03

simulatedFiles = opt$simulated
#simulatedFiles = "simulation_1/pooja_UTR_annotation_examples_sample_1_0min_utrsummary.csv,simulation_1/pooja_UTR_annotation_examples_sample_2_15min_utrsummary.csv,simulation_1/pooja_UTR_annotation_examples_sample_3_30min_utrsummary.csv,simulation_1/pooja_UTR_annotation_examples_sample_4_60min_utrsummary.csv,simulation_1/pooja_UTR_annotation_examples_sample_5_180min_utrsummary.csv,simulation_1/pooja_UTR_annotation_examples_sample_6_360min_utrsummary.csv,simulation_1/pooja_UTR_annotation_examples_sample_7_720min_utrsummary.csv,simulation_1/pooja_UTR_annotation_examples_sample_8_1440min_utrsummary.csv"
slamDunkFiles = opt$slamdunk
#slamDunkFiles = "simulation_1/slamdunk/count/pooja_UTR_annotation_examples_sample_1_0min_reads_slamdunk_mapped_filtered_tcount.csv,simulation_1/slamdunk/count/pooja_UTR_annotation_examples_sample_2_15min_reads_slamdunk_mapped_filtered_tcount.csv,simulation_1/slamdunk/count/pooja_UTR_annotation_examples_sample_3_30min_reads_slamdunk_mapped_filtered_tcount.csv,simulation_1/slamdunk/count/pooja_UTR_annotation_examples_sample_4_60min_reads_slamdunk_mapped_filtered_tcount.csv,simulation_1/slamdunk/count/pooja_UTR_annotation_examples_sample_5_180min_reads_slamdunk_mapped_filtered_tcount.csv,simulation_1/slamdunk/count/pooja_UTR_annotation_examples_sample_6_360min_reads_slamdunk_mapped_filtered_tcount.csv,simulation_1/slamdunk/count/pooja_UTR_annotation_examples_sample_7_720min_reads_slamdunk_mapped_filtered_tcount.csv,simulation_1/slamdunk/count/pooja_UTR_annotation_examples_sample_8_1440min_reads_slamdunk_mapped_filtered_tcount.csv"
filesSimulated = as.character(ordered(strsplit(simulatedFiles, ",")[[1]]))
filesSlamDunk = as.character(ordered(strsplit(slamDunkFiles, ",")[[1]]))
outputFile = opt$output
conversionRate = opt$conversionrate

pdf(outputFile, width = 15, height = 6)
for(timepoint in 1:length(filesSimulated)) {
  #timepoint = 5
  simDataFile = filesSimulated[timepoint]
  slamDunkFile = filesSlamDunk[timepoint]
  name = basename(simDataFile)
  
  simulation = read.table(simDataFile)
  colnames(simulation) = c("chr", "start", "stop", "name", "strand", "conversionRate", "readsCPM", "tCount", "tcCount", "readCount", "convertedReads", "multiMapCount")
  #simulation$coverage = (simulation$simulatedReads / (simulation$stop - simulation$start) * 50)
  simulation$convertedReadsRate = simulation$convertedReads / simulation$readCount
  
  slamdunk = read.table(slamDunkFile)
  colnames(slamdunk) = c("chr", "start", "stop", "name", "strand", "conversionRate", "readsCPM", "tCount", "tcCount", "readCount", "convertedReads", "multiMapCount")
  slamdunk$readsCPM_sim = simulation$readsCPM
  slamdunk$log2diff = log2((simulation$conversionRate + 0.0000001) / (slamdunk$conversionRate + 0.0000001))
  slamdunk$diff = (simulation$conversionRate - slamdunk$conversionRate)
  #slamdunk$convertedReadsRate = slamdunk$convertedReads / slamdunk$readCount
  #slamdunk$diffconvertedReadsRate = (simulation$convertedReads - slamdunk$convertedReads)
  #slamdunk$diffconvertedReadsRatediff = ((simulation$convertedReads / simulation$readCount) - (slamdunk$convertedReads / slamdunk$readCount))
  #plot(simulation$V11 ~ slamdunk$V6, xlim=c(0, 0.01), ylim=c(0, 0.01))
  
  par(mfrow=c(1,2))
  #yLim = max(abs(slamdunk$diffconvertedReadsRate))
  yLim = as.numeric(conversionRate)
  #boxplot(slamdunk$log2diff)
  slamDunkUniq = slamdunk[slamdunk$multiMapCount <= 0, ]
  slamDunkMulti = slamdunk[slamdunk$multiMapCount > 0, ]
  plot(slamDunkUniq$readsCPM_sim, slamDunkUniq$diff, main=name, pch=4, ylim=c(-yLim, yLim), ylab="conversion (sim) - conversion (slamdunk)", xlab="read counts per million")
  points(slamDunkMulti$readsCPM_sim, slamDunkMulti$diff, pch=4, col="red")
  abline(h=0, lty=2, col="grey")
  
  yLim = 4
  plot(slamDunkUniq$readsCPM_sim, slamDunkUniq$log2diff, main=name, pch=4, ylim=c(-yLim, yLim), ylab="log2(conversion (sim) / conversion (slamdunk))", xlab="read counts per million")
  points(slamDunkMulti$readsCPM_sim, slamDunkMulti$log2diff, pch=4, col="red")
  abline(h=0, lty=2, col="grey")
} 
dev.off()
