#!/usr/bin/env Rscript
#
# Script to merge SlamDunk count files
# 
# Author: Philipp Rescheneder
# Email: philipp.rescheneder@univie.ac.at 
###############################################################################

# Load packages only from local Rslamdunk library 
libLoc = .libPaths()[grep("Rslamdunk",.libPaths())]

# Check if libraries are available, install otherwise
source(paste(libLoc,'/../checkLibraries.R',sep=""))

checkLib(libLoc)

library(getopt, lib.loc = libLoc)

#library(getopt)

spec = matrix(c(
  'help'      , 'h', 0, "logical","print the usage of the command",
  'slamdunk', "f", 2,"character","Comma seperated list of SlamDunk results",
  'output', "o", 2,"character","Output tsv",
  'column', "c", 2,"character","Column or Expression used to summarize files",
  'columnname', "n", 2,"character","Index of meta data field to use as column name"
),ncol = 5,byrow=T)

opt = getopt(spec)

if ( !is.null(opt$help) || length(opt)==2 ) {
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
if ( is.null(opt$column) ) { opt$column = "TcReadCount / ReadCount" }
if ( is.null(opt$column) ) { opt$columnname = 2 }

slamDunkFiles = opt$slamdunk
#slamDunkFiles = "/project/ngs/philipp/slamdunk-analysis/simulation/simulation_1/slamdunk/count/pooja_UTR_annotation_examples_7_720min_reads_slamdunk_mapped_filtered_tcount.tsv,/project/ngs/philipp/slamdunk-analysis/simulation/simulation_1/slamdunk/count/pooja_UTR_annotation_examples_1_0min_reads_slamdunk_mapped_filtered_tcount.tsv,/project/ngs/philipp/slamdunk-analysis/simulation/simulation_1/slamdunk/count/pooja_UTR_annotation_examples_6_360min_reads_slamdunk_mapped_filtered_tcount.tsv,/project/ngs/philipp/slamdunk-analysis/simulation/simulation_1/slamdunk/count/pooja_UTR_annotation_examples_2_15min_reads_slamdunk_mapped_filtered_tcount.tsv,/project/ngs/philipp/slamdunk-analysis/simulation/simulation_1/slamdunk/count/pooja_UTR_annotation_examples_8_1440min_reads_slamdunk_mapped_filtered_tcount.tsv,/project/ngs/philipp/slamdunk-analysis/simulation/simulation_1/slamdunk/count/pooja_UTR_annotation_examples_3_30min_reads_slamdunk_mapped_filtered_tcount.tsv,/project/ngs/philipp/slamdunk-analysis/simulation/simulation_1/slamdunk/count/pooja_UTR_annotation_examples_5_180min_reads_slamdunk_mapped_filtered_tcount.tsv,/project/ngs/philipp/slamdunk-analysis/simulation/simulation_1/slamdunk/count/pooja_UTR_annotation_examples_4_60min_reads_slamdunk_mapped_filtered_tcount.tsv"
#slamDunkFiles = "ngm-20161027/count/34330_An312_wt-2n_mRNA-slamseq-autoquant_0h-R1.fq_slamdunk_mapped_filtered_tcount.tsv,ngm-20161027/count/34331_An312_wt-2n_mRNA-slamseq-autoquant_0.25h-R1.fq_slamdunk_mapped_filtered_tcount.tsv,ngm-20161027/count/34332_An312_wt-2n_mRNA-slamseq-autoquant_0.5h-R1.fq_slamdunk_mapped_filtered_tcount.tsv,ngm-20161027/count/34333_An312_wt-2n_mRNA-slamseq-autoquant_1h-R1.fq_slamdunk_mapped_filtered_tcount.tsv,ngm-20161027/count/34334_An312_wt-2n_mRNA-slamseq-autoquant_3h-R1.fq_slamdunk_mapped_filtered_tcount.tsv,ngm-20161027/count/34335_An312_wt-2n_mRNA-slamseq-autoquant_6h-R1.fq_slamdunk_mapped_filtered_tcount.tsv,ngm-20161027/count/34336_An312_wt-2n_mRNA-slamseq-autoquant_12h-R1.fq_slamdunk_mapped_filtered_tcount.tsv,ngm-20161027/count/34343_An312_wt-2n_mRNA-slamseq-autoquant_0h-R2.fq_slamdunk_mapped_filtered_tcount.tsv,ngm-20161027/count/34344_An312_wt-2n_mRNA-slamseq-autoquant_0.25h-R2.fq_slamdunk_mapped_filtered_tcount.tsv,ngm-20161027/count/34347_An312_wt-2n_mRNA-slamseq-autoquant_3h-R2.fq_slamdunk_mapped_filtered_tcount.tsv,ngm-20161027/count/34348_An312_wt-2n_mRNA-slamseq-autoquant_6h-R2.fq_slamdunk_mapped_filtered_tcount.tsv,ngm-20161027/count/34349_An312_wt-2n_mRNA-slamseq-autoquant_12h-R2.fq_slamdunk_mapped_filtered_tcount.tsv,ngm-20161027/count/34354_An312_wt-2n_mRNA-slamseq-autoquant_24h-R2.fq_slamdunk_mapped_filtered_tcount.tsv,ngm-20161027/count/34356_An312_wt-2n_mRNA-slamseq-autoquant_0h-R3.fq_slamdunk_mapped_filtered_tcount.tsv,ngm-20161027/count/34357_An312_wt-2n_mRNA-slamseq-autoquant_0.25h-R3.fq_slamdunk_mapped_filtered_tcount.tsv,ngm-20161027/count/34358_An312_wt-2n_mRNA-slamseq-autoquant_0.5h-R3.fq_slamdunk_mapped_filtered_tcount.tsv,ngm-20161027/count/34359_An312_wt-2n_mRNA-slamseq-autoquant_1h-R3.fq_slamdunk_mapped_filtered_tcount.tsv,ngm-20161027/count/34361_An312_wt-2n_mRNA-slamseq-autoquant_6h-R3.fq_slamdunk_mapped_filtered_tcount.tsv,ngm-20161027/count/34362_An312_wt-2n_mRNA-slamseq-autoquant_12h-R3.fq_slamdunk_mapped_filtered_tcount.tsv,ngm-20161027/count/34367_An312_wt-2n_mRNA-slamseq-autoquant_24h-R3.fq_slamdunk_mapped_filtered_tcount.tsv,ngm-20161027/count/34503_An312_wt-2n_mRNA-slamseq-autoquant_0.5h-R2.fq_slamdunk_mapped_filtered_tcount.tsv,ngm-20161027/count/34504_An312_wt-2n_mRNA-slamseq-autoquant_1h-R2.fq_slamdunk_mapped_filtered_tcount.tsv,ngm-20161027/count/34505_An312_wt-2n_mRNA-slamseq-autoquant_3h-R3.fq_slamdunk_mapped_filtered_tcount.tsv,ngm-20161027/count/34506_An312_wt-2n_mRNA-slamseq-autoquant_24h-R1.fq_slamdunk_mapped_filtered_tcount.tsv,ngm-20161027/count/34507_An312_wt-2n_mRNA-slamseq-autoquant_24h-star-R1.fq_slamdunk_mapped_filtered_tcount.tsv,ngm-20161027/count/34508_An312_wt-2n_mRNA-slamseq-autoquant_24h-star-R2.fq_slamdunk_mapped_filtered_tcount.tsv,ngm-20161027/count/34509_An312_wt-2n_mRNA-slamseq-autoquant_24h-star-R3.fq_slamdunk_mapped_filtered_tcount.tsv,ngm-20161027/count/34510_An312_wt-2n_mRNA-slamseq-autoquant_24h-star-0.25h-R1.fq_slamdunk_mapped_filtered_tcount.tsv,ngm-20161027/count/34511_An312_wt-2n_mRNA-slamseq-autoquant_24h-star-0.25h-R2.fq_slamdunk_mapped_filtered_tcount.tsv,ngm-20161027/count/34512_An312_wt-2n_mRNA-slamseq-autoquant_24h-star-0.25h-R3.fq_slamdunk_mapped_filtered_tcount.tsv,ngm-20161027/count/34513_An312_wt-2n_mRNA-slamseq-autoquant_24h-star-0.5h-R1.fq_slamdunk_mapped_filtered_tcount.tsv,ngm-20161027/count/34514_An312_wt-2n_mRNA-slamseq-autoquant_24h-star-0.5h-R2.fq_slamdunk_mapped_filtered_tcount.tsv,ngm-20161027/count/34515_An312_wt-2n_mRNA-slamseq-autoquant_24h-star-0.5h-R3.fq_slamdunk_mapped_filtered_tcount.tsv,ngm-20161027/count/34516_An312_wt-2n_mRNA-slamseq-autoquant_24h-star-1h-R1.fq_slamdunk_mapped_filtered_tcount.tsv,ngm-20161027/count/34517_An312_wt-2n_mRNA-slamseq-autoquant_24h-star-1h-R2.fq_slamdunk_mapped_filtered_tcount.tsv,ngm-20161027/count/34518_An312_wt-2n_mRNA-slamseq-autoquant_24h-star-1h-R3.fq_slamdunk_mapped_filtered_tcount.tsv,ngm-20161027/count/34519_An312_wt-2n_mRNA-slamseq-autoquant_24h-star-3h-R1.fq_slamdunk_mapped_filtered_tcount.tsv,ngm-20161027/count/34520_An312_wt-2n_mRNA-slamseq-autoquant_24h-star-3h-R2.fq_slamdunk_mapped_filtered_tcount.tsv,ngm-20161027/count/34521_An312_wt-2n_mRNA-slamseq-autoquant_24h-star-3h-R3.fq_slamdunk_mapped_filtered_tcount.tsv,ngm-20161027/count/34522_An312_wt-2n_mRNA-slamseq-autoquant_24h-star-6h-R1.fq_slamdunk_mapped_filtered_tcount.tsv,ngm-20161027/count/34523_An312_wt-2n_mRNA-slamseq-autoquant_24h-star-6h-R2.fq_slamdunk_mapped_filtered_tcount.tsv,ngm-20161027/count/34524_An312_wt-2n_mRNA-slamseq-autoquant_24h-star-6h-R3.fq_slamdunk_mapped_filtered_tcount.tsv,ngm-20161027/count/34525_An312_wt-2n_mRNA-slamseq-autoquant_24h-star-12h-R1.fq_slamdunk_mapped_filtered_tcount.tsv,ngm-20161027/count/34526_An312_wt-2n_mRNA-slamseq-autoquant_24h-star-12h-R2.fq_slamdunk_mapped_filtered_tcount.tsv,ngm-20161027/count/34527_An312_wt-2n_mRNA-slamseq-autoquant_24h-star-12h-R3.fq_slamdunk_mapped_filtered_tcount.tsv,ngm-20161027/count/34528_An312_wt-2n_mRNA-slamseq-autoquant_24h-star-24h-R1.fq_slamdunk_mapped_filtered_tcount.tsv,ngm-20161027/count/34529_An312_wt-2n_mRNA-slamseq-autoquant_24h-star-24h-R2.fq_slamdunk_mapped_filtered_tcount.tsv,ngm-20161027/count/34530_An312_wt-2n_mRNA-slamseq-autoquant_24h-star-24h-R3.fq_slamdunk_mapped_filtered_tcount.tsv"
filesSlamDunk = as.character(ordered(strsplit(slamDunkFiles, ",")[[1]]))
outputFile = opt$output
#outputFile = "/project/ngs/philipp/slamdunk-analysis/simulation/simulation_1/eval/halflife_per_gene_eval_plots.tsv"
evalExpression = opt$column
#evalExpression = "TcReadCount / ReadCount"
columnName = as.integer(opt$columnname)
#columnName = 2

readMeatInfo <- function(fileName) {
  #fileName = filesSlamDunk[1]
  sampleInfo = read.table(fileName, nrows = 1, comment.char = "")
  version = paste(lapply(sampleInfo[1,1:3], as.character), collapse = '\t')
  sampleID = as.character(sampleInfo[1, ]$V7)
  sampleName = as.character(sampleInfo[1, ]$V6)
  sampleType = as.character(sampleInfo[1, ]$V8)
  sampleTime = as.numeric(sampleInfo[1, ]$V9)
  sampleInfo = read.table(fileName, nrows = 1, skip = 1, comment.char = "")
  annotationMD5 = as.character(sampleInfo[1, ]$V3)
  annotationName = as.character(sampleInfo[1, ]$V2)
  c(sampleID, sampleName, sampleType, sampleTime, annotationName, annotationMD5, version)  
}
  
sampleNumber = length(filesSlamDunk)
mergedRates = data.frame()

annotationName = ""
annotationMD5 = ""
version = ""
IDs = c()

# Merge rates from all samples
for(i in 1:length(filesSlamDunk)) {
  #i = 1
  file = filesSlamDunk[i]  
  meta = readMeatInfo(file)
  sampleName = meta[columnName]
  
  if(i == 1) {
    version = meta[7]
    annotationName = meta[5]
    annotationMD5 = meta[6]
  } else {
    if(annotationMD5 != meta[6]) {
      
    }
  }
  
  IDs = c(IDs, as.numeric(meta[1]))
  data = read.table(file, header = T)
  if(i == 1) {
    mergedRates = data[, c(1:6)]
    mergedRates$avgReadsCPM = data$ReadsCPM
    mergedRates$avgMultimapper = data$multimapCount
    mergedRates$avgTcontent = data$Tcontent
    mergedRates$avgCoverageOnTs = data$CoverageOnTs
  } else {
    mergedRates$avgReadsCPM = mergedRates$avgReadsCPM + data$ReadsCPM
    mergedRates$avgMultimapper = mergedRates$avgMultimapper + data$multimapCount
    mergedRates$avgTcontent = mergedRates$avgTcontent + data$Tcontent
    mergedRates$avgCoverageOnTs = mergedRates$avgCoverageOnTs + data$CoverageOnTs
  }
  #if(perRead == T) {
  attach(data)
  #mergedRates[,sampleName] = data$TcReadCount / data$ReadCount
  mergedRates[,sampleName] = eval(parse(text=evalExpression))
  detach(data)
  mergedRates[data$ReadCount == 0,sampleName] = 0
  #} else {
  #  mergedRates[,sampleName] = data$ConversionRate
  #}
}
# compute average CPM and multimapper per UTR
mergedRates$avgReadsCPM = mergedRates$avgReadsCPM / sampleNumber
mergedRates$avgMultimapper = mergedRates$avgMultimapper / sampleNumber
mergedRates$avgTcontent = mergedRates$avgTcontent / sampleNumber
mergedRates$avgCoverageOnTs = mergedRates$avgCoverageOnTs / sampleNumber

#head(mergedRates)
# Sort columns by sample name
colNumber = length(colnames(mergedRates))
firstSampleColumn = (colNumber - sampleNumber + 1)
sampleNames = colnames(mergedRates)[firstSampleColumn:colNumber]
sampleColumnOrder = order(IDs)
mergedRates = mergedRates[, c(1:(firstSampleColumn - 1), (sampleColumnOrder + firstSampleColumn - 1))]

#head(mergedRates)

# Write to output file
con <- file(outputFile, open="wt") 
writeLines(version, con)
writeLines(paste0("#Annotation:\t", annotationName, "\t", annotationMD5), con)
writeLines(paste0("#Expression:\t", evalExpression), con)
write.table(mergedRates, con, sep = "\t", quote = F, row.names = F, col.names = T)
close(con) 