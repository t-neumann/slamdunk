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
  'output', "o", 2,"character","Output tsv",
  'alternativecounting', "a", 2,"character","Use alternative counting not percentage of T->C reads"
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
if ( is.null(opt$alternativecounting) ) { opt$alternativecounting = 0 }

slamDunkFiles = opt$slamdunk
#slamDunkFiles = "/project/ngs/philipp/slamdunk-analysis/simulation/simulation_1/slamdunk/count/pooja_UTR_annotation_examples_1_0min_reads_slamdunk_mapped_filtered_tcount.csv,/project/ngs/philipp/slamdunk-analysis/simulation/simulation_1/slamdunk/count/pooja_UTR_annotation_examples_2_15min_reads_slamdunk_mapped_filtered_tcount.csv,/project/ngs/philipp/slamdunk-analysis/simulation/simulation_1/slamdunk/count/pooja_UTR_annotation_examples_3_30min_reads_slamdunk_mapped_filtered_tcount.csv,/project/ngs/philipp/slamdunk-analysis/simulation/simulation_1/slamdunk/count/pooja_UTR_annotation_examples_4_60min_reads_slamdunk_mapped_filtered_tcount.csv,/project/ngs/philipp/slamdunk-analysis/simulation/simulation_1/slamdunk/count/pooja_UTR_annotation_examples_5_180min_reads_slamdunk_mapped_filtered_tcount.csv,/project/ngs/philipp/slamdunk-analysis/simulation/simulation_1/slamdunk/count/pooja_UTR_annotation_examples_6_360min_reads_slamdunk_mapped_filtered_tcount.csv,/project/ngs/philipp/slamdunk-analysis/simulation/simulation_1/slamdunk/count/pooja_UTR_annotation_examples_7_720min_reads_slamdunk_mapped_filtered_tcount.csv,/project/ngs/philipp/slamdunk-analysis/simulation/simulation_1/slamdunk/count/pooja_UTR_annotation_examples_8_1440min_reads_slamdunk_mapped_filtered_tcount.csv"
#slamDunkFiles = "/project/ngs/philipp/slamdunk-analysis/simulation/simulation_1/slamdunk/count/pooja_UTR_annotation_examples_7_720min_reads_slamdunk_mapped_filtered_tcount.tsv,/project/ngs/philipp/slamdunk-analysis/simulation/simulation_1/slamdunk/count/pooja_UTR_annotation_examples_1_0min_reads_slamdunk_mapped_filtered_tcount.tsv,/project/ngs/philipp/slamdunk-analysis/simulation/simulation_1/slamdunk/count/pooja_UTR_annotation_examples_6_360min_reads_slamdunk_mapped_filtered_tcount.tsv,/project/ngs/philipp/slamdunk-analysis/simulation/simulation_1/slamdunk/count/pooja_UTR_annotation_examples_2_15min_reads_slamdunk_mapped_filtered_tcount.tsv,/project/ngs/philipp/slamdunk-analysis/simulation/simulation_1/slamdunk/count/pooja_UTR_annotation_examples_8_1440min_reads_slamdunk_mapped_filtered_tcount.tsv,/project/ngs/philipp/slamdunk-analysis/simulation/simulation_1/slamdunk/count/pooja_UTR_annotation_examples_3_30min_reads_slamdunk_mapped_filtered_tcount.tsv,/project/ngs/philipp/slamdunk-analysis/simulation/simulation_1/slamdunk/count/pooja_UTR_annotation_examples_5_180min_reads_slamdunk_mapped_filtered_tcount.tsv,/project/ngs/philipp/slamdunk-analysis/simulation/simulation_1/slamdunk/count/pooja_UTR_annotation_examples_4_60min_reads_slamdunk_mapped_filtered_tcount.tsv"
filesSlamDunk = as.character(ordered(strsplit(slamDunkFiles, ",")[[1]]))
outputFile = opt$output
#outputFile = "/project/ngs/philipp/slamdunk-analysis/simulation/simulation_1/eval/halflife_per_gene_eval_plots.tsv"
perRead = T
if(opt$alternativecounting > 0) {
  cat("Using alternative counting\n")
  perRead = F
}

readMeatInfo <- function(fileName) {
  #fileName = filesSlamDunk[1]
  sampleInfo = read.table(fileName, nrows = 1, comment.char = "")
  version = paste(lapply(sampleInfo[1,1:3], as.character), collapse = '\t')
  sampleName = as.character(sampleInfo[1, ]$V6)
  sampleType = as.character(sampleInfo[1, ]$V7)
  sampleTime = as.numeric(sampleInfo[1, ]$V8)
  sampleInfo = read.table(fileName, nrows = 1, skip = 1, comment.char = "")
  annotationMD5 = as.character(sampleInfo[1, ]$V3)
  annotationName = as.character(sampleInfo[1, ]$V2)
  c(sampleName, sampleType, sampleTime, annotationName, annotationMD5, version)  
}
  
sampleNumber = length(filesSlamDunk)
mergedRates = data.frame()

annotationName = ""
annotationMD5 = ""
version = ""

# Merge rates from all samples
for(i in 1:length(filesSlamDunk)) {
  #i = 1
  file = filesSlamDunk[i]  
  meta = readMeatInfo(file)
  sampleName = meta[1]
  
  if(i == 1) {
    version = meta[6]
    annotationName = meta[4]
    annotationMD5 = meta[5]
  } else {
    if(annotationMD5 != meta[5]) {
      
    }
  }
  
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
  if(perRead == T) {
    mergedRates[,sampleName] = data$TcReadCount / data$ReadCount
    mergedRates[data$ReadCount == 0,sampleName] = 0
  } else {
    mergedRates[,sampleName] = data$ConversionRate
  }
}
# compute average CPM and multimapper per UTR
mergedRates$avgReadsCPM = mergedRates$avgReadsCPM / sampleNumber
mergedRates$avgMultimapper = mergedRates$avgMultimapper / sampleNumber
mergedRates$avgTcontent = mergedRates$avgTcontent / sampleNumber
mergedRates$avgCoverageOnTs = mergedRates$avgCoverageOnTs / sampleNumber

# Sort columns by sample name
#colNumber = length(colnames(mergedRates))
#firstSampleColumn = (colNumber - sampleNumber + 1)
#sampleNames = colnames(mergedRates)[firstSampleColumn:colNumber]
#sampleColumnOrder = order(sampleNames)
#mergedRates = mergedRates[, c(1:(firstSampleColumn - 1), (sampleColumnOrder + firstSampleColumn - 1))]

#head(mergedRates)

# Write to output file
con <- file(outputFile, open="wt") 
writeLines(version, con)
writeLines(paste0("#Annotation:\t", annotationName, "\t", annotationMD5), con) 
write.table(mergedRates, con, sep = "\t", quote = F, row.names = F, col.names = T)
close(con) 