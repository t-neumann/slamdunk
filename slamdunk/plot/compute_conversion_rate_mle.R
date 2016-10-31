#!/usr/bin/env Rscript
#
# 
# 
# Author: Philipp Rescheneder
# Email: philipp.rescheneder@univie.ac.at 
###############################################################################

# Load packages only from local Rslamdunk library 
#libLoc = .libPaths()[grep("Rslamdunk",.libPaths())]

# Check if libraries are available, install otherwise
#source(paste(libLoc,'/../checkLibraries.R',sep=""))

#checkLib(libLoc)

#library(getopt, lib.loc = libLoc)
library(getopt)
library(bbmle)

spec = matrix(c(
  'help'      , 'h', 0, "logical","print the usage of the command",
  'file', "f", 2,"character","",
  'rate', "r", 2,"character", "",
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

if ( is.null(opt$file) ) stop("arg file must be specified")
if ( is.null(opt$output) ) stop("arg output must be specified")
if ( is.null(opt$rate) ) stop("arg rate must be specified")

# a: percentage of converted transcripts
# b: conversion rate
LL <- function(a, b) {
  R = a * (( 1 - b ) ^ ( sample$n - sample$k )) * (b ^ sample$k) * choose(sample$n, sample$k) + ( 1 - a ) * as.numeric(sample$k == 0)
  -sum(log(R))
}

# Estimate a with fixed b
LL2 <- function(as) {
  b = estb
  rs = c()
  for(a in as)  {
    R = a * (( 1 - b ) ^ ( sample$n - sample$k )) * (b ^ sample$k) * choose(sample$n, sample$k) + ( 1 - a ) * as.numeric(sample$k == 0)
    rs = c(rs, -sum(log(R)))
  }
  rs
}

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



file = opt$file
#file = "/project/ngs/philipp/slamdunk-analysis/veronika/ngm-20161027/count-examples/34504_An312_wt-2n_mRNA-slamseq-autoquant_1h-R2.fq_slamdunk_mapped_filtered_tcount_perread.tsv,/project/ngs/philipp/slamdunk-analysis/veronika/ngm-20161027/count-examples/34504_An312_wt-2n_mRNA-slamseq-autoquant_1h-R2.fq_slamdunk_mapped_filtered_tcount_perread.tsv,/project/ngs/philipp/slamdunk-analysis/veronika/ngm-20161027/count-examples/34359_An312_wt-2n_mRNA-slamseq-autoquant_1h-R3.fq_slamdunk_mapped_filtered_tcount_perread.tsv"
files = as.character(ordered(strsplit(file, ",")[[1]]))

output = opt$output
#output = "/project/ngs/philipp/slamdunk-analysis/veronika/ngm-20161027/count-examples/34504_An312_wt-2n_mRNA-slamseq-autoquant_1h-R2.fq_slamdunk_mapped_filtered_tcount_mle.tsv"

estb = as.numeric(opt$rate)
#estb = 0.024

meta = readMeatInfo(files[1])
id = meta[1]
type = meta[3]
time = meta[4]

all = data.frame()
for(file in files) {
  part = read.table(file, header = T)
  all = rbind(all, part)
}
#head(all)

result = c()
names = as.character(unique(all$utr))

for(name in names) {
  #name = names[2]
  sample = all[all$utr == name, ]  
  
  fit = mle2(minuslogl = LL2, start = list(as = 0.29), method = "L-BFGS-B", lower = c(as = 0.000001), upper = c(as = 0.99))
  
  confinv = confint(fit)
  result = rbind(result, c(id, type, time, name, fit@coef[[1]], confinv[[1]], confinv[[2]]))
}

write.table(result, output, sep = "\t", quote = F, row.names = F, col.names = F)
