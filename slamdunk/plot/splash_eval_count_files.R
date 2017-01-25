#!/usr/bin/env Rscript
# Script to evaluate Slamdunk count results
# 
# Author: Philipp Rescheneder
# Email: philipp.rescheneder@gmail.com
###############################################################################

# Load packages only from local Rslamdunk library 
libLoc = .libPaths()[grep("Rslamdunk",.libPaths())]

# Check if libraries are available, install otherwise
source(paste(libLoc,'/../checkLibraries.R',sep=""))

checkLib(libLoc)

library(getopt, lib.loc = libLoc)


spec = matrix(c(
  'help'      , 'h', 0, "logical","print the usage of the command",
  'simulated', "s", 2,"character","Summarized count file",
  'slamdunk', "d", 2,"character","Summarized count file",
  'output', "o", 2,"character","Output pdf"
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


if ( is.null(opt$simulated) ) stop("arg simulated must be specified")
if ( is.null(opt$slamdunk) ) stop("arg slamdunk must be specified")
if ( is.null(opt$output) ) stop("arg output must be specified")


rsme <- function(model, measure) {
  sqrt( mean( (model-measure)^2 , na.rm = TRUE ) )
}

#folder = "/project/ngs/philipp/slamdunk-analysis/simulation/simulation_6/"
#version = "slamdunk/"
#tcRatePerPosition = 0.024
#readLength = 50 - 12
#sampleNumber = 21
#cfactor = 1 - dbinom(0, round(readLength / 4), tcRatePerPosition)

simulatedFileRates = opt$simulated
#simulatedFileRates = "/project/libby/slamdunk-analysis/simulation/data/test_rates_full_cov50_rl88_1/utrsummary_all_samples_rates_reads.tsv"

slamDunkFile = opt$slamdunk
#slamDunkFile = "/project/libby/slamdunk-analysis/simulation/data/test_rates_full_cov50_rl88_1/slamdunk/count/tcounts_all_samples_rates.tsv"

outputFile = opt$output
outputFileCSV = paste0(outputFile, ".tsv")

simulatedRates = read.table(simulatedFileRates, header=T, sep="\t", stringsAsFactors = F)
slamDunkRates = read.table(slamDunkFile, header=T, sep="\t", stringsAsFactors = F)

# Should not be neccessary, but for large datasets some entries are lost.
# Keep all that is found in both
inBoth = intersect(simulatedRates$Name, slamDunkRates$Name)
simulatedRates = simulatedRates[simulatedRates$Name %in% inBoth,]
slamDunkRates = slamDunkRates[slamDunkRates$Name %in% inBoth,]

fixedColumns = 11
sampleNumber = ncol(simulatedRates) - fixedColumns + 1

sampleNames = colnames(simulatedRates)[fixedColumns:(fixedColumns + sampleNumber - 1)]
simulatedSamples = simulatedRates[, fixedColumns:(fixedColumns + sampleNumber - 1)]
slamDunkSamples = slamDunkRates[, fixedColumns:(fixedColumns + sampleNumber - 1)]

pdf(outputFile)
par(mfrow=c(2,1))
boxplot(simulatedSamples - slamDunkSamples, ylim=c(-1,1), names = sampleNames, ylab="Simulated - Slamdunk", xlab="Labeled Transcripts [%]", main="", las=2)
abline(h=0, lty=2, col="grey")

boxplot(log2((simulatedSamples + 0.001) / (slamDunkSamples  + 0.001)), names = sampleNames, ylab="log2(Simulated / Slamdunk)", xlab="Labeled Transcripts [%]", main="", las=2)
abline(h=0, lty=2, col="grey")

boxplot(simulatedSamples - slamDunkSamples, ylim=c(-0.1,0.1), names = sampleNames, ylab="Simulated - Slamdunk", xlab="Labeled Transcripts [%]", main="", las=2)
abline(h=0, lty=2, col="grey")

boxplot(log2((simulatedSamples + 0.001) / (slamDunkSamples  + 0.001)), ylim=c(-1,1), names = sampleNames, ylab="log2(Simulated / Slamdunk)", xlab="Labeled Transcripts [%]", main="", las=2)
abline(h=0, lty=2, col="grey")

merged = data.frame()
#rsmeTab = data.frame(File=character(), Rate=character(), RSME=character(), stringsAsFactors=F)
rsmeTab = matrix("", ncol=3, nrow=0)
for(currentSample in 0:(sampleNumber - 1)) {
  #currentSample = 0
  current = cbind(slamDunkRates[, c(1:fixedColumns - 1, fixedColumns + currentSample)], simulatedRates[, fixedColumns + currentSample])
  colnames(current) = c(colnames(slamDunkRates[, c(1:fixedColumns - 1)]), "Simulate", "Slamdunk")
  merged = rbind(merged, current)
  
  rsmeTab = rbind(rsmeTab, c(as.character(simulatedFileRates), as.character(substring(sampleNames[currentSample + 1], 2)), as.character(rsme(current$Simulate, current$Slamdunk))))
}

par(mfrow=c(1,1))
perr = round(rsme(merged$Simulate, merged$Slamdunk), digits = 4)
pcorr = round(cor(merged$Simulate, merged$Slamdunk), digits = 4)
plot(merged$Slamdunk,  merged$Simulate, xlim=c(0,1), ylim=c(0,1), pch=4, xlab="Simulated", ylab="Slamdunk", main=paste("Cor: ", pcorr, ", RMSE: ", perr))
abline(a = 0, b = 1, col="grey", lty=2)

plot(merged$avgTcontent, merged$Slamdunk - merged$Simulate, ylim=c(-1,1), pch=4)

plot(merged$avgReadsCPM, merged$Slamdunk - merged$Simulate, ylim=c(-1,1), pch=4)

plot(merged$avgMultimapper, merged$Slamdunk - merged$Simulate, ylim=c(-1,1), pch=4)

dev.off()

rsmeTab = rbind(rsmeTab, c(as.character(simulatedFileRates), as.character(-1), as.character(rsme(merged$Simulate, merged$Slamdunk))))

write.table(rsmeTab, outputFileCSV, sep = "\t", quote = F, row.names = F, col.names = T)

