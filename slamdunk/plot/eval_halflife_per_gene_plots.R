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
  'slamdunk', "f", 2,"character","Comma seperated list of SlamDunk results",
  'timepoints', "t", 2,"character","Comma seperated list of time points",
  'bed', "b", 2,"character","BED file containing half lifes",
  'output', "o", 2,"character","Output pdf",
  'conversionrate', "c", 2,"character","Simulated conversion rate"
),ncol = 5,byrow=T)

opt = getopt(spec)

if ( !is.null(opt$help) || length(opt)==4 ) {
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
if ( is.null(opt$timepoints) ) stop("arg timepoints must be specified")
if ( is.null(opt$bed) ) stop("arg bed specified")
if ( is.null(opt$conversionrate) ) { opt$conversionrate = 0.03 }

simulatedFiles = opt$simulated
#simulatedFiles = "/project/ngs/philipp/slamdunk-analysis/simulation/simulation_1/pooja_UTR_annotation_examples_1_0min_utrsummary.csv,/project/ngs/philipp/slamdunk-analysis/simulation/simulation_1/pooja_UTR_annotation_examples_2_15min_utrsummary.csv,/project/ngs/philipp/slamdunk-analysis/simulation/simulation_1/pooja_UTR_annotation_examples_3_30min_utrsummary.csv,/project/ngs/philipp/slamdunk-analysis/simulation/simulation_1/pooja_UTR_annotation_examples_4_60min_utrsummary.csv,/project/ngs/philipp/slamdunk-analysis/simulation/simulation_1/pooja_UTR_annotation_examples_5_180min_utrsummary.csv,/project/ngs/philipp/slamdunk-analysis/simulation/simulation_1/pooja_UTR_annotation_examples_6_360min_utrsummary.csv,/project/ngs/philipp/slamdunk-analysis/simulation/simulation_1/pooja_UTR_annotation_examples_7_720min_utrsummary.csv,/project/ngs/philipp/slamdunk-analysis/simulation/simulation_1/pooja_UTR_annotation_examples_8_1440min_utrsummary.csv"
#simulatedFiles = "/project/ngs/philipp/slamdunk-analysis/simulation/simulation_2/pooja_UTR_annotation_examples_sample_1_0min_utrsummary.csv,/project/ngs/philipp/slamdunk-analysis/simulation/simulation_2/pooja_UTR_annotation_examples_sample_2_15min_utrsummary.csv,/project/ngs/philipp/slamdunk-analysis/simulation/simulation_2/pooja_UTR_annotation_examples_sample_3_30min_utrsummary.csv,/project/ngs/philipp/slamdunk-analysis/simulation/simulation_2/pooja_UTR_annotation_examples_sample_4_60min_utrsummary.csv,/project/ngs/philipp/slamdunk-analysis/simulation/simulation_2/pooja_UTR_annotation_examples_sample_5_180min_utrsummary.csv,/project/ngs/philipp/slamdunk-analysis/simulation/simulation_2/pooja_UTR_annotation_examples_sample_6_360min_utrsummary.csv,/project/ngs/philipp/slamdunk-analysis/simulation/simulation_2/pooja_UTR_annotation_examples_sample_7_720min_utrsummary.csv,/project/ngs/philipp/slamdunk-analysis/simulation/simulation_2/pooja_UTR_annotation_examples_sample_8_1440min_utrsummary.csv"

slamDunkFiles = opt$slamdunk
#slamDunkFiles = "/project/ngs/philipp/slamdunk-analysis/simulation/simulation_1/slamdunk/count/pooja_UTR_annotation_examples_1_0min_reads_slamdunk_mapped_filtered_tcount.csv,/project/ngs/philipp/slamdunk-analysis/simulation/simulation_1/slamdunk/count/pooja_UTR_annotation_examples_2_15min_reads_slamdunk_mapped_filtered_tcount.csv,/project/ngs/philipp/slamdunk-analysis/simulation/simulation_1/slamdunk/count/pooja_UTR_annotation_examples_3_30min_reads_slamdunk_mapped_filtered_tcount.csv,/project/ngs/philipp/slamdunk-analysis/simulation/simulation_1/slamdunk/count/pooja_UTR_annotation_examples_4_60min_reads_slamdunk_mapped_filtered_tcount.csv,/project/ngs/philipp/slamdunk-analysis/simulation/simulation_1/slamdunk/count/pooja_UTR_annotation_examples_5_180min_reads_slamdunk_mapped_filtered_tcount.csv,/project/ngs/philipp/slamdunk-analysis/simulation/simulation_1/slamdunk/count/pooja_UTR_annotation_examples_6_360min_reads_slamdunk_mapped_filtered_tcount.csv,/project/ngs/philipp/slamdunk-analysis/simulation/simulation_1/slamdunk/count/pooja_UTR_annotation_examples_7_720min_reads_slamdunk_mapped_filtered_tcount.csv,/project/ngs/philipp/slamdunk-analysis/simulation/simulation_1/slamdunk/count/pooja_UTR_annotation_examples_8_1440min_reads_slamdunk_mapped_filtered_tcount.csv"
#slamDunkFiles = "/project/ngs/philipp/slamdunk-analysis/simulation/simulation_2/slamdunk/count/pooja_UTR_annotation_examples_sample_1_0min_reads_slamdunk_mapped_filtered_tcount.csv,/project/ngs/philipp/slamdunk-analysis/simulation/simulation_2/slamdunk/count/pooja_UTR_annotation_examples_sample_2_15min_reads_slamdunk_mapped_filtered_tcount.csv,/project/ngs/philipp/slamdunk-analysis/simulation/simulation_2/slamdunk/count/pooja_UTR_annotation_examples_sample_3_30min_reads_slamdunk_mapped_filtered_tcount.csv,/project/ngs/philipp/slamdunk-analysis/simulation/simulation_2/slamdunk/count/pooja_UTR_annotation_examples_sample_4_60min_reads_slamdunk_mapped_filtered_tcount.csv,/project/ngs/philipp/slamdunk-analysis/simulation/simulation_2/slamdunk/count/pooja_UTR_annotation_examples_sample_5_180min_reads_slamdunk_mapped_filtered_tcount.csv,/project/ngs/philipp/slamdunk-analysis/simulation/simulation_2/slamdunk/count/pooja_UTR_annotation_examples_sample_6_360min_reads_slamdunk_mapped_filtered_tcount.csv,/project/ngs/philipp/slamdunk-analysis/simulation/simulation_2/slamdunk/count/pooja_UTR_annotation_examples_sample_7_720min_reads_slamdunk_mapped_filtered_tcount.csv,/project/ngs/philipp/slamdunk-analysis/simulation/simulation_2/slamdunk/count/pooja_UTR_annotation_examples_sample_8_1440min_reads_slamdunk_mapped_filtered_tcount.csv"

filesSimulated = as.character(ordered(strsplit(simulatedFiles, ",")[[1]]))
filesSlamDunk = as.character(ordered(strsplit(slamDunkFiles, ",")[[1]]))
outputFile = opt$output
#outputFile = "/project/ngs/philipp/slamdunk-analysis/simulation/simulation_1/eval/halflife_per_gene_eval_plots.pdf"
#timesParameter = "0,15,30,60,180,360,720,1440"
timesParameter = opt$timepoints
times = as.numeric(strsplit(timesParameter, ",")[[1]])
times = times / 60
#bedFile = "/project/ngs/philipp/slamdunk-analysis/simulation/simulation_1/pooja_UTR_annotation_examples_original_utrs.bed"
#bedFile = "/project/ngs/philipp/slamdunk-analysis/simulation/simulation_2/finalAnnotation_test_cut_chrM_correct_100_original_utrs.bed"
bedFile = opt$bed
#conversionRate = 0.03
conversionRate = as.numeric(opt$conversionrate)

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
  
  c(halfLifePred, C, k)
}

mergeRates <- function(times, files, perRead) {
  mergedRates = data.frame()
  for(i in 1:length(times)) {
    time = times[i]
    #print(time)
    simDataFile = files[i]
    simulation = read.table(simDataFile)
    colnames(simulation) = c("chr", "start", "stop", "name", "strand", "conversionRate", "readsCPM", "tCount", "tcCount", "readCount", "convertedReads", "multiMapCount")
    if(nrow(mergedRates) == 0) {
      mergedRates = simulation[, c("chr", "start", "stop", "name", "strand")]
      if(perRead == TRUE) {
        mergedRates$conversionRate = simulation$convertedReads / simulation$readCount
      } else {
        mergedRates$conversionRate = simulation$conversionRate
      }
    } else {
      if(perRead == TRUE) {
        mergedRates = cbind(mergedRates, simulation$convertedReads / simulation$readCount)
      } else {
        mergedRates = cbind(mergedRates, simulation$conversionRate)
      }
    }
  }
  colnames(mergedRates) = c("chr", "start", "stop", "name", "strand", times)
  mergedRates
}

bed = read.table(bedFile)
colnames(bed) = c("char", "start", "stop", "name", "score", "strand")



perRead = F
slamDunkMergedRates = mergeRates(times, filesSlamDunk, perRead)
simMergedRates = mergeRates(times, filesSimulated, perRead)

pageNumber = 1
#outputFile = "~/test.pdf"
#pdf(outputFile, height=6, width=9)
pdf(paste0(outputFile, "_page_", pageNumber, ".pdf"), height=6, width=9)
for(utr in 1:nrow(slamDunkMergedRates)) {
  #utr = 1
  pulseSlamDunk = data.frame(y = as.numeric(t(slamDunkMergedRates[utr, 6:(5 + length(times))])[,1]), x = times)
  pulseSimulated = data.frame(y = as.numeric(t(simMergedRates[utr, 6:(5 + length(times))])[,1]), x = times)
  #yLim = max(max(pulseSlamDunk$y), max(pulseSimulated$y))
  yLim = conversionRate * 1.25
  yLab = "conversion rate"
  if(perRead) {
    yLab = "% of T->C reads"
  }
  
  # Infere half life from data
  halfLifeResultSlamDunk = computeHalfLife(pulseSlamDunk$y, pulseSlamDunk$x)
  halfLifePred = halfLifeResultSlamDunk[1]
  halfLifeResultSimulated = computeHalfLife(pulseSimulated$y, pulseSimulated$x)
  halfLifeSim = halfLifeResultSimulated[1]
  halfLifeTruth = bed[utr, ]$score
  
  
  plot(0, type="n", main=paste0(slamDunkMergedRates[utr, ]$name, "\n half life: ", round(halfLifeTruth, digits = 0), " (truth), ", round(halfLifeSim, digits = 0), " (sim), ", round(halfLifePred, digits = 0)," (slamDunk)"), xlab="Time (hours)", ylab=yLab, ylim=c(0, yLim), xlim=c(times[1], times[length(times)]), pch=4)
  lines(pulseSimulated$x, pulseSimulated$y, type="p", col="green", lty=1, pch=4)     
  lines(pulseSlamDunk$x, pulseSlamDunk$y, type="p", col="blue", lty=1, pch=4)     
  legend("bottomright", c("rates (slamDunk)", "rates (simulated)", "slamDunk", "simulated", "truth"), col=c("blue", "green", "blue", "green", "grey"), lty=c(1, 1, 2, 2, 2), bty="n")
  
  
  t = 0:max(times)
  # Print truth
  lambda = log(2) / (halfLifeTruth / 60)
  lines((1 - exp(-lambda*t)) * conversionRate ~ t, type="l", lty=2, col="grey")
  # Print simulated
  lines((1 - exp(-halfLifeResultSimulated[3]*t)) * halfLifeResultSimulated[2] ~ t, type="l", lty=2, col="green")
  # Print slamDunk
  lines((1 - exp(-halfLifeResultSlamDunk[3]*t)) * halfLifeResultSlamDunk[2] ~ t, type="l", lty=2, col="blue")
  
  if(utr %% 100 == 0) {
    dev.off()
    pageNumber = pageNumber + 1
    pdf(paste0(outputFile, "_page_", pageNumber, ".pdf"), height=6, width=9)
    print(paste0(outputFile, "_page_", pageNumber, ".pdf"))
  }
  
}  
dev.off()
