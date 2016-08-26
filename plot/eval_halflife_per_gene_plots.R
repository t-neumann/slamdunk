#!/usr/bin/env Rscript
library(getopt)

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
if ( is.null(opt$conversionrate) ) { opt$conversionrate = 0.2 }

simulatedFiles = opt$simulated
#simulatedFiles = "/project/ngs/philipp/slamdunk-analysis/simulation/simulation_1/pooja_UTR_annotation_examples_1_0min_utrsummary.csv,/project/ngs/philipp/slamdunk-analysis/simulation/simulation_1/pooja_UTR_annotation_examples_2_15min_utrsummary.csv,/project/ngs/philipp/slamdunk-analysis/simulation/simulation_1/pooja_UTR_annotation_examples_3_30min_utrsummary.csv,/project/ngs/philipp/slamdunk-analysis/simulation/simulation_1/pooja_UTR_annotation_examples_4_60min_utrsummary.csv,/project/ngs/philipp/slamdunk-analysis/simulation/simulation_1/pooja_UTR_annotation_examples_5_180min_utrsummary.csv,/project/ngs/philipp/slamdunk-analysis/simulation/simulation_1/pooja_UTR_annotation_examples_6_360min_utrsummary.csv,/project/ngs/philipp/slamdunk-analysis/simulation/simulation_1/pooja_UTR_annotation_examples_7_720min_utrsummary.csv,/project/ngs/philipp/slamdunk-analysis/simulation/simulation_1/pooja_UTR_annotation_examples_8_1440min_utrsummary.csv"
slamDunkFiles = opt$slamdunk
#slamDunkFiles = "/project/ngs/philipp/slamdunk-analysis/simulation/simulation_1/slamdunk/count/pooja_UTR_annotation_examples_1_0min_reads_slamdunk_mapped_filtered_tcount.csv,/project/ngs/philipp/slamdunk-analysis/simulation/simulation_1/slamdunk/count/pooja_UTR_annotation_examples_2_15min_reads_slamdunk_mapped_filtered_tcount.csv,/project/ngs/philipp/slamdunk-analysis/simulation/simulation_1/slamdunk/count/pooja_UTR_annotation_examples_3_30min_reads_slamdunk_mapped_filtered_tcount.csv,/project/ngs/philipp/slamdunk-analysis/simulation/simulation_1/slamdunk/count/pooja_UTR_annotation_examples_4_60min_reads_slamdunk_mapped_filtered_tcount.csv,/project/ngs/philipp/slamdunk-analysis/simulation/simulation_1/slamdunk/count/pooja_UTR_annotation_examples_5_180min_reads_slamdunk_mapped_filtered_tcount.csv,/project/ngs/philipp/slamdunk-analysis/simulation/simulation_1/slamdunk/count/pooja_UTR_annotation_examples_6_360min_reads_slamdunk_mapped_filtered_tcount.csv,/project/ngs/philipp/slamdunk-analysis/simulation/simulation_1/slamdunk/count/pooja_UTR_annotation_examples_7_720min_reads_slamdunk_mapped_filtered_tcount.csv,/project/ngs/philipp/slamdunk-analysis/simulation/simulation_1/slamdunk/count/pooja_UTR_annotation_examples_8_1440min_reads_slamdunk_mapped_filtered_tcount.csv"
filesSimulated = as.character(ordered(strsplit(simulatedFiles, ",")[[1]]))
filesSlamDunk = as.character(ordered(strsplit(slamDunkFiles, ",")[[1]]))
outputFile = opt$output
#outputFile = "/project/ngs/philipp/slamdunk-analysis/simulation/simulation_1/eval/halflife_per_gene_eval_plots.pdf"
#timesParameter = "0,15,30,60,180,360,720,1440"
timesParameter = opt$timepoints
times = as.numeric(strsplit(timesParameter, ",")[[1]])
times = times / 60
#bedFile = "/project/ngs/philipp/slamdunk-analysis/simulation/simulation_1/pooja_UTR_annotation_examples_original_utrs.bed"
bedFile = opt$bed
#conversionRate = 0.03
conversionRate = opt$conversionrate


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

pdf(outputFile, height=6, width=9)
for(utr in 1:nrow(slamDunkMergedRates)) {
  #utr = 8
  pulseSlamDunk = data.frame(y = as.numeric(t(slamDunkMergedRates[utr, 6:(5 + length(times))])[,1]), x = times)
  pulseSimulated = data.frame(y = as.numeric(t(simMergedRates[utr, 6:(5 + length(times))])[,1]), x = times)
  #yLim = max(max(pulseSlamDunk$y), max(pulseSimulated$y))
  yLim = conversionRate
  yLab = "conversion rate"
  if(perRead) {
    yLab = "% of T->C reads"
  }
  
  # Infere half life from data
  y = pulseSlamDunk$y
  timepoints = pulseSimulated$x
  
  a_start<-max(y) #param a is the y value when x=0
  k_start = log(2, base = exp(1))/5
  
  fit = nls(y ~ a*(1-exp(-k*(timepoints))), start=list(a=a_start,k=k_start))
  halfLifePred = log(2, base = exp(1))/coef(fit)[2]
  
  halfLifeSim = bed[utr, ]$score / 60
  plot(0, type="n", main=paste0(slamDunkMergedRates[utr, ]$name, "\n half life: ", round(halfLifeSim, digits = 2), " (sim), ", round(halfLifePred, digits = 2)," (pred)"), xlab="Time (hours)", ylab=yLab, ylim=c(0, yLim), xlim=c(times[1], times[length(times)]), pch=4)
  lines(pulseSlamDunk$x, pulseSlamDunk$y, type="p", col="blue", lty=1, pch=4)     
  lines(pulseSimulated$x, pulseSimulated$y, type="p", col="green", lty=1, pch=4)     
  legend("bottomright", c("rates (slamDunk)", "rates (simulated)", "theoretical (slamDunk)", "theoretical (simulated)"), col=c("blue", "green", "red", "grey"), lty=c(1, 1, 2, 2), bty="n")
  
  lambda = log(2) / (bed[utr, ]$score / 60)
  t = 0:max(times)
  lines((1 - exp(-lambda*t)) * conversionRate ~ t, type="l", lty=2, col="grey")
  
  #lines((1 - exp(-coef(fit)[2]*timepoints)) * coef(fit)[1] ~ timepoints, type="l", lty=2, col="red")
  lines((1 - exp(-coef(fit)[2]*t)) * coef(fit)[1] ~ t, type="l", lty=2, col="red")
  
}  
dev.off()

# y = pulseSlamDunk$y
# timepoints = pulseSimulated$x
# 
# a_start<-max(y) #param a is the y value when x=0
# k_start<-2*log(2)/a_start #b is the decay rate
# k_start = log(2, base = exp(1))/5
# k_start = lambda
# 
# fit = nls(y ~ a*(1-exp(-k*(timepoints))), start=list(a=a_start,k=k_start))
# lines(timepoints,predict(fit), col="red",lty=2)
# 
# lines((1 - exp(-coef(fit)[2]*t)) * coef(fit)[1] ~ t, type="l", lty=2, col="black")
# halfLives_pulse = log(2, base = exp(1))/coef(fit)[2]
# 
# 
