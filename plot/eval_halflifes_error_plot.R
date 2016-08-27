#!/usr/bin/env Rscript
library(getopt)

spec = matrix(c(
  'help'      , 'h', 0, "logical","print the usage of the command",
  'simulated', "s", 2,"character","Half-lifes from SlamDunk",
  'predicted', "p", 2,"character","Simulated Half-Lifes",
  'output', "o", 2,"character","Output pdf",
  'missing', "m", 2,"character","List of utrs with missing half-life"
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
if ( is.null(opt$predicted) ) stop("arg slamdunk must be specified")
if ( is.null(opt$output) ) stop("arg output must be specified")
if ( is.null(opt$missing) ) stop("arg missing must be specified")

simHLFile = opt$simulated
#simHLFile = "/project/ngs/philipp/slamdunk-analysis/simulation/simulation_2/finalAnnotation_test_cut_chrM_correct_100_original_utrs.bed"
#simHLFile = "/project/ngs/philipp/slamdunk-analysis/simulation/simulation_3/finalAnnotation_test_cut_chrM_correct_original_utrs.bed"
predHLFile = opt$predicted
#predHLFile = "/project/ngs/philipp/slamdunk-analysis/simulation/simulation_2/slamdunk/halflifes/pooja_UTR_annotation_examples_sample_1_0min_reads_slamdunk_mapped_filtered_tcount_halflifes.tsv"
#predHLFile = "/project/ngs/philipp/slamdunk-analysis/simulation/simulation_3/slamdunk/halflifes/finalAnnotation_test_1_0min_reads_slamdunk_mapped_filtered_tcount_halflifes.tsv"
output = opt$output
missing = opt$missing

simHL = read.csv(simHLFile, sep="\t", header = F)
colnames(simHL) = c("chr", "start", "stop", "name", "halflife", "strand")
predHL = read.csv(predHLFile, sep="\t")

predHL$simulated_hl = simHL$halflife

tmp = predHL[is.na(predHL$score), ]
predHL = predHL[!is.na(predHL$score), ]

predHL$log2Diff = log2((predHL$score + 0.1) / (predHL$simulated_hl + 0.1))
rmse = sqrt(sum((predHL$simulated_hl - predHL$score) ^ 2) / nrow(predHL))
avgHL = mean(predHL$simulated_hl)

predHLUniq = predHL[predHL$multiMapCount == 0,]
predHLMulti = predHL[predHL$multiMapCount > 0,]

pdf(output, height = 6, width = 9)
plot(predHLUniq$readsCPM, predHLUniq$log2Diff, main=paste0("Avg. HL: ", round(avgHL), "\nRMSE: ", round(rmse)), pch=4, ylim=c(-8,8), xlab="couts per million (reads)", ylab="log2 FC (half-lifes)")
points(predHLMulti$readsCPM, predHLMulti$log2Diff, pch=4, col="red")
abline(h=0, lty=2, col="grey")

plot(predHLUniq$simulated_hl / 60, predHLUniq$log2Diff, main=paste0("Avg. HL: ", round(avgHL), "\nRMSE: ", round(rmse)), pch=4, ylim=c(-8,8), xlab="Simulated half-life (hours)", ylab="log2 FC (half-lifes)")
points(predHLMulti$simulated_hl / 60, predHLMulti$log2Diff, pch=4, col="red")
abline(h=0, lty=2, col="grey")
dev.off()

write.table(tmp, missing, quote = F, row.names = F, col.names = T)
