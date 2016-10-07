#!/usr/bin/env Rscript

# Load packages only from local Rslamdunk library 
libLoc = .libPaths()[grep("Rslamdunk",.libPaths())]

# Check if libraries are available, install otherwise
source(paste(libLoc,'/../checkLibraries.R',sep=""))

checkLib(libLoc)

library(getopt, lib.loc = libLoc)

spec = matrix(c(
  'help'      , 'h', 0, "logical","print the usage of the command",
  'simulated', "s", 2,"character","Half-lifes inferred from SlamDunk results",
  'predicted', "p", 2,"character","Simulated Half-Lifes",
  'truth', "t", 2,"character","True Half-Lifes",
  'output', "o", 2,"character","Output pdf",
  'missing', "m", 2,"character","List of utrs with missing half-life"
),ncol = 5,byrow=T)

opt = getopt(spec)

if ( !is.null(opt$help) || length(opt)==5 ) {
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
if ( is.null(opt$truth) ) stop("arg truth must be specified")
if ( is.null(opt$output) ) stop("arg output must be specified")
if ( is.null(opt$missing) ) stop("arg missing must be specified")

truthFile = opt$truth
#truthFile = "/project/ngs/philipp/slamdunk-analysis/simulation/simulation_2/finalAnnotation_test_cut_chrM_correct_100_original_utrs.bed"
#truthFile = "/project/ngs/philipp/slamdunk-analysis/simulation/simulation_3/finalAnnotation_test_cut_chrM_correct_original_utrs.bed"
simHLFile = opt$simulated
#simHLFile = "/project/ngs/philipp/slamdunk-analysis/simulation/simulation_2/pooja_UTR_annotation_examples_sample_1_0min_utrsummary_halflifes.tsv"
#simHLFile = "/project/ngs/philipp/slamdunk-analysis/simulation/simulation_3/finalAnnotation_test_1_0min_utrsummary_halflifes.tsv"
predHLFile = opt$predicted
#predHLFile = "/project/ngs/philipp/slamdunk-analysis/simulation/simulation_2/slamdunk/halflifes/pooja_UTR_annotation_examples_sample_1_0min_reads_slamdunk_mapped_filtered_tcount_halflifes.tsv"
#predHLFile = "/project/ngs/philipp/slamdunk-analysis/simulation/simulation_3/slamdunk/halflifes/finalAnnotation_test_1_0min_reads_slamdunk_mapped_filtered_tcount_halflifes.tsv"
output = opt$output
missing = opt$missing

trueHL = read.csv(truthFile, sep="\t", header = F)
colnames(trueHL) = c("chr", "start", "stop", "name", "halflife", "strand")
simHL = read.csv(simHLFile, sep="\t", header = T)
predHL = read.csv(predHLFile, sep="\t")

predHL$simulated_hl = simHL$score
predHL$true_hl = trueHL$halflife
predHL$multiperc = predHL$multiMapCount / predHL$readsCPM
head(simHL)

tmp = predHL[is.na(predHL$score) | is.na(predHL$simulated_hl), ]
predHL = predHL[!is.na(predHL$score) & !is.na(predHL$simulated_hl), ]

predHL$log2DiffSim = log2((predHL$score + 0.1) / (predHL$simulated_hl + 0.1))
predHL$log2DiffTrue = log2((predHL$score + 0.1) / (predHL$true_hl + 0.1))
rmseSim = sqrt(sum((predHL$simulated_hl - predHL$score) ^ 2) / nrow(predHL))
rmseTrue = sqrt(sum((predHL$true_hl - predHL$score) ^ 2) / nrow(predHL))
avgHL = mean(predHL$true_hl)

predHLUniq = predHL[predHL$multiMapCount == 0,]
predHLMulti = predHL[predHL$multiMapCount > 0,]

#ggplot(predHLUniq,aes(x=readsCPM,y=log2Diff)) + stat_binhex()
head(predHLUniq)

pdf(output, height = 6, width = 9)
plot(0, main=paste0("Unique UTRs (", nrow(predHLUniq), ")\nAvg. HL: ", round(avgHL), ", RMSE: ", round(rmseSim)), pch=4, xlim=c(min(predHLUniq$readsCPM), max(predHLUniq$readsCPM)), ylim=c(-6,6), xlab="couts per million (reads)", ylab="log2 FC (half-lifes)", type="n")
points(predHLUniq$readsCPM, predHLUniq$log2DiffSim, pch=1, col="#00000033")
abline(h=0, lty=2, col="grey")
plot(0, main=paste0("Multimapper UTRs (", nrow(predHLMulti) ,")\nAvg. HL: ", round(avgHL), ", RMSE: ", round(rmseSim)), pch=4, xlim=c(min(predHLUniq$readsCPM), max(predHLUniq$readsCPM)), ylim=c(-6,6), xlab="couts per million (reads)", ylab="log2 FC (half-lifes)", type="n")
points(predHLMulti$readsCPM, predHLMulti$log2DiffSim, pch=1, col="#00000033")
abline(h=0, lty=2, col="grey")

plot(0, main=paste0("Unique UTRs (", nrow(predHLUniq), ")\nAvg. HL: ", round(avgHL), ", RMSE: ", round(rmseSim)), pch=4, xlim=c(min(predHLUniq$true_hl), max(predHLUniq$true_hl)), ylim=c(-6,6), xlab="couts per million (reads)", ylab="log2 FC (half-lifes)", type="n")
points(predHLUniq$true_hl, predHLUniq$log2DiffSim, pch=1, col="#00000033")
abline(h=0, lty=2, col="grey")
plot(0, main=paste0("Multimapper UTRs (", nrow(predHLMulti) ,")\nAvg. HL: ", round(avgHL), ", RMSE: ", round(rmseSim)), pch=4, xlim=c(min(predHLUniq$true_hl), max(predHLUniq$true_hl)), ylim=c(-6,6), xlab="couts per million (reads)", ylab="log2 FC (half-lifes)", type="n")
points(predHLMulti$true_hl, predHLMulti$log2DiffSim, pch=1, col="#00000033")
abline(h=0, lty=2, col="grey")


#plot(0, main=paste0("Multimapper UTRs (", nrow(predHLMulti) ,")\nAvg. HL: ", round(avgHL), ", RMSE: ", round(rmseSim)), pch=4, xlim=c(min(predHLUniq$multiperc), max(predHLUniq$multiperc)), ylim=c(-6,6), xlab="couts per million (reads)", ylab="log2 FC (half-lifes)", type="n")
#points(predHLMulti$multiperc, predHLMulti$log2DiffSim, pch=1, col="#00000033")
#abline(h=0, lty=2, col="grey")


#lim = max(predHLUniq$simulated_hl) * 1.25
lim = 1440
corr = cor(predHLUniq$simulated_hl, predHLUniq$score)
plot(predHLUniq$simulated_hl ~ predHLUniq$score, main=paste0("Simulated vs. SlamDunk (unique)\nPearson: ", round(corr, digits = 4)), xlim=c(0, lim), ylim=c(0, lim), ylab="Half-life (simulated)", xlab="Half-Life (slamDunk)", col="#00000033")
abline(a = 0, b = 1, col="grey", lty=2)

corr = cor(predHLMulti$simulated_hl, predHLMulti$score)
plot(predHLMulti$simulated_hl ~ predHLMulti$score, main=paste0("Simulated vs. SlamDunk (multimapper)\nPearson: ", round(corr, digits = 4)), xlim=c(0, lim), ylim=c(0, lim), ylab="Half-life (simulated)", xlab="Half-Life (slamDunk)", col="#00000033")
abline(a = 0, b = 1, col="grey", lty=2)


corr = cor(predHLUniq$true_hl, predHLUniq$score)
plot(predHLUniq$true_hl ~ predHLUniq$score, main=paste0("True vs. SlamDunk (unique)\nPearson: ", round(corr, digits = 4)), xlim=c(0, lim), ylim=c(0, lim), ylab="Half-life (true)", xlab="Half-Life (slamDunk)", col="#00000033")
abline(a = 0, b = 1, col="grey", lty=2)

corr = cor(predHLMulti$true_hl, predHLMulti$score)
plot(predHLMulti$true_hl ~ predHLMulti$score, main=paste0("True vs. SlamDunk (multimapper)\nPearson: ", round(corr, digits = 4)), xlim=c(0, lim), ylim=c(0, lim), ylab="Half-life (true)", xlab="Half-Life (slamDunk)", col="#00000033")
abline(a = 0, b = 1, col="grey", lty=2)


dev.off()

#ggplot(predHLUniq,aes(x=simulated_hl,y=log2Diff)) + geom_point(alpha = 0.3)
#p = p + stat_binhex(bins=100)
#p = p 
#ggplot(predHLMulti,aes(x=simulated_hl,y=log2Diff)) + stat_binhex(bins=100)



write.table(tmp, missing, quote = F, row.names = F, col.names = T)
