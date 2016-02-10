#!/software/R/R-3.1.2/bin/Rscript
#!/usr/bin/Rscript
#!/software/R/R-3.1.2/bin/Rscript
#!/biosw/debian7-x86_64/R/3.0.2/bin/Rscript

library(getopt)


spec = matrix(c(
				'help'      , 'h', 0, "logical","print the usage of the command",
				'inputFile', "i", 2,"character","tsv table of mutations per position",
				'outputFile', "o", 2,"character","output pdf file name"
		),ncol = 5,byrow=T)

opt = getopt(spec)

if ( !is.null(opt$help) || length(opt)==1 ) {
	#get the script name
	cmd = commandArgs(FALSE)
	self = strsplit(cmd[grep("--file",cmd)],"=")[[1]][2]
	cat(basename(self),": Create mismatches per read position plots.\n\n")
	#print a friendly message and exit with a non-zero error code
	cat(getopt(spec,command = self,usage=T))
	q(status=1);
}


if ( is.null(opt$inputFile) ) stop("arg input must be specified")
if ( is.null(opt$outputFile) ) { opt$outputFile = paste(opt$inputFile, ".pdf", sep="") }


mut = read.table(opt$inputFile)


#mut = read.table("test_mut_bowtie.csv")

#totalFwd = mut[1,1]
#totalRev = mut[1,2]
#tcFwd = mut[1,3]
#tcRev = mut[1,4]

#mut = mut[-1,]

counts = rbind(c(mut$V1)/c(mut$V5) * 100, c(mut$V2)/c(mut$V6) * 100)
countsTC = rbind(c(mut$V3)/c(mut$V5) * 100, c(mut$V4)/c(mut$V6) * 100)

pdf(opt$outputFile, width=10, height=10)
par(mfrow=c(2,1))

barplot(counts, beside=T, names.arg=1:nrow(mut), main="All mutations", ylim=c(0,10), xlab="Position on read", ylab="% of reads with mutation", legend=c("forward", "reverse"))
barplot(countsTC, beside=T, names.arg=1:nrow(mut), main="T->C on fwd, A->G on rev", ylim=c(0,1), xlab="Position on read", ylab="% of reads with mutation", legend=c("forward", "reverse"))

dev.off()