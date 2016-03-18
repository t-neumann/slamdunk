#!/usr/bin/env python

#########################################################################
# Main routine for the SLAMdunk analyzer
#########################################################################
# Imports
#########################################################################

from __future__ import print_function
import sys, os

from argparse import ArgumentParser, RawDescriptionHelpFormatter
    
from os.path import basename

from joblib import Parallel, delayed
from dunks import tcounter, mapper, filter, deduplicator, stats, snps, dump
from utils.misc import replaceExtension, readSampleNames

########################################################################
# Global variables
########################################################################

printOnly = False
verbose = False

mainOutput = sys.stderr

logToMainOutput = False

########################################################################
# Routine definitions
########################################################################

def getLogFile(path):
    if(logToMainOutput):
        return mainOutput
    else:
        log = open(path, "a")
        return log
    
def closeLogFile(log):
    if(not logToMainOutput):
        log.close()

def message(msg):
    print(msg, file=mainOutput)

def error(msg, code=-1):
    print(msg, file=mainOutput)
    sys.exit(code)
    
def stepFinished():
    print(".", end="", file=mainOutput)

def dunkFinished():
    print("", file=mainOutput)
    
def createDir(directory):
    if not os.path.exists(directory):
        message("Creating output directory: " + directory)
        os.makedirs(directory)

def runMap(tid, inputBAM, referneceFile, threads, trim5p, outputDirectory) :
    createDir(outputDirectory)
    outputSAM = os.path.join(outputDirectory, replaceExtension(basename(inputBAM), ".sam", "_slamdunk_mapped"))
    outputLOG = os.path.join(outputDirectory, replaceExtension(basename(inputBAM), ".log", "_slamdunk_mapped"))
    mapper.Map(inputBAM, referneceFile, outputSAM, getLogFile(outputLOG), threads=threads, trim5p=trim5p, printOnly=printOnly, verbose=verbose)
    stepFinished()

def runSort(tid, bam, outputDirectory):
    createDir(outputDirectory)
    inputSAM = os.path.join(outputDirectory, replaceExtension(basename(bam), ".sam", "_slamdunk_mapped"))
    outputBAM = os.path.join(outputDirectory, replaceExtension(basename(bam), ".bam", "_slamdunk_mapped"))
    outputLOG = os.path.join(outputDirectory, replaceExtension(basename(bam), ".log", "_slamdunk_mapped"))
    mapper.sort(inputSAM, outputBAM, getLogFile(outputLOG), False, printOnly, verbose)
    stepFinished()

def runDedup(tid, bam, outputDirectory) :
    createDir(outputDirectory)
    outputBAM = os.path.join(outputDirectory, replaceExtension(basename(bam), ".bam", "_dedup"))
    outputLOG = os.path.join(outputDirectory, replaceExtension(basename(bam), ".log", "_dedup"))
    log = getLogFile(outputLOG)
    deduplicator.Dedup(bam, outputBAM, log)
    closeLogFile(log)
    stepFinished()
        
def runFilter(tid, bam, outputDirectory):
    createDir(outputDirectory)
    outputBAM = os.path.join(outputDirectory, replaceExtension(basename(bam), ".bam", "_filtered"))
    outputLOG = os.path.join(outputDirectory, replaceExtension(basename(bam), ".log", "_filtered"))
    filter.Filter(bam, outputBAM, getLogFile(outputLOG), args.mq, 0.85, -1, printOnly, verbose)
    stepFinished()

def runSnp(tid, referenceFile, minCov, minVarFreq, inputBAM, outputDirectory) :
    createDir(outputDirectory)
    outputSNP = os.path.join(outputDirectory, replaceExtension(basename(inputBAM), ".vcf", "_snp"))
    outputLOG = os.path.join(outputDirectory, replaceExtension(basename(inputBAM), ".log", "_snp"))
    snps.SNPs(inputBAM, outputSNP, referenceFile, minVarFreq, minCov, getLogFile(outputLOG), printOnly, verbose, True)
    stepFinished()
                
def runCount(tid, bam, outputDirectory, snpDirectory) :
    createDir(outputDirectory)
    outputCSV = os.path.join(outputDirectory, replaceExtension(basename(bam), ".csv", "_tcount"))
    outputBedgraphPlus = os.path.join(outputDirectory, replaceExtension(basename(bam), ".bedgraph", "_tcount_plus"))
    outputBedgraphMinus = os.path.join(outputDirectory, replaceExtension(basename(bam), ".bedgraph", "_tcount_mins"))
    outputLOG = os.path.join(outputDirectory, replaceExtension(basename(bam), ".log", "_tcount"))
    if(snpDirectory != None):
        inputSNP = os.path.join(snpDirectory, replaceExtension(basename(bam), ".vcf", "_snp"))
    else:
        inputSNP = None
    tcounter.computeTconversions(args.ref, args.bed, inputSNP, bam, args.maxLength, args.minQual, outputCSV, outputBedgraphPlus, outputBedgraphMinus, getLogFile(outputLOG))
    stepFinished()
    return outputCSV

def runCountCombine(bams, sampleNames, outputPrefix, outputDirectory):
    createDir(outputDirectory)
    samples = readSampleNames(sampleNames, bams)
    
    NON_TC_READ_COUNT = 4
    NON_TC_NORM_READ_COUNT = 5
    TC_READ_COUNT = 6
    TC_READ_NORM_COUNT = 7
    TC_READ_PERC = 8
    outputFile = os.path.join(outputDirectory, outputPrefix + "_non_tc_counts.csv")
    tcounter.summary(bams, samples, outputFile, NON_TC_READ_COUNT)
    outputFile = os.path.join(outputDirectory, outputPrefix + "_non_tc_norm_counts.csv")
    tcounter.summary(bams, samples, outputFile, NON_TC_NORM_READ_COUNT)
    outputFile = os.path.join(outputDirectory, outputPrefix + "_tc_counts.csv")
    tcounter.summary(bams, samples, outputFile, TC_READ_COUNT)
    outputFile = os.path.join(outputDirectory, outputPrefix + "_tc_norm_counts.csv")
    tcounter.summary(bams, samples, outputFile, TC_READ_NORM_COUNT)
    outputFile = os.path.join(outputDirectory, outputPrefix + "_tc_percentage.csv")
    tcounter.summary(bams, samples, outputFile, TC_READ_PERC)
    stepFinished()
        
        
def runStatsRates(tid, bam, referenceFile, minMQ, outputDirectory) :
    createDir(outputDirectory)
    outputCSV = os.path.join(outputDirectory, replaceExtension(basename(bam), ".csv", "_tcount_overallrates"))
    outputPDF = os.path.join(outputDirectory, replaceExtension(basename(bam), ".pdf", "_tcount_overallrates"))
    outputLOG = os.path.join(outputDirectory, replaceExtension(basename(bam), ".log", "_tcount_overallrates"))
    log = getLogFile(outputLOG)
    stats.statsComputeOverallRates(referenceFile, bam, minMQ, outputCSV, outputPDF, log)
    closeLogFile(log)
    stepFinished()
    
def runSTcPerReadPos(tid, bam, referenceFile, minMQ, maxReadLength, outputDirectory, snpDirectory):
    createDir(outputDirectory)
    outputCSV = os.path.join(outputDirectory, replaceExtension(basename(bam), ".csv", "_tcperreadpos"))
    outputPDF = os.path.join(outputDirectory, replaceExtension(basename(bam), ".pdf", "_tcperreadpos"))
    outputLOG = os.path.join(outputDirectory, replaceExtension(basename(bam), ".log", "_tcperreadpos"))
    if(snpDirectory != None):
        inputSNP = os.path.join(snpDirectory, replaceExtension(basename(bam), ".vcf", "_snp"))
    else:
        inputSNP = None
    log = getLogFile(outputLOG)
    
    stats.tcPerReadPos(referenceFile, bam, minMQ, maxReadLength, outputCSV, outputPDF, inputSNP, log)
    
    closeLogFile(log)
    stepFinished()
    
def runSTcPerUtr(tid, bam, referenceFile, minMQ, maxReadLength, outputDirectory, snpDirectory):
    createDir(outputDirectory)
    outputCSV = os.path.join(outputDirectory, replaceExtension(basename(bam), ".csv", "_tcperutr"))
    outputPDF = os.path.join(outputDirectory, replaceExtension(basename(bam), ".pdf", "_tcperutr"))
    outputLOG = os.path.join(outputDirectory, replaceExtension(basename(bam), ".log", "_tcperutr"))
    if(snpDirectory != None):
        inputSNP = os.path.join(snpDirectory, replaceExtension(basename(bam), ".vcf", "_snp"))
    else:
        inputSNP = None
    log = getLogFile(outputLOG)
    
    stats.tcPerUtr(referenceFile, args.bed, bam, minMQ, maxReadLength, outputCSV, outputPDF, inputSNP, log, False, True, True)
    
    closeLogFile(log)
    stepFinished()

def runUtrCoverage(tid, bam, minMQ, outputDirectory):
    createDir(outputDirectory)
    outputCSV = os.path.join(outputDirectory, replaceExtension(basename(bam), ".csv", "_utrcoverage"))
    outputPDF = os.path.join(outputDirectory, replaceExtension(basename(bam), ".pdf", "_utrcoverage"))
    outputLOG = os.path.join(outputDirectory, replaceExtension(basename(bam), ".log", "_utrcoverage"))
    log = getLogFile(outputLOG)
    
    #stats.coveragePerUtr(args.bed, bam, minMQ, outputCSV, outputPDF, log, False, True, True)
    
    closeLogFile(log)
    stepFinished()

    
def runDumpReadInfo(tid, bam, referenceFile, minMQ, outputDirectory, snpDirectory):
    createDir(outputDirectory)
    outputCSV = os.path.join(outputDirectory, replaceExtension(basename(bam), ".sdunk", "_readinfo"))
    outputLOG = os.path.join(outputDirectory, replaceExtension(basename(bam), ".log", "_readinfo"))
    if(snpDirectory != None):
        inputSNP = os.path.join(snpDirectory, replaceExtension(basename(bam), ".vcf", "_snp"))
    else:
        inputSNP = None
    log = getLogFile(outputLOG)
    
    dump.dumpReadInfo(referenceFile, bam, minMQ, outputCSV, inputSNP, log)
    
    closeLogFile(log)
    stepFinished()


def runAll() :
    message("slamdunk all")
    # TODO

########################################################################
# Argument parsing
########################################################################

# Info
usage = "SLAMdunk software for analyzing SLAM-seq data"
version = "1.0"

# Main Parsers
parser = ArgumentParser(description=usage, formatter_class=RawDescriptionHelpFormatter, version=version)

# Initialize Subparsers
subparsers = parser.add_subparsers(help="", dest="command")

# map command

mapparser = subparsers.add_parser('map', help='Map SLAM-seq read data')
mapparser.add_argument('bam', action='store', help='Bam file(s)' , nargs="+")
mapparser.add_argument("-r", "--reference", type=str, required=True, dest="referenceFile", help="Reference fasta file")
mapparser.add_argument("-o", "--outputDir", type=str, required=True, dest="outputDir", help="Output directory for mapped BAM files.")
mapparser.add_argument("-5", "--trim-5p", type=int, required=False, dest="trim5", help="Number of bp removed from 5' end of all reads.")
mapparser.add_argument("-t", "--threads", type=int, required=False, dest="threads", help="Thread number")

# filter command

filterparser = subparsers.add_parser('filter', help='Filter SLAM-seq aligned data')
filterparser.add_argument('bam', action='store', help='Bam file(s)' , nargs="+")
filterparser.add_argument("-o", "--outputDir", type=str, required=True, dest="outputDir", help="Output directory for mapped BAM files.")
filterparser.add_argument("-mq", "--min-mq", type=int, required=False, default=2, dest="mq", help="Minimal mapping quality")
filterparser.add_argument("-t", "--threads", type=int, required=False, dest="threads", help="Thread number")

# snp command

snpparser = subparsers.add_parser('snp', help='Call SNPs on SLAM-seq aligned data')
snpparser.add_argument('bam', action='store', help='Bam file(s)' , nargs="+")
snpparser.add_argument("-o", "--outputDir", type=str, required=True, dest="outputDir", help="Output directory for mapped BAM files.")
snpparser.add_argument("-f", "--fasta", required=True, dest="fasta", type=str, help="Reference fasta file")
snpparser.add_argument("-c", "--min-coverage", required=False, dest="cov", type=int, help="Minimimum coverage to call variant", default=10)
snpparser.add_argument("-a", "--var-fraction", required=False, dest="var", type=float, help="Minimimum variant fraction variant", default=0.8)
snpparser.add_argument("-t", "--threads", type=int, required=False, default=1, dest="threads", help="Thread number")

# dedup command

dedupparser = subparsers.add_parser('dedup', help='Deduplicate SLAM-seq aligned data')
dedupparser.add_argument("-o", "--outputDir", type=str, required=True, dest="outputDir", help="Output directory for mapped BAM files.")
dedupparser.add_argument("-t", "--threads", type=int, required=False, default=1, dest="threads", help="Thread number")
dedupparser.add_argument('bam', action='store', help='Bam file(s)' , nargs="+")

# count command

countparser = subparsers.add_parser('count', help='Count T/C conversions in SLAM-seq aligned data')
countparser.add_argument('bam', action='store', help='Bam file(s)' , nargs="+")
countparser.add_argument("-p", "--output-prefix", type=str, required=False, default="summary", dest="outputPrefix", help="Name of output file.")
countparser.add_argument("-n", "--sample-names", type=str, required=False, dest="sampleNames", help="CSV file containing name for all samples.")
countparser.add_argument("-o", "--outputDir", type=str, required=True, dest="outputDir", help="Output directory for mapped BAM files.")
countparser.add_argument("-s", "--snp-directory", type=str, required=False, dest="snpDir", help="Directory containing SNP files.")
countparser.add_argument("-r", "--reference", type=str, required=True, dest="ref", help="Reference fasta file")
countparser.add_argument("-b", "--bed", type=str, required=True, dest="bed", help="BED file")
countparser.add_argument("-l", "--max-read-length", type=int, required=True, dest="maxLength", help="Max read length in BAM file")
countparser.add_argument("-q", "--min-base-qual", type=int, default=0, required=False, dest="minQual", help="Min base quality for T -> C conversions")
countparser.add_argument("-t", "--threads", type=int, required=False, default=1, dest="threads", help="Thread number")

# stats command

statsparser = subparsers.add_parser('stats.rates', help='Calculate stats on SLAM-seq datasets')
statsparser.add_argument('bam', action='store', help='Bam file(s)' , nargs="+")
statsparser.add_argument("-o", "--outputDir", type=str, required=True, dest="outputDir", help="Output directory for mapped BAM files.")
statsparser.add_argument("-r", "--reference", type=str, required=True, dest="referenceFile", help="Reference fasta file")
statsparser.add_argument("-mq", "--min-basequality", type=int, required=False, default=0, dest="mq", help="Minimal base quality for SNPs")
#statsparser.add_argument('-R', "--compute-rates", dest="overallRates", action='store_true', help="Compute overall conversion rates.")
statsparser.add_argument("-t", "--threads", type=int, required=False, default=1, dest="threads", help="Thread number")

# stats summary command

statsSumParser = subparsers.add_parser('stats.summary', help='Prints a CSV file containing the number of sequenced, mapped and filtered reads for all samples')
statsSumParser.add_argument("-o", "--outputPrefix", type=str, required=True, dest="outputPrefix", help="Prefix for output files")
statsSumParser.add_argument("-n", "--sample-names", type=str, required=True, dest="sampleNames", help="CSV file containing name for all samples.")
statsSumParser.add_argument("-r", "--read-counts", type=str, required=True, dest="readCounts", help="CSV file containing read counts.")
statsSumParser.add_argument("-s", "--snp-files", type=str, nargs="+", required=True, dest="snpFiles", help="SNP files for all samples")
statsSumParser.add_argument("-m", "--mapped-files", type=str, nargs="+", required=True, dest="mappedFiles", help="BAM files for all samples")
statsSumParser.add_argument("-f", "--filtered-files", type=str, nargs="+", required=False, dest="filteredFiles", help="Filtered BAM files for all samples")
statsSumParser.add_argument("-d", "--deduplicated-files", type=str, nargs="+", required=False, dest="dedupFiles", help="Deduplicated BAM files for all samples")

# stats read info command

conversionRateParser = subparsers.add_parser('stats.tcperreadpos', help='Get SlamSeq info per read')
conversionRateParser.add_argument('bam', action='store', help='Bam file(s)' , nargs="+")
conversionRateParser.add_argument("-r", "--reference", type=str, required=True, dest="referenceFile", help="Reference fasta file")
conversionRateParser.add_argument("-s", "--snp-directory", type=str, required=False, dest="snpDir", help="Directory containing SNP files.")
conversionRateParser.add_argument("-l", "--max-read-length", type=int, required=True, dest="maxLength", help="Max read length in BAM file")
conversionRateParser.add_argument("-o", "--outputDir", type=str, required=True, dest="outputDir", help="Output directory for mapped BAM files.")#conversionRateParser.add_argument("-5", "--trim-5p", type=int, required=False, dest="trim5", help="Number of bp removed from 5' end of all reads.")
conversionRateParser.add_argument("-mq", "--min-basequality", type=int, required=False, default=0, dest="mq", help="Minimal base quality for SNPs")
conversionRateParser.add_argument("-t", "--threads", type=int, required=False, dest="threads", help="Thread number")

# stats utr info command

utrRateParser = subparsers.add_parser('stats.tcperutrpos', help='Get SlamSeq info per utr')
utrRateParser.add_argument('bam', action='store', help='Bam file(s)' , nargs="+")
utrRateParser.add_argument("-r", "--reference", type=str, required=True, dest="referenceFile", help="Reference fasta file")
utrRateParser.add_argument("-b", "--bed", type=str, required=True, dest="bed", help="BED file")
utrRateParser.add_argument("-s", "--snp-directory", type=str, required=False, dest="snpDir", help="Directory containing SNP files.")
utrRateParser.add_argument("-l", "--max-read-length", type=int, required=True, dest="maxLength", help="Max read length in BAM file")
utrRateParser.add_argument("-o", "--outputDir", type=str, required=True, dest="outputDir", help="Output directory for mapped BAM files.")#conversionRateParser.add_argument("-5", "--trim-5p", type=int, required=False, dest="trim5", help="Number of bp removed from 5' end of all reads.")
utrRateParser.add_argument("-mq", "--min-basequality", type=int, required=False, default=0, dest="mq", help="Minimal base quality for SNPs")
utrRateParser.add_argument("-t", "--threads", type=int, required=False, dest="threads", help="Thread number")

# stats mean coverage for all utrs

utrCoverageParser = subparsers.add_parser('stats.utrcoverage', help='Get SlamSeq info per utr')
utrCoverageParser.add_argument('bam', action='store', help='Bam file(s)' , nargs="+")
# utrCoverageParser.add_argument("-r", "--reference", type=str, required=True, dest="referenceFile", help="Reference fasta file")
utrCoverageParser.add_argument("-b", "--bed", type=str, required=True, dest="bed", help="BED file")
# utrCoverageParser.add_argument("-s", "--snp-directory", type=str, required=False, dest="snpDir", help="Directory containing SNP files.")
# utrCoverageParser.add_argument("-l", "--max-read-length", type=int, required=True, dest="maxLength", help="Max read length in BAM file")
utrCoverageParser.add_argument("-o", "--outputDir", type=str, required=True, dest="outputDir", help="Output directory for mapped BAM files.")#conversionRateParser.add_argument("-5", "--trim-5p", type=int, required=False, dest="trim5", help="Number of bp removed from 5' end of all reads.")
utrCoverageParser.add_argument("-mq", "--min-basequality", type=int, required=False, default=0, dest="mq", help="Minimal base quality for SNPs")
utrCoverageParser.add_argument("-t", "--threads", type=int, required=False, dest="threads", help="Thread number")


# dump read info command

dumpReadInfo = subparsers.add_parser('dump.reads', help='Print all info available for reads')
dumpReadInfo.add_argument('bam', action='store', help='Bam file(s)' , nargs="+")
dumpReadInfo.add_argument("-r", "--reference", type=str, required=True, dest="referenceFile", help="Reference fasta file")
dumpReadInfo.add_argument("-s", "--snp-directory", type=str, required=False, dest="snpDir", help="Directory containing SNP files.")
dumpReadInfo.add_argument("-o", "--outputDir", type=str, required=True, dest="outputDir", help="Output directory for mapped BAM files.")#conversionRateParser.add_argument("-5", "--trim-5p", type=int, required=False, dest="trim5", help="Number of bp removed from 5' end of all reads.")
dumpReadInfo.add_argument("-mq", "--min-basequality", type=int, required=False, default=0, dest="mq", help="Minimal base quality for SNPs")
dumpReadInfo.add_argument("-t", "--threads", type=int, required=False, dest="threads", help="Thread number")


# all command

countparser = subparsers.add_parser('all', help='Run entire SLAMdunk analysis')
countparser.add_argument('bam', action='store', help='Bam file(s)' , nargs="+")

args = parser.parse_args()

########################################################################
# Routine selection
########################################################################

command = args.command

if (command == "map") :
    outputDirectory = args.outputDir
    n = args.threads
    referenceFile = args.referenceFile
    message("Running slamDunk map for " + str(len(args.bam)) + " files (" + str(n) + " threads)")
    for bam in args.bam:
        runMap(0, bam, referenceFile, n, args.trim5, outputDirectory)
    dunkFinished()
    message("Running slamDunk sort for " + str(len(args.bam)) + " files (" + str(n) + " threads)")
    results = Parallel(n_jobs=n, verbose=verbose)(delayed(runSort)(tid, args.bam[tid], outputDirectory) for tid in range(0, len(args.bam)))
    dunkFinished()
     
elif (command == "filter") :
    outputDirectory = args.outputDir
    n = args.threads
    message("Running slamDunk filter for " + str(len(args.bam)) + " files (" + str(n) + " threads)")
    results = Parallel(n_jobs=n, verbose=verbose)(delayed(runFilter)(tid, args.bam[tid], outputDirectory) for tid in range(0, len(args.bam)))
    dunkFinished()
    
elif (command == "snp") :
    outputDirectory = args.outputDir
    fasta = args.fasta
    minCov = args.cov
    minVarFreq = args.var
    n = args.threads
    if(n > 1):
        n = n / 2
    message("Running slamDunk SNP for " + str(len(args.bam)) + " files (" + str(n) + " threads)")
    results = Parallel(n_jobs=n, verbose=verbose)(delayed(runSnp)(tid, fasta, minCov, minVarFreq, args.bam[tid], outputDirectory) for tid in range(0, len(args.bam)))
    dunkFinished()

elif (command == "dedup") :
    outputDirectory = args.outputDir
    n = args.threads
    message("Running slamDunk dedup for " + str(len(args.bam)) + " files (" + str(n) + " threads)")
    results = Parallel(n_jobs=n, verbose=verbose)(delayed(runDedup)(tid, args.bam[tid], outputDirectory) for tid in range(0, len(args.bam)))
    dunkFinished()
    
elif (command == "count") :
    outputDirectory = args.outputDir
    snpDirectory = args.snpDir
    n = args.threads
    message("Running slamDunk tcount for " + str(len(args.bam)) + " files (" + str(n) + " threads)")
    results = Parallel(n_jobs=n, verbose=verbose)(delayed(runCount)(tid, args.bam[tid], outputDirectory, snpDirectory) for tid in range(0, len(args.bam)))
    #runCountCombine(results, args.sampleNames, args.outputPrefix, outputDirectory)
    dunkFinished()
    
elif (command == "stats.rates") :  
    outputDirectory = args.outputDir
    n = args.threads
    referenceFile = args.referenceFile
    minMQ = args.mq
    message("Running slamDunk stats for " + str(len(args.bam)) + " files (" + str(n) + " threads)")
    results = Parallel(n_jobs=n, verbose=verbose)(delayed(runStatsRates)(tid, args.bam[tid], referenceFile, minMQ, outputDirectory) for tid in range(0, len(args.bam)))
    dunkFinished() 
    
elif (command == "stats.summary") :
    samples = readSampleNames(args.sampleNames, None)
    n = 1
    message("Running slamDunk stats read summary for " + str(len(args.mappedFiles)) + " files (" + str(n) + " threads)")
    outputLog = args.outputPrefix + ".log"
    stats.readSummary(args.mappedFiles, args.filteredFiles, args.dedupFiles, args.snpFiles, samples, args.outputPrefix, getLogFile(outputLog))
    stats.sampleSummary(args.readCounts, args.outputPrefix, getLogFile(outputLog))
    dunkFinished() 

elif (command == "stats.tcperreadpos") :
    outputDirectory = args.outputDir
    n = args.threads
    snpDirectory = args.snpDir
    referenceFile = args.referenceFile
    minMQ = args.mq
    message("Running slamDunk stats.tcperreadpos for " + str(len(args.bam)) + " files (" + str(n) + " threads)")
    results = Parallel(n_jobs=n, verbose=verbose)(delayed(runSTcPerReadPos)(tid, args.bam[tid], referenceFile, minMQ, args.maxLength, outputDirectory, snpDirectory) for tid in range(0, len(args.bam)))
    dunkFinished()
    
elif (command == "stats.tcperutrpos") :
    
    outputDirectory = args.outputDir
    n = args.threads
    snpDirectory = args.snpDir
    referenceFile = args.referenceFile
    minMQ = args.mq
    snpDirectory = args.snpDir
    message("Running slamDunk stats.tcperutrpos for " + str(len(args.bam)) + " files (" + str(n) + " threads)")
    results = Parallel(n_jobs=n, verbose=verbose)(delayed(runSTcPerUtr)(tid, args.bam[tid], referenceFile, minMQ, args.maxLength, outputDirectory, snpDirectory) for tid in range(0, len(args.bam)))
    dunkFinished()

elif (command == "stats.utrcoverage"):

    outputDirectory = args.outputDir
    n = args.threads
#     snpDirectory = args.snpDir
#     referenceFile = args.referenceFile
    minMQ = args.mq
#     snpDirectory = args.snpDir
    message("Running slamDunk stats.utrcoverage for " + str(len(args.bam)) + " files (" + str(n) + " threads)")
    #results = Parallel(n_jobs=n, verbose=verbose)(delayed(runUtrCoverage)(tid, args.bam[tid], minMQ, outputDirectory) for tid in range(0, len(args.bam)))
    dunkFinished()

elif (command == "dump.reads") :
    outputDirectory = args.outputDir
    n = args.threads
    snpDirectory = args.snpDir
    referenceFile = args.referenceFile
    minMQ = args.mq
    message("Running slamDunk dump for " + str(len(args.bam)) + " files (" + str(n) + " threads)")
    results = Parallel(n_jobs=n, verbose=verbose)(delayed(runDumpReadInfo)(tid, args.bam[tid], referenceFile, minMQ, outputDirectory, snpDirectory) for tid in range(0, len(args.bam)))
    dunkFinished()

    
elif (command == "all") :
    runAll()
    dunkFinished()
    
#########################################################################
# Cleanup
########################################################################
    
sys.exit(0)
