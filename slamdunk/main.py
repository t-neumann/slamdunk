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
import csv

from joblib import Parallel, delayed
from dunks import tcounter, mapper, filter, stats, snps
from dunks.utils import replaceExtension, removeExtension

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

def runMap(tid, inputBAM, referneceFile, threads, trim5p, outputDirectory) :
    outputSAM = os.path.join(outputDirectory, replaceExtension(basename(inputBAM), ".sam", "_slamdunk_mapped"))
    outputLOG = os.path.join(outputDirectory, replaceExtension(basename(inputBAM), ".log", "_slamdunk_mapped"))
    mapper.Map(inputBAM, referneceFile, outputSAM, getLogFile(outputLOG), threads=threads, trim5p=trim5p, printOnly=printOnly, verbose=verbose)
    stepFinished()

def runSort(tid, bam, outputDirectory):
    inputSAM = os.path.join(outputDirectory, replaceExtension(basename(bam), ".sam", "_slamdunk_mapped"))
    outputBAM = os.path.join(outputDirectory, replaceExtension(basename(bam), ".bam", "_slamdunk_mapped"))
    outputLOG = os.path.join(outputDirectory, replaceExtension(basename(bam), ".log", "_slamdunk_mapped"))
    mapper.sort(inputSAM, outputBAM, getLogFile(outputLOG), False, printOnly, verbose)
    stepFinished()

def runDedup() :
    message("slamdunk dedup")
    # TODO
        
def runFilter(tid, bam, outputDirectory):
    outputBAM = os.path.join(outputDirectory, replaceExtension(basename(bam), ".bam", "_filtered"))
    outputLOG = os.path.join(outputDirectory, replaceExtension(basename(bam), ".log", "_filtered"))
    filter.Filter(bam, outputBAM, getLogFile(outputLOG), args.mq, printOnly, verbose)
    stepFinished()

def runSnp(tid, referenceFile, minCov, minVarFreq, inputBAM, outputDirectory) :
    outputSNP = os.path.join(outputDirectory, replaceExtension(basename(bam), ".txt", "_snp"))
    outputLOG = os.path.join(outputDirectory, replaceExtension(basename(bam), ".log", "_snp"))
    snps.SNPs(inputBAM, outputSNP, referenceFile, minVarFreq, minCov, getLogFile(outputLOG), printOnly, verbose, True)
    stepFinished()
                
def runCount(tid, bam, outputDirectory, snpDirectory) :
    outputCSV = os.path.join(outputDirectory, replaceExtension(basename(bam), ".csv", "_tcount"))
    outputLOG = os.path.join(outputDirectory, replaceExtension(basename(bam), ".log", "_tcount"))
    inputSNP = os.path.join(snpDirectory, replaceExtension(basename(bam), ".txt", "_snp"))
    #tcounter.count(args.ref, args.bed, inputSNP, bam, args.maxLength, args.minQual, outputCSV, getLogFile(outputLOG))
    stepFinished()
    return outputCSV
    
def readSampleNames(sampleNames):
    samples = {}
    
    with open(sampleNames, "r") as sampleFile:
        samplesReader = csv.reader(sampleFile, delimiter='\t')
        for row in samplesReader:
            samples[removeExtension(row[0])] = row[1]
    return samples

def runCountCombine(bams, sampleNames, outputPrefix, outputDirectory):
    
    samples = readSampleNames(sampleNames)
    
    NON_TC_READ_COUNT = 4
    TC_READ_COUNT = 5
    TC_READ_PERC = 6
    outputFile = os.path.join(outputDirectory, outputPrefix + "_non_tc_counts.csv")
    tcounter.summary(bams, samples, outputFile, NON_TC_READ_COUNT)
    outputFile = os.path.join(outputDirectory, outputPrefix + "_tc_counts.csv")
    tcounter.summary(bams, samples, outputFile, TC_READ_COUNT)
    outputFile = os.path.join(outputDirectory, outputPrefix + "_tc_percentage.csv")
    tcounter.summary(bams, samples, outputFile, TC_READ_PERC)
    stepFinished()
        
        
def runStats(tid, bam, referenceFile, minMQ, outputDirectory, computeOverallRates) :
    if(computeOverallRates):
        outputCSV = os.path.join(outputDirectory, replaceExtension(basename(bam), ".csv", "_tcount_overallrates"))
        outputPDF = os.path.join(outputDirectory, replaceExtension(basename(bam), ".pdf", "_tcount_overallrates"))
        outputLOG = os.path.join(outputDirectory, replaceExtension(basename(bam), ".log", "_tcount_overallrates"))
        log = getLogFile(outputLOG)
        stats.statsComputeOverallRates(referenceFile, bam, minMQ, outputCSV, outputPDF, log)
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
filterparser.add_argument("-mq", "--min-mq", type=int, required=False, default=0, dest="mq", help="Minimal mapping quality")
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
dedupparser.add_argument('bam', action='store', help='Bam file(s)' , nargs="+")

# count command

countparser = subparsers.add_parser('count', help='Count T/C conversions in SLAM-seq aligned data')
countparser.add_argument('bam', action='store', help='Bam file(s)' , nargs="+")
countparser.add_argument("-p", "--ouput-prefix", type=str, required=False, default="summary", dest="outputPrefix", help="Name of output file.")
countparser.add_argument("-n", "--sample-names", type=str, required=False, dest="sampleNames", help="CSV file containing name for all samples.")
countparser.add_argument("-o", "--outputDir", type=str, required=True, dest="outputDir", help="Output directory for mapped BAM files.")
countparser.add_argument("-s", "--snp-directory", type=str, required=False, dest="snpDir", help="Directory containing SNP files.")
countparser.add_argument("-r", "--reference", type=str, required=True, dest="ref", help="Reference fasta file")
countparser.add_argument("-b", "--bed", type=str, required=True, dest="bed", help="BED file")
countparser.add_argument("-l", "--max-read-length", type=int, required=True, dest="maxLength", help="Max read length in BAM file")
countparser.add_argument("-q", "--min-base-qual", type=int, default=0, required=False, dest="minQual", help="Min base quality for T -> C conversions")
countparser.add_argument("-t", "--threads", type=int, required=False, default=1, dest="threads", help="Thread number")

# stats command

statsparser = subparsers.add_parser('stats', help='Calculate stats on SLAM-seq datasets')
statsparser.add_argument('bam', action='store', help='Bam file(s)' , nargs="+")
statsparser.add_argument("-o", "--outputDir", type=str, required=True, dest="outputDir", help="Output directory for mapped BAM files.")
statsparser.add_argument("-r", "--reference", type=str, required=True, dest="referenceFile", help="Reference fasta file")
statsparser.add_argument("-mq", "--min-mq", type=int, required=False, default=2, dest="mq", help="Minimal mapping quality")
statsparser.add_argument('-R', "--compute-rates", dest="overallRates", action='store_true', help="Compute overall conversion rates.")
statsparser.add_argument("-t", "--threads", type=int, required=False, default=1, dest="threads", help="Thread number")

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
    message("Running slamDunk map for " + str(len(args.bam)) + " files (" + str(n) + " threads)")
    for bam in args.bam:
        runMap(0, bam, outputDirectory)
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
    runDedup()
    dunkFinished()
    
elif (command == "count") :
    outputDirectory = args.outputDir
    snpDirectory = args.snpDir
    n = args.threads
    message("Running slamDunk tcount for " + str(len(args.bam)) + " files (" + str(n) + " threads)")
    results = Parallel(n_jobs=n, verbose=verbose)(delayed(runCount)(tid, args.bam[tid], outputDirectory, snpDirectory) for tid in range(0, len(args.bam)))
    runCountCombine(results, args.sampleNames, args.outputPrefix, outputDirectory)
    dunkFinished()
    
elif (command == "stats") :  
    outputDirectory = args.outputDir
    n = args.threads
    referenceFile = args.referenceFile
    minMQ = args.mq
    computeOverallRates = args.overallRates
    message("Running slamDunk stats for " + str(len(args.bam)) + " files (" + str(n) + " threads)")
    results = Parallel(n_jobs=n, verbose=verbose)(delayed(runStats)(tid, args.bam[tid], referenceFile, minMQ, outputDirectory, computeOverallRates) for tid in range(0, len(args.bam)))
    dunkFinished()  
        
elif (command == "all") :
    runAll()
    dunkFinished()
    
#########################################################################
# Cleanup
########################################################################
    
sys.exit(0)
