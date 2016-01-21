#!/usr/bin/env python

#########################################################################
# Main routine for the SLAMdunk analyzer
#########################################################################
# Imports
#########################################################################

from __future__ import print_function
import sys, os, subprocess
import tempfile

from argparse import ArgumentParser, RawDescriptionHelpFormatter
    
from os.path import basename

from joblib import Parallel, delayed
from dunks import tcounter, mapper, filter, stats
from dunks.utils import files_exist, replaceExtension, run
########################################################################
# Global variables
########################################################################

logFile = "slamdunk.log"

# Open log
log = open(logFile, 'w')

projectPath = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
varScanPath = os.path.join(projectPath, "bin", "VarScan.v2.4.1.jar")


printOnly = False
verbose = False

########################################################################
# Routine definitions
########################################################################

def error(msg, code=-1):
    print(msg)
    sys.exit(code)

def runMap() :
    print("slamdunk map", end="")
    
    outputDirectory = args.outputDir
    if(files_exist(outputDirectory)):
        # Map
        for bam in args.bam:
            outputSAM = os.path.join(outputDirectory, replaceExtension(basename(bam), ".sam", "slamdunk_mapped"))
            mapper.Map(bam, args.referenceFile, outputSAM, threads=args.threads, trim5p=args.trim5, printOnly=printOnly, verbose=verbose)
            
        # Sort
        # n   args.threads
        # results = Parallel(n_jobs=n, verbose=True)(delayed(slamSeqAnalysisMt)(tid, ref, bedUTR, files[tid], outputBase) for tid in range(0, len(files)))
        for bam in args.bam:
            inputSAM = os.path.join(outputDirectory, replaceExtension(basename(bam), ".sam", "slamdunk_mapped"))
            outputBAM = os.path.join(outputDirectory, replaceExtension(basename(bam), ".bam", "slamdunk_mapped"))
            mapper.sort(inputSAM, outputBAM, False, printOnly, verbose)
    else:
        error("Output directory doesn't exist.")
    print()
        
def runFilter(tid, bam, outputDirectory):
    outputBAM = os.path.join(outputDirectory, replaceExtension(basename(bam), ".bam", "_filtered"))
    filter.Filter(bam, outputBAM, args.mq, printOnly, verbose)

def runSnp(tid, fasta, minCov, minVarFreq, bam, outputDirectory) :
            
    outputSNP = os.path.join(outputDirectory, replaceExtension(basename(bam), ".txt", "_snp"))
    
    fileSNP = open(outputSNP, 'w')
    
    mpileupCmd = "samtools mpileup -f" + fasta + " " + bam
    mpileup = subprocess.Popen(mpileupCmd, shell=True, stdout=subprocess.PIPE, stderr=sys.stderr)
        
    varscanCmd = "java -jar " + varScanPath + " mpileup2snp  --strand-filter 0 --min-var-freq " + str(minVarFreq) + " --min-coverage " + str(minCov) + " --variants 1"
    varscan = subprocess.Popen(varscanCmd, shell=True, stdin=mpileup.stdout, stdout=fileSNP, stderr=sys.stderr)
    varscan.wait()
    
    fileSNP.close()
        
def runDedup() :
    print("slamdunk dedup", end="")
    for bam in args.bam :
        print(" " + bam, end="")
    print()
        
def runCount(tid, bam, outputDirectory) :
    outputCSV = os.path.join(outputDirectory, replaceExtension(basename(bam), ".csv", "_tcount"))
    inputSNP = os.path.join(outputDirectory, replaceExtension(basename(bam), ".txt", "_snp"))
    tcounter.count(args.ref, args.bed, inputSNP, bam, args.maxLength, args.minQual, outputCSV)
        
def runStats(tid, bam, referenceFile, minMQ, outputDirectory, computeOverallRates) :
    if(computeOverallRates):
        outputCSV = os.path.join(outputDirectory, replaceExtension(basename(bam), ".csv", "_tcount_overallrates"))
        outputPDF = os.path.join(outputDirectory, replaceExtension(basename(bam), ".pdf", "_tcount_overallrates"))
        stats.statsComputeOverallRates(referenceFile, bam, minMQ, outputCSV, outputPDF)
        
def runAll() :
    print("slamdunk all", end="")
    for bam in args.bam :
        print(" " + bam, end="")
    print()

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

# mapparser.add_argument("-i", "--input", type=str, required=True, dest="bam", help="BAM file" )

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
countparser.add_argument("-o", "--outputDir", type=str, required=True, dest="outputDir", help="Output directory for mapped BAM files.")
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
    runMap()
elif (command == "filter") :
    outputDirectory = args.outputDir
    n = args.threads
    print("Running slamDunk filter for " + str(len(args.bam)) + " files (" + str(n) + " threads)")
    results = Parallel(n_jobs=n, verbose=True)(delayed(runFilter)(tid, args.bam[tid], outputDirectory) for tid in range(0, len(args.bam)))
elif (command == "snp") :
    outputDirectory = args.outputDir
    fasta = args.fasta
    minCov = args.cov
    minVarFreq = args.var
    n = args.threads
    print("Running slamDunk SNP for " + str(len(args.bam)) + " files (" + str(n) + " threads)")
    results = Parallel(n_jobs=n, verbose=True)(delayed(runSnp)(tid, fasta, minCov, minVarFreq, args.bam[tid], outputDirectory) for tid in range(0, len(args.bam)))

elif (command == "dedup") :
    runDedup()
elif (command == "count") :
    outputDirectory = args.outputDir
    n = args.threads
    print("Running slamDunk tcount for " + str(len(args.bam)) + " files (" + str(n) + " threads)")
    results = Parallel(n_jobs=n, verbose=True)(delayed(runCount)(tid, args.bam[tid], outputDirectory) for tid in range(0, len(args.bam)))
elif (command == "stats") :
  
    outputDirectory = args.outputDir
    n = args.threads
    referenceFile = args.referenceFile
    minMQ = args.mq
    computeOverallRates = args.overallRates
    print("Running slamDunk stats for " + str(len(args.bam)) + " files (" + str(n) + " threads)")
    results = Parallel(n_jobs=n, verbose=True)(delayed(runStats)(tid, args.bam[tid], referenceFile, minMQ, outputDirectory, computeOverallRates) for tid in range(0, len(args.bam)))  
        
elif (command == "all") :
    runAll()
    
#########################################################################
# Cleanup
########################################################################

log.close()
    
sys.exit(0)
