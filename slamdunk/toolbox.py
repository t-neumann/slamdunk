#!/usr/bin/env python

#ALLEYOOP:

#Additional sLamdunk heLpEr tools for anY diagnOstics Or Plots


#########################################################################
# Main routine for the SLAMdunk simulation
#########################################################################
# Imports
#########################################################################

from __future__ import print_function
import sys, os, random

from argparse import ArgumentParser, RawDescriptionHelpFormatter
    
from os.path import basename

from joblib import Parallel, delayed
from dunks import tcounter, mapper, filter, deduplicator, stats, snps, dump
from utils.misc import replaceExtension, readSampleNames

########################################################################
# Global variables
########################################################################

printOnly = False
verbose = True

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
            
def runDedup(tid, bam, outputDirectory) :
    outputBAM = os.path.join(outputDirectory, replaceExtension(basename(bam), ".bam", "_dedup"))
    outputLOG = os.path.join(outputDirectory, replaceExtension(basename(bam), ".log", "_dedup"))
    log = getLogFile(outputLOG)
    deduplicator.Dedup(bam, outputBAM, log)
    closeLogFile(log)
    stepFinished()
    
def runHalfLifes(bams, timepoints, outputDirectory) :
    outputCSV = os.path.join(outputDirectory, replaceExtension(basename(bams[0]), ".tsv", "_halflifes"))
    outputLOG = os.path.join(outputDirectory, replaceExtension(basename(bams[0]), ".log", "_halflifes"))
    log = getLogFile(outputLOG)
    stats.halflifes(",".join(bams), outputCSV, timepoints, log)
    closeLogFile(log)
    stepFinished()
    
def runStatsRates(tid, bam, referenceFile, minMQ, outputDirectory) :
    outputCSV = os.path.join(outputDirectory, replaceExtension(basename(bam), ".csv", "_overallrates"))
    outputPDF = os.path.join(outputDirectory, replaceExtension(basename(bam), ".pdf", "_overallrates"))
    outputLOG = os.path.join(outputDirectory, replaceExtension(basename(bam), ".log", "_overallrates"))
    log = getLogFile(outputLOG)
    stats.statsComputeOverallRates(referenceFile, bam, minMQ, outputCSV, outputPDF, log)
    closeLogFile(log)
    stepFinished()
    
def runStatsTCContext(tid, bam, referenceFile, minMQ, outputDirectory) :
    outputCSV = os.path.join(outputDirectory, replaceExtension(basename(bam), ".csv", "_tccontext"))
    outputPDF = os.path.join(outputDirectory, replaceExtension(basename(bam), ".pdf", "_tccontext"))
    outputLOG = os.path.join(outputDirectory, replaceExtension(basename(bam), ".log", "_tccontext"))
    log = getLogFile(outputLOG)
    stats.statsComputeTCContext(referenceFile, bam, minMQ, outputCSV, outputPDF, log)
    closeLogFile(log)
    stepFinished()

def runStatsRatesUTR(tid, bam, referenceFile, minMQ, outputDirectory, utrFile, maxReadLength) :
    outputCSV = os.path.join(outputDirectory, replaceExtension(basename(bam), ".csv", "_mutationrates_utr"))
    outputPDF = os.path.join(outputDirectory, replaceExtension(basename(bam), ".pdf", "_mutationrates_utr"))
    outputLOG = os.path.join(outputDirectory, replaceExtension(basename(bam), ".log", "_mutationrates_utr"))
    log = getLogFile(outputLOG)
    stats.statsComputeOverallRatesPerUTR(referenceFile, bam, minMQ, outputCSV, outputPDF, utrFile, maxReadLength, log)
    closeLogFile(log)
    stepFinished()

    
def runSTcPerReadPos(tid, bam, referenceFile, minMQ, maxReadLength, outputDirectory, snpDirectory):
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
    
def runSTcPerUtr(tid, bam, referenceFile, bed, minMQ, maxReadLength, outputDirectory, snpDirectory):
    outputCSV = os.path.join(outputDirectory, replaceExtension(basename(bam), ".csv", "_tcperutr"))
    outputPDF = os.path.join(outputDirectory, replaceExtension(basename(bam), ".pdf", "_tcperutr"))
    outputLOG = os.path.join(outputDirectory, replaceExtension(basename(bam), ".log", "_tcperutr"))
    if(snpDirectory != None):
        inputSNP = os.path.join(snpDirectory, replaceExtension(basename(bam), ".vcf", "_snp"))
    else:
        inputSNP = None
    log = getLogFile(outputLOG)
    
    stats.tcPerUtr(referenceFile, bed, bam, minMQ, maxReadLength, outputCSV, outputPDF, inputSNP, log, False, True, True)
    
    closeLogFile(log)
    stepFinished()

def runUtrCoverage(tid, bam, minMQ, outputDirectory):
    outputCSV = os.path.join(outputDirectory, replaceExtension(basename(bam), ".csv", "_utrcoverage"))
    outputPDF = os.path.join(outputDirectory, replaceExtension(basename(bam), ".pdf", "_utrcoverage"))
    outputLOG = os.path.join(outputDirectory, replaceExtension(basename(bam), ".log", "_utrcoverage"))
    log = getLogFile(outputLOG)
    
    #stats.coveragePerUtr(args.bed, bam, minMQ, outputCSV, outputPDF, log, False, True, True)
    
    closeLogFile(log)
    stepFinished()

    
def runDumpReadInfo(tid, bam, referenceFile, minMQ, outputDirectory, snpDirectory):
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


def run():
    
    ########################################################################
    # Argument parsing
    ########################################################################
    
    # Info
    usage = "AlleyOop utility tools and diagnostics for SLAMSeq data"
    version = "0.1.0"
    
    # Main Parsers
    parser = ArgumentParser(description=usage, formatter_class=RawDescriptionHelpFormatter, version=version)
    
    # Initialize Subparsers
    subparsers = parser.add_subparsers(help="", dest="command")
    
    # dedup command
    dedupparser = subparsers.add_parser('dedup', help='Deduplicate SLAM-seq aligned data')
    dedupparser.add_argument("-o", "--outputDir", type=str, required=True, dest="outputDir", help="Output directory for mapped BAM files.")
    dedupparser.add_argument("-t", "--threads", type=int, required=False, default=1, dest="threads", help="Thread number")
    dedupparser.add_argument('bam', action='store', help='Bam file(s)' , nargs="+")
    
    # stats command
    statsparser = subparsers.add_parser('stats.rates', help='Calculate overall conversion rates on SLAM-seq datasets')
    statsparser.add_argument('bam', action='store', help='Bam file(s)' , nargs="+")
    statsparser.add_argument("-o", "--outputDir", type=str, required=True, dest="outputDir", help="Output directory for mapped BAM files.")
    statsparser.add_argument("-r", "--reference", type=str, required=True, dest="referenceFile", help="Reference fasta file")
    statsparser.add_argument("-mq", "--min-basequality", type=int, required=False, default=0, dest="mq", help="Minimal base quality for SNPs")
    #statsparser.add_argument('-R', "--compute-rates", dest="overallRates", action='store_true', help="Compute overall conversion rates.")
    statsparser.add_argument("-t", "--threads", type=int, required=False, default=1, dest="threads", help="Thread number")
    
    # context command
    tccontextparser = subparsers.add_parser('stats.TCcontext', help='Calculate T->C conversion context on SLAM-seq datasets')
    tccontextparser.add_argument('bam', action='store', help='Bam file(s)' , nargs="+")
    #tccontextparser.add_argument("-b", "--bed", type=str, required=True, dest="bed", help="BED file")
    tccontextparser.add_argument("-o", "--outputDir", type=str, required=True, dest="outputDir", help="Output directory for mapped BAM files.")
    tccontextparser.add_argument("-r", "--reference", type=str, required=True, dest="referenceFile", help="Reference fasta file")
    tccontextparser.add_argument("-mq", "--min-basequality", type=int, required=False, default=0, dest="mq", help="Minimal base quality for SNPs")
    tccontextparser.add_argument("-t", "--threads", type=int, required=False, default=1, dest="threads", help="Thread number")

    # stats rates utr command
    statsutrrateparser = subparsers.add_parser('stats.utrrates', help='Calculate conversion rates per UTR on SLAM-seq datasets')
    statsutrrateparser.add_argument('bam', action='store', help='Bam file(s)' , nargs="+")
    statsutrrateparser.add_argument("-o", "--outputDir", type=str, required=True, dest="outputDir", help="Output directory for mapped BAM files.")
    statsutrrateparser.add_argument("-r", "--reference", type=str, required=True, dest="referenceFile", help="Reference fasta file")
    statsutrrateparser.add_argument("-mq", "--min-basequality", type=int, required=False, default=0, dest="mq", help="Minimal base quality for SNPs")
    statsutrrateparser.add_argument("-t", "--threads", type=int, required=False, default=1, dest="threads", help="Thread number")
    statsutrrateparser.add_argument("-b", "--bed", type=str, required=True, dest="bed", help="BED file")
    statsutrrateparser.add_argument("-l", "--max-read-length", type=int, required=True, dest="maxLength", help="Max read length in BAM file")
    
    # stats summary command
    statsSumParser = subparsers.add_parser('stats.summary', help='Display summary information and statistics on read numbers')
    statsSumParser.add_argument("-o", "--outputPrefix", type=str, required=True, dest="outputPrefix", help="Prefix for output files")
    statsSumParser.add_argument("-n", "--sample-names", type=str, required=True, dest="sampleNames", help="CSV file containing name for all samples.")
    statsSumParser.add_argument("-r", "--read-counts", type=str, required=True, dest="readCounts", help="CSV file containing read counts.")
    statsSumParser.add_argument("-s", "--snp-files", type=str, nargs="+", required=True, dest="snpFiles", help="SNP files for all samples")
    statsSumParser.add_argument("-m", "--mapped-files", type=str, nargs="+", required=True, dest="mappedFiles", help="BAM files for all samples")
    statsSumParser.add_argument("-f", "--filtered-files", type=str, nargs="+", required=False, dest="filteredFiles", help="Filtered BAM files for all samples")
    statsSumParser.add_argument("-d", "--deduplicated-files", type=str, nargs="+", required=False, dest="dedupFiles", help="Deduplicated BAM files for all samples")
    
    # stats read info command
    conversionRateParser = subparsers.add_parser('stats.tcperreadpos', help='Calculate conversion rates per read position on SLAM-seq datasets')
    conversionRateParser.add_argument('bam', action='store', help='Bam file(s)' , nargs="+")
    conversionRateParser.add_argument("-r", "--reference", type=str, required=True, dest="referenceFile", help="Reference fasta file")
    conversionRateParser.add_argument("-s", "--snp-directory", type=str, required=False, dest="snpDir", help="Directory containing SNP files.")
    conversionRateParser.add_argument("-l", "--max-read-length", type=int, required=True, dest="maxLength", help="Max read length in BAM file")
    conversionRateParser.add_argument("-o", "--outputDir", type=str, required=True, dest="outputDir", help="Output directory for mapped BAM files.")#conversionRateParser.add_argument("-5", "--trim-5p", type=int, required=False, dest="trim5", help="Number of bp removed from 5' end of all reads.")
    conversionRateParser.add_argument("-mq", "--min-basequality", type=int, required=False, default=0, dest="mq", help="Minimal base quality for SNPs")
    conversionRateParser.add_argument("-t", "--threads", type=int, required=False, dest="threads", help="Thread number")
    
    # stats utr info command
    utrRateParser = subparsers.add_parser('stats.tcperutrpos', help='Calculate conversion rates per UTR position on SLAM-seq datasets')
    utrRateParser.add_argument('bam', action='store', help='Bam file(s)' , nargs="+")
    utrRateParser.add_argument("-r", "--reference", type=str, required=True, dest="referenceFile", help="Reference fasta file")
    utrRateParser.add_argument("-b", "--bed", type=str, required=True, dest="bed", help="BED file")
    utrRateParser.add_argument("-s", "--snp-directory", type=str, required=False, dest="snpDir", help="Directory containing SNP files.")
    utrRateParser.add_argument("-l", "--max-read-length", type=int, required=True, dest="maxLength", help="Max read length in BAM file")
    utrRateParser.add_argument("-o", "--outputDir", type=str, required=True, dest="outputDir", help="Output directory for mapped BAM files.")#conversionRateParser.add_argument("-5", "--trim-5p", type=int, required=False, dest="trim5", help="Number of bp removed from 5' end of all reads.")
    utrRateParser.add_argument("-mq", "--min-basequality", type=int, required=False, default=0, dest="mq", help="Minimal base quality for SNPs")
    utrRateParser.add_argument("-t", "--threads", type=int, required=False, dest="threads", help="Thread number")
    
    # stats mean coverage for all utrs
    utrCoverageParser = subparsers.add_parser('stats.utrcoverage', help='Calculate UTR coverage on SLAM-seq datasets')
    utrCoverageParser.add_argument('bam', action='store', help='Bam file(s)' , nargs="+")
    # utrCoverageParser.add_argument("-r", "--reference", type=str, required=True, dest="referenceFile", help="Reference fasta file")
    utrCoverageParser.add_argument("-b", "--bed", type=str, required=True, dest="bed", help="BED file")
    # utrCoverageParser.add_argument("-s", "--snp-directory", type=str, required=False, dest="snpDir", help="Directory containing SNP files.")
    # utrCoverageParser.add_argument("-l", "--max-read-length", type=int, required=True, dest="maxLength", help="Max read length in BAM file")
    utrCoverageParser.add_argument("-o", "--outputDir", type=str, required=True, dest="outputDir", help="Output directory for mapped BAM files.")#conversionRateParser.add_argument("-5", "--trim-5p", type=int, required=False, dest="trim5", help="Number of bp removed from 5' end of all reads.")
    utrCoverageParser.add_argument("-mq", "--min-basequality", type=int, required=False, default=0, dest="mq", help="Minimal base quality for SNPs")
    utrCoverageParser.add_argument("-t", "--threads", type=int, required=False, dest="threads", help="Thread number")
    
    # dump read info command
    dumpReadInfo = subparsers.add_parser('dump.reads', help='Print all info available for slamdunk reads')
    dumpReadInfo.add_argument('bam', action='store', help='Bam file(s)' , nargs="+")
    dumpReadInfo.add_argument("-r", "--reference", type=str, required=True, dest="referenceFile", help="Reference fasta file")
    dumpReadInfo.add_argument("-s", "--snp-directory", type=str, required=False, dest="snpDir", help="Directory containing SNP files.")
    dumpReadInfo.add_argument("-o", "--outputDir", type=str, required=True, dest="outputDir", help="Output directory for mapped BAM files.")#conversionRateParser.add_argument("-5", "--trim-5p", type=int, required=False, dest="trim5", help="Number of bp removed from 5' end of all reads.")
    dumpReadInfo.add_argument("-mq", "--min-basequality", type=int, required=False, default=0, dest="mq", help="Minimal base quality for SNPs")
    dumpReadInfo.add_argument("-t", "--threads", type=int, required=False, dest="threads", help="Thread number")
    
    args = parser.parse_args()
    
    ########################################################################
    # Routine selection
    ########################################################################
    
    command = args.command

    if (command == "dedup") :
        outputDirectory = args.outputDir
        createDir(outputDirectory)
        n = args.threads
        message("Running alleyoop dedup for " + str(len(args.bam)) + " files (" + str(n) + " threads)")
        results = Parallel(n_jobs=n, verbose=verbose)(delayed(runDedup)(tid, args.bam[tid], outputDirectory) for tid in range(0, len(args.bam)))
        dunkFinished()
        
    elif (command == "half-lifes") :
        
        outputDirectory = args.outputDir
        createDir(outputDirectory)
        
        timepoints = args.timepoints
        
        message("Running alleyoop half-lifes for " + str(len(args.bam)) + " files")
        runHalfLifes(args.bam, timepoints, outputDirectory)
        #results = Parallel(n_jobs=n, verbose=verbose)(delayed(runHalfLifes)(tid, args.bam[tid], timepoints, outputDirectory) for tid in range(0, len(args.bam)))
        dunkFinished()
        
    elif (command == "stats.rates") :  
        outputDirectory = args.outputDir
        createDir(outputDirectory)
        n = args.threads
        referenceFile = args.referenceFile
        minMQ = args.mq
        message("Running alleyoop stats for " + str(len(args.bam)) + " files (" + str(n) + " threads)")
        results = Parallel(n_jobs=n, verbose=verbose)(delayed(runStatsRates)(tid, args.bam[tid], referenceFile, minMQ, outputDirectory) for tid in range(0, len(args.bam)))
        dunkFinished()
    
    elif (command == "stats.TCcontext") :  
        outputDirectory = args.outputDir
        createDir(outputDirectory)
        n = args.threads
        referenceFile = args.referenceFile
        minMQ = args.mq
        message("Running alleyoop TC context for " + str(len(args.bam)) + " files (" + str(n) + " threads)")
        results = Parallel(n_jobs=n, verbose=verbose)(delayed(runStatsTCContext)(tid, args.bam[tid], referenceFile, minMQ, outputDirectory) for tid in range(0, len(args.bam)))
        dunkFinished()
    
    elif (command == "stats.utrrates") :  
        outputDirectory = args.outputDir
        createDir(outputDirectory)
        n = args.threads
        referenceFile = args.referenceFile
        minMQ = args.mq
        
        message("Running alleyoop stats for " + str(len(args.bam)) + " files (" + str(n) + " threads)")
        results = Parallel(n_jobs=n, verbose=verbose)(delayed(runStatsRatesUTR)(tid, args.bam[tid], referenceFile, minMQ, outputDirectory, args.bed, args.maxLength) for tid in range(0, len(args.bam)))
        dunkFinished()
    
    elif (command == "stats.summary") :
        samples = readSampleNames(args.sampleNames, None)
        n = 1
        message("Running alleyoop stats read summary for " + str(len(args.mappedFiles)) + " files (" + str(n) + " threads)")
        outputLog = args.outputPrefix + ".log"
        stats.readSummary(args.mappedFiles, args.filteredFiles, args.dedupFiles, args.snpFiles, samples, args.outputPrefix, getLogFile(outputLog))
        stats.sampleSummary(args.readCounts, args.outputPrefix, getLogFile(outputLog))
        dunkFinished() 
    
    elif (command == "stats.tcperreadpos") :
        outputDirectory = args.outputDir
        createDir(outputDirectory)
        n = args.threads
        snpDirectory = args.snpDir
        referenceFile = args.referenceFile
        minMQ = args.mq
        message("Running alleyoop stats.tcperreadpos for " + str(len(args.bam)) + " files (" + str(n) + " threads)")
        results = Parallel(n_jobs=n, verbose=verbose)(delayed(runSTcPerReadPos)(tid, args.bam[tid], referenceFile, minMQ, args.maxLength, outputDirectory, snpDirectory) for tid in range(0, len(args.bam)))
        dunkFinished()
        
    elif (command == "stats.tcperutrpos") :
        outputDirectory = args.outputDir
        createDir(outputDirectory)
        n = args.threads
        snpDirectory = args.snpDir
        referenceFile = args.referenceFile
        minMQ = args.mq
        snpDirectory = args.snpDir
        message("Running alleyoop stats.tcperutrpos for " + str(len(args.bam)) + " files (" + str(n) + " threads)")
        results = Parallel(n_jobs=n, verbose=verbose)(delayed(runSTcPerUtr)(tid, args.bam[tid], referenceFile, args.bed, minMQ, args.maxLength, outputDirectory, snpDirectory) for tid in range(0, len(args.bam)))
        dunkFinished()
    
    elif (command == "stats.utrcoverage"):
        outputDirectory = args.outputDir
        createDir(outputDirectory)
        n = args.threads
    #     snpDirectory = args.snpDir
    #     referenceFile = args.referenceFile
        minMQ = args.mq
    #     snpDirectory = args.snpDir
        message("Running alleyoop stats.utrcoverage for " + str(len(args.bam)) + " files (" + str(n) + " threads)")
        #results = Parallel(n_jobs=n, verbose=verbose)(delayed(runUtrCoverage)(tid, args.bam[tid], minMQ, outputDirectory) for tid in range(0, len(args.bam)))
        dunkFinished()
    
    elif (command == "dump.reads") :
        outputDirectory = args.outputDir
        createDir(outputDirectory)
        n = args.threads
        snpDirectory = args.snpDir
        referenceFile = args.referenceFile
        minMQ = args.mq
        message("Running alleyoop dump.reads for " + str(len(args.bam)) + " files (" + str(n) + " threads)")
        results = Parallel(n_jobs=n, verbose=verbose)(delayed(runDumpReadInfo)(tid, args.bam[tid], referenceFile, minMQ, outputDirectory, snpDirectory) for tid in range(0, len(args.bam)))
        dunkFinished()
    
if __name__ == '__main__':
    run()