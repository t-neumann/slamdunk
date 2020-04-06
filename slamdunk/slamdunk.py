#!/usr/bin/env python

# Copyright (c) 2015 Tobias Neumann, Philipp Rescheneder.
#
# This file is part of Slamdunk.
#
# Slamdunk is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as
# published by the Free Software Foundation, either version 3 of the
# License, or (at your option) any later version.
#
# Slamdunk is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Affero General Public License for more details.
#
# You should have received a copy of the GNU Affero General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

#########################################################################
# Main routine for the SLAMdunk analyzer
#########################################################################
# Imports
#########################################################################
from __future__ import print_function
import sys, os, random

from time import sleep
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter, SUPPRESS

from os.path import basename

from joblib import Parallel, delayed
from slamdunk.dunks import tcounter, mapper, filter, deduplicator, snps
from slamdunk.utils.misc import replaceExtension, estimateMaxReadLength
from slamdunk.version import __version__

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

def readSampleFile(fileName):
    samples = []
    infos = []

    with open(fileName, "r") as ins:
        for line in ins:
            line = line.strip()
            if(len(line) > 1):
                if(fileName.endswith(".tsv")):
                    cols = line.split("\t")
                elif(fileName.endswith(".csv")):
                    cols = line.split(",")
                else:
                    raise RuntimeError("Unknown file extension found: " + fileName)
                if(len(cols) < 4):
                    raise RuntimeError("Invalid sample file found: " + fileName)
                samples.append(cols[0])
                infos.append(cols[1] + ":" + cols[2] + ":" + cols[3])

    return samples, infos

def getSamples(bams, runOnly=-1):
    samples = []
    samplesInfos = []
    if len(bams) == 1 and (bams[0].endswith(".tsv") or bams[0].endswith(".csv")):
        # Sample file specified
        samples, samplesInfos = readSampleFile(bams[0])
    else:
        # List of BAM files specified
        samples = bams
        samplesInfos = [""] * len(samples)

    if(runOnly > 0):
        if(runOnly > len(samples)):
            raise RuntimeError("Sample index out of range. " + str(runOnly) + " > " + str(len(samples)) + ". Check -i/--sample-index")
        message("Running only job " + str(runOnly))
        samples = [ samples[runOnly - 1] ]
        samplesInfos = [ samplesInfos[runOnly - 1] ]
    elif(runOnly == 0):
        raise RuntimeError("Sample index (" + str(runOnly) + ") out of range. Starts with 1. Check -i/--sample-index")

    return samples, samplesInfos

def runMap(tid, inputBAM, referenceFile, threads, trim5p, maxPolyA, quantseqMapping, endtoendMapping, topn, sampleDescription, outputDirectory, skipSAM) :
    if skipSAM:
        outputSAM = os.path.join(outputDirectory, replaceExtension(basename(inputBAM), ".bam", "_slamdunk_mapped"))
    else:
        outputSAM = os.path.join(outputDirectory, replaceExtension(basename(inputBAM), ".sam", "_slamdunk_mapped"))
    outputLOG = os.path.join(outputDirectory, replaceExtension(basename(inputBAM), ".log", "_slamdunk_mapped"))

    #sampleName = "sample_" + str(tid)
    sampleName = replaceExtension(basename(inputBAM), ".bam", "")
    sampleType = "NA"
    sampleTime = "-1"
    if(sampleDescription != ""):
        sampleDescriptions = sampleDescription.split(":")
        if(len(sampleDescriptions) >= 1):
            sampleName = sampleDescriptions[0]
        if(len(sampleDescriptions) >= 2):
            typeDict = { 'p': 'pulse', 'c':'chase', 'pulse':'pulse', 'chase':'chase', '':'NA'}
            if sampleDescriptions[1] in typeDict:
                sampleType = typeDict[sampleDescriptions[1]]
            else:
                sampleType = sampleDescriptions[1]
        if(len(sampleDescriptions) >= 3):
            sampleTime = sampleDescriptions[2]

    mapper.Map(inputBAM, referenceFile, outputSAM, getLogFile(outputLOG), quantseqMapping, endtoendMapping, threads=threads, trim5p=trim5p, maxPolyA=maxPolyA, topn=topn, sampleId=tid, sampleName=sampleName, sampleType=sampleType, sampleTime=sampleTime, printOnly=printOnly, verbose=verbose)
    stepFinished()

def runSam2Bam(tid, bam, threads, outputDirectory):
    inputSAM = os.path.join(outputDirectory, replaceExtension(basename(bam), ".sam", "_slamdunk_mapped"))
    outputBAM = os.path.join(outputDirectory, replaceExtension(basename(bam), ".bam", "_slamdunk_mapped"))
    outputLOG = os.path.join(outputDirectory, replaceExtension(basename(bam), ".log", "_slamdunk_mapped"))
    mapper.sort(inputSAM, outputBAM, getLogFile(outputLOG), threads, False, printOnly, verbose)
    stepFinished()

def runDedup(tid, bam, outputDirectory) :
    outputBAM = os.path.join(outputDirectory, replaceExtension(basename(bam), ".bam", "_dedup"))
    outputLOG = os.path.join(outputDirectory, replaceExtension(basename(bam), ".log", "_dedup"))
    log = getLogFile(outputLOG)
    deduplicator.Dedup(bam, outputBAM, log)
    closeLogFile(log)
    stepFinished()

def runFilter(tid, bam, bed, mq, minIdentity, maxNM, outputDirectory):
    outputBAM = os.path.join(outputDirectory, replaceExtension(basename(bam), ".bam", "_filtered"))
    outputLOG = os.path.join(outputDirectory, replaceExtension(basename(bam), ".log", "_filtered"))
    filter.Filter(bam, outputBAM, getLogFile(outputLOG), bed, mq, minIdentity, maxNM, printOnly, verbose)
    stepFinished()

def runSnp(tid, referenceFile, minCov, minVarFreq, minQual, inputBAM, outputDirectory) :
    outputSNP = os.path.join(outputDirectory, replaceExtension(basename(inputBAM), ".vcf", "_snp"))
    outputLOG = os.path.join(outputDirectory, replaceExtension(basename(inputBAM), ".log", "_snp"))
    snps.SNPs(inputBAM, outputSNP, referenceFile, minVarFreq, minCov, minQual, getLogFile(outputLOG), printOnly, verbose, False)
    stepFinished()

def runCount(tid, bam, ref, bed, maxLength, minQual, conversionThreshold, outputDirectory, snpDirectory, vcfFile) :
    outputCSV = os.path.join(outputDirectory, replaceExtension(basename(bam), ".tsv", "_tcount"))
    outputBedgraphPlus = os.path.join(outputDirectory, replaceExtension(basename(bam), ".bedgraph", "_tcount_plus"))
    outputBedgraphMinus = os.path.join(outputDirectory, replaceExtension(basename(bam), ".bedgraph", "_tcount_mins"))
    outputLOG = os.path.join(outputDirectory, replaceExtension(basename(bam), ".log", "_tcount"))

    if (vcfFile != None) :
        inputSNP = vcfFile
    elif(snpDirectory != None):
        inputSNP = os.path.join(snpDirectory, replaceExtension(basename(bam), ".vcf", "_snp"))
    else:
        inputSNP = None

    if (maxLength == None) :
        maxLength = estimateMaxReadLength(bam)
    if (maxLength < 0) :
        print("Difference between minimum and maximum read length is > 10. Please specify --max-read-length parameter.")
        sys.exit(0)

    log = getLogFile(outputLOG)

    print("Using " + str(maxLength) + " as maximum read length.",file=log)

    tcounter.computeTconversions(ref, bed, inputSNP, bam, maxLength, minQual, outputCSV, outputBedgraphPlus, outputBedgraphMinus, conversionThreshold, log)
    stepFinished()
    return outputCSV

def runAll(args) :
    message("slamdunk all")

    if args.sampleIndex > -1:
        sec = random.randrange(200,2000) / 1000.0
        message("Waiting " + str(sec) + " seconds")
        sleep(sec)

    # Setup slamdunk run folder

    outputDirectory = args.outputDir
    createDir(outputDirectory)

    n = args.threads
    referenceFile = args.referenceFile

    # Run mapper dunk

    dunkPath = os.path.join(outputDirectory, "map")
    createDir(dunkPath)

    samples, samplesInfos = getSamples(args.files, runOnly=args.sampleIndex)

    message("Running slamDunk map for " + str(len(samples)) + " files (" + str(n) + " threads)")

    for i in range(0, len(samples)):
        bam = samples[i]

        if not args.sampleName or len(samples) > 1:
            sampleName = replaceExtension(basename(bam), "", "")
        else :
            sampleName = args.sampleName

        sampleInfo = samplesInfos[i]
        if sampleInfo == "":
            sampleInfo = sampleName + ":" + args.sampleType + ":" + str(args.sampleTime)
        tid = i
        if args.sampleIndex > -1:
            tid = args.sampleIndex
        runMap(tid, bam, referenceFile, n, args.trim5, args.maxPolyA, args.quantseq, args.endtoend, args.topn, sampleInfo, dunkPath, args.skipSAM)

    dunkFinished()

    if(not args.skipSAM):
        message("Running slamDunk sam2bam for " + str(len(samples)) + " files (" + str(n) + " threads)")
        results = Parallel(n_jobs=1, verbose=verbose)(delayed(runSam2Bam)(tid, samples[tid], n, dunkPath) for tid in range(0, len(samples)))
        dunkFinished()

    dunkbufferIn = []

    for file in samples :
        dunkbufferIn.append(os.path.join(dunkPath, replaceExtension(basename(file), ".bam", "_slamdunk_mapped")))

    # Run filter dunk

    bed = args.bed

    if args.filterbed:
        bed = args.filterbed
        args.multimap = True

    if (not args.multimap) :
        bed = None

    dunkPath = os.path.join(outputDirectory, "filter")
    createDir(dunkPath)

    message("Running slamDunk filter for " + str(len(samples)) + " files (" + str(n) + " threads)")
    results = Parallel(n_jobs=n, verbose=verbose)(delayed(runFilter)(tid, dunkbufferIn[tid], bed, args.mq, args.identity, args.nm, dunkPath) for tid in range(0, len(samples)))

    dunkFinished()

    # Run filter dunk

    dunkbufferOut = []

    for file in dunkbufferIn :
        dunkbufferOut.append(os.path.join(dunkPath, replaceExtension(basename(file), ".bam", "_filtered")))

    dunkbufferIn = dunkbufferOut

    dunkbufferOut = []

    dunkFinished()

    # Run snps dunk only if vcf not specified

    snpDirectory = None
    vcfFile = None

    if not "vcfFile" in args:

        dunkPath = os.path.join(outputDirectory, "snp")
        createDir(dunkPath)

        minCov = args.cov
        minVarFreq = args.var

        snpThread = n
        if(snpThread > 1):
            snpThread = int(snpThread / 2)

        #if (args.minQual == 0) :
        #    snpqual = 13
        #else :
        snpqual = args.minQual

        message("Running slamDunk SNP for " + str(len(samples)) + " files (" + str(snpThread) + " threads)")
        results = Parallel(n_jobs=snpThread, verbose=verbose)(delayed(runSnp)(tid, referenceFile, minCov, minVarFreq, snpqual, dunkbufferIn[tid], dunkPath) for tid in range(0, len(samples)))

        snpDirectory = os.path.join(outputDirectory, "snp")

        dunkFinished()
    else :
        vcfFile = args.vcfFile

    # Run count dunk

    dunkPath = os.path.join(outputDirectory, "count")
    createDir(dunkPath)

    message("Running slamDunk tcount for " + str(len(samples)) + " files (" + str(n) + " threads)")
    results = Parallel(n_jobs=n, verbose=verbose)(delayed(runCount)(tid, dunkbufferIn[tid], referenceFile, args.bed, args.maxLength, args.minQual, args.conversionThreshold, dunkPath, snpDirectory, vcfFile) for tid in range(0, len(samples)))

    dunkFinished()

def run():
    ########################################################################
    # Argument parsing
    ########################################################################

    # Info
    usage = "SLAMdunk software for analyzing SLAM-seq data"

    # Main Parsers
    parser = ArgumentParser(description=usage, formatter_class=ArgumentDefaultsHelpFormatter)
    parser.add_argument('--version', action='version', version='%(prog)s ' + __version__)

    # Initialize Subparsers
    subparsers = parser.add_subparsers(help="", dest="command")

    # map command

    mapparser = subparsers.add_parser('map', help='Map SLAM-seq read data', formatter_class=ArgumentDefaultsHelpFormatter)
    mapparser.add_argument('files', action='store', help='Single csv/tsv file (recommended) containing all sample files and sample info or a list of all sample BAM/FASTA(gz)/FASTQ(gz) files' , nargs="+")
    mapparser.add_argument("-r", "--reference", type=str, required=True, dest="referenceFile", default=SUPPRESS, help="Reference fasta file")
    mapparser.add_argument("-o", "--outputDir", type=str, required=True, dest="outputDir", default=SUPPRESS, help="Output directory for mapped BAM files.")
    mapparser.add_argument("-5", "--trim-5p", type=int, required=False, dest="trim5", default=12, help="Number of bp removed from 5' end of all reads.")
    mapparser.add_argument("-n", "--topn", type=int, required=False, dest="topn", default=1, help="Max. number of alignments to report per read")
    mapparser.add_argument("-a", "--max-polya", type=int, required=False, dest="maxPolyA", default=4, help="Max number of As at the 3' end of a read.")
    mapparser.add_argument("-t", "--threads", type=int, required=False, dest="threads", default = 1, help="Thread number")
    mapparser.add_argument("-q", "--quantseq", dest="quantseq", action='store_true', required=False, help="Run plain Quantseq alignment without SLAM-seq scoring")
    mapparser.add_argument('-e', "--endtoend", action='store_true', dest="endtoend", help="Use a end to end alignment algorithm for mapping.")
    mapparser.add_argument('-sn', "--sampleName", type=str, dest="sampleName", required = False, help="Use this sample name for all supplied samples")
    mapparser.add_argument('-sy', "--sampleType", type=str, dest="sampleType", required = False, default = "pulse", help="Use this sample type for all supplied samples")
    mapparser.add_argument('-st', "--sampleTime", type=int, dest="sampleTime", required = False, default = 0, help="Use this sample time for all supplied samples")
    mapparser.add_argument("-i", "--sample-index", type=int, required=False, default=-1, dest="sampleIndex", help="Run analysis only for sample <i>. Use for distributing slamdunk analysis on a cluster (index is 1-based).")
    mapparser.add_argument('-ss', "--skip-sam", action='store_true', dest="skipSAM", help="Output BAM while mapping. Slower but, uses less hard disk.")

    # filter command

    filterparser = subparsers.add_parser('filter', help='Filter SLAM-seq aligned data')
    filterparser.add_argument('bam', action='store', help='Bam file(s)' , nargs="+")
    filterparser.add_argument("-o", "--outputDir", type=str, required=True, dest="outputDir", help="Output directory for mapped BAM files.")
    filterparser.add_argument("-b", "--bed", type=str, required=False, dest="bed", help="BED file, overrides MQ filter to 0")
    filterparser.add_argument("-mq", "--min-mq", type=int, required=False, default=2, dest="mq", help="Minimum mapping quality (default: %(default)d)")
    filterparser.add_argument("-mi", "--min-identity", type=float, required=False, default=0.95, dest="identity", help="Minimum alignment identity (default: %(default)s)")
    filterparser.add_argument("-nm", "--max-nm", type=int, required=False, default=-1, dest="nm", help="Maximum NM for alignments (default: %(default)d)")
    filterparser.add_argument("-t", "--threads", type=int, required=False, dest="threads", default=1, help="Thread number (default: %(default)d)")

    # snp command

    snpparser = subparsers.add_parser('snp', help='Call SNPs on SLAM-seq aligned data', formatter_class=ArgumentDefaultsHelpFormatter)
    snpparser.add_argument('bam', action='store', help='Bam file(s)' , nargs="+")
    snpparser.add_argument("-o", "--outputDir", type=str, required=True, dest="outputDir", default=SUPPRESS, help="Output directory for mapped BAM files.")
    snpparser.add_argument("-r", "--reference", required=True, dest="fasta", type=str, default=SUPPRESS, help="Reference fasta file")
    snpparser.add_argument("-c", "--min-coverage", required=False, dest="cov", type=int, help="Minimimum coverage to call variant", default=10)
    #snpparser.add_argument("-q", "--min-base-qual", type=int, default=13, required=False, dest="minQual", help="Min base quality for T -> C conversions (default: %(default)d)")
    snpparser.add_argument("-f", "--var-fraction", required=False, dest="var", type=float, help="Minimimum variant fraction to call variant", default=0.8)
    snpparser.add_argument("-t", "--threads", type=int, required=False, default=1, dest="threads", help="Thread number")

    # count command

    countparser = subparsers.add_parser('count', help='Count T/C conversions in SLAM-seq aligned data')
    countparser.add_argument('bam', action='store', help='Bam file(s)' , nargs="+")
    countparser.add_argument("-o", "--outputDir", type=str, required=True, dest="outputDir", default=SUPPRESS, help="Output directory for mapped BAM files.")
    countparser.add_argument("-s", "--snp-directory", type=str, required=False, dest="snpDir", default=SUPPRESS, help="Directory containing SNP files.")
    countparser.add_argument("-v", "--vcf", type=str, required=False, dest="vcfFile", default=SUPPRESS, help="Externally provided custom variant file.")
    countparser.add_argument("-r", "--reference", type=str, required=True, dest="ref", default=SUPPRESS, help="Reference fasta file")
    countparser.add_argument("-b", "--bed", type=str, required=True, dest="bed", default=SUPPRESS, help="BED file")
    countparser.add_argument("-c", "--conversion-threshold", type=int, dest="conversionThreshold", required=False, default=1,help="Number of T>C conversions required to count read as T>C read (default: %(default)d)")
    countparser.add_argument("-l", "--max-read-length", type=int, required=False, dest="maxLength", help="Max read length in BAM file")
    countparser.add_argument("-q", "--min-base-qual", type=int, default=27, required=False, dest="minQual", help="Min base quality for T -> C conversions (default: %(default)d)")
    countparser.add_argument("-t", "--threads", type=int, required=False, default=1, dest="threads", help="Thread number (default: %(default)d)")


    # all command

    allparser = subparsers.add_parser('all', help='Run entire SLAMdunk analysis')
    allparser.add_argument('files', action='store', help='Single csv/tsv file (recommended) containing all sample files and sample info or a list of all sample BAM/FASTA(gz)/FASTQ(gz) files' , nargs="+")
    allparser.add_argument("-r", "--reference", type=str, required=True, dest="referenceFile", help="Reference fasta file")
    allparser.add_argument("-b", "--bed", type=str, required=True, dest="bed", help="BED file with 3'UTR coordinates")
    allparser.add_argument("-fb", "--filterbed", type=str, required=False, dest="filterbed", help="BED file with 3'UTR coordinates to filter multimappers (activates -m)")
    allparser.add_argument("-v", "--vcf", type=str, required=False, dest="vcfFile", default=SUPPRESS, help="Skip SNP step and provide custom variant file.")
    allparser.add_argument("-o", "--outputDir", type=str, required=True, dest="outputDir", help="Output directory for slamdunk run.")
    allparser.add_argument("-5", "--trim-5p", type=int, required=False, dest="trim5", default=12, help="Number of bp removed from 5' end of all reads (default: %(default)s)")
    allparser.add_argument("-a", "--max-polya", type=int, required=False, dest="maxPolyA", default=4, help="Max number of As at the 3' end of a read (default: %(default)s)")
    allparser.add_argument("-n", "--topn", type=int, required=False, dest="topn", default=1, help="Max. number of alignments to report per read (default: %(default)s)")
    allparser.add_argument("-t", "--threads", type=int, required=False, default=1, dest="threads", help="Thread number (default: %(default)s)")
    allparser.add_argument("-q", "--quantseq", dest="quantseq", action='store_true', required=False, help="Run plain Quantseq alignment without SLAM-seq scoring")
    allparser.add_argument('-e', "--endtoend", action='store_true', dest="endtoend", help="Use a end to end alignment algorithm for mapping.")
    allparser.add_argument('-m', "--multimap", action='store_true', dest="multimap", help="Use reference to resolve multimappers (requires -n > 1).")
    allparser.add_argument("-mq", "--min-mq", type=int, required=False, default=2, dest="mq", help="Minimum mapping quality (default: %(default)s)")
    allparser.add_argument("-mi", "--min-identity", type=float, required=False, default=0.95, dest="identity", help="Minimum alignment identity (default: %(default)s)")
    allparser.add_argument("-nm", "--max-nm", type=int, required=False, default=-1, dest="nm", help="Maximum NM for alignments (default: %(default)s)")
    allparser.add_argument("-mc", "--min-coverage", required=False, dest="cov", type=int, help="Minimimum coverage to call variant (default: %(default)s)", default=10)
    allparser.add_argument("-mv", "--var-fraction", required=False, dest="var", type=float, help="Minimimum variant fraction to call variant (default: %(default)s)", default=0.8)
    allparser.add_argument("-c", "--conversion-threshold", type=int, dest="conversionThreshold", required=False, default=1,help="Number of T>C conversions required to count read as T>C read (default: %(default)d)")
    allparser.add_argument("-rl", "--max-read-length", type=int, required=False, dest="maxLength", help="Max read length in BAM file")
    allparser.add_argument("-mbq", "--min-base-qual", type=int, default=27, required=False, dest="minQual", help="Min base quality for T -> C conversions (default: %(default)d)")
    allparser.add_argument('-sn', "--sampleName", type=str, dest="sampleName", required = False, help="Use this sample name for all supplied samples")
    allparser.add_argument('-sy', "--sampleType", type=str, dest="sampleType", required = False, default = "pulse", help="Use this sample type for all supplied samples")
    allparser.add_argument('-st', "--sampleTime", type=int, dest="sampleTime", required = False, default = 0, help="Use this sample time for all supplied samples")
    allparser.add_argument("-i", "--sample-index", type=int, required=False, default=-1, dest="sampleIndex", help="Run analysis only for sample <i>. Use for distributing slamdunk analysis on a cluster (index is 1-based).")
    allparser.add_argument("-ss", "--skip-sam", action='store_true', dest="skipSAM", help="Output BAM while mapping. Slower but, uses less hard disk.")

    args = parser.parse_args()

    ########################################################################
    # Routine selection
    ########################################################################

    command = args.command

    if (command == "map") :
        mapper.checkNextGenMapVersion()

        outputDirectory = args.outputDir

        if args.sampleIndex > -1:
            sec = random.randrange(0,2000) / 1000
            message("Waiting " + str(sec) + " seconds")
            sleep(sec)

        createDir(outputDirectory)
        n = args.threads
        referenceFile = args.referenceFile

        samples, samplesInfos = getSamples(args.files, runOnly=args.sampleIndex)

        message("Running slamDunk map for " + str(len(samples)) + " files (" + str(n) + " threads)")
        for i in range(0, len(samples)):
            bam = samples[i]

            if not args.sampleName or len(samples) > 1:
                sampleName = replaceExtension(basename(bam), "", "")
            else :
                sampleName = args.sampleName

            sampleInfo = samplesInfos[i]
            if sampleInfo == "":
                sampleInfo = sampleName + ":" + args.sampleType + ":" + str(args.sampleTime)
            tid = i
            if args.sampleIndex > -1:
                tid = args.sampleIndex
            runMap(tid, bam, referenceFile, n, args.trim5, args.maxPolyA, args.quantseq, args.endtoend, args.topn, sampleInfo, outputDirectory, args.skipSAM)

        dunkFinished()

        if not args.skipSAM:
            message("Running slamDunk sam2bam for " + str(len(samples)) + " files (" + str(n) + " threads)")
            results = Parallel(n_jobs=1, verbose=verbose)(delayed(runSam2Bam)(tid, samples[tid], n, outputDirectory) for tid in range(0, len(samples)))
            dunkFinished()

    elif (command == "filter") :
        outputDirectory = args.outputDir
        createDir(outputDirectory)
        n = args.threads
        message("Running slamDunk filter for " + str(len(args.bam)) + " files (" + str(n) + " threads)")
        results = Parallel(n_jobs=n, verbose=verbose)(delayed(runFilter)(tid, args.bam[tid], args.bed, args.mq, args.identity, args.nm, outputDirectory) for tid in range(0, len(args.bam)))
        dunkFinished()

    elif (command == "snp") :
        outputDirectory = args.outputDir
        createDir(outputDirectory)
        fasta = args.fasta
        minCov = args.cov
        minVarFreq = args.var
        #minQual = args.minQual
        minQual = 15
        n = args.threads
        if(n > 1):
            n = int(n / 2)
        message("Running slamDunk SNP for " + str(len(args.bam)) + " files (" + str(n) + " threads)")
        results = Parallel(n_jobs=n, verbose=verbose)(delayed(runSnp)(tid, fasta, minCov, minVarFreq, minQual, args.bam[tid], outputDirectory) for tid in range(0, len(args.bam)))
        dunkFinished()

    elif (command == "count") :
        outputDirectory = args.outputDir
        createDir(outputDirectory)
        if "snpDir" in args:
            snpDirectory = args.snpDir
        else :
            snpDirectory = None
        if "vcfFile" in args:
            vcfFile = args.vcfFile
        else :
            vcfFile = None
        n = args.threads
        message("Running slamDunk tcount for " + str(len(args.bam)) + " files (" + str(n) + " threads)")
        results = Parallel(n_jobs=n, verbose=verbose)(delayed(runCount)(tid, args.bam[tid], args.ref, args.bed, args.maxLength, args.minQual, args.conversionThreshold, outputDirectory, snpDirectory, vcfFile) for tid in range(0, len(args.bam)))
        dunkFinished()

    elif (command == "all") :
        runAll(args)
        dunkFinished()

    else:
        parser.error("Too few arguments.")

if __name__ == '__main__':
    run()
