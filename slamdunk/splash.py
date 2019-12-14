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

#SPLASH:

#Systematic Performance evaLuAtion of Slamdunk beHaviour (accomplishemnt, achievement)


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

from slamdunk.dunks import simulator
from slamdunk.utils.misc import replaceExtension, removeExtension, SampleInfo
from slamdunk.version import __version__
from shutil import copyfile

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

def message(msg):
    print(msg, file=mainOutput)

def createDir(directory):
    if directory:
        if not os.path.exists(directory):
            message("Creating output directory: " + directory)
            os.makedirs(directory)

def reads(outputDirectory, bed, sampleName, readLenght, readNumber, readCoverage, seqError, pulseTimePoint, chaseTimePoint, conversionRate, sampleInfo, labledTranscripots = -1.0):
    message("Simulating read sample: " + sampleName)

    bed12File = replaceExtension(bed, ".bed12")
    bed12FastaFile = replaceExtension(bed, ".fa")
    explvFile = replaceExtension(bed, ".eplv")

    bedReads = os.path.join(outputDirectory, sampleName + "_reads_tmp.bed")
    faReads = os.path.join(outputDirectory, sampleName + "_reads_tmp.fa")

    totalUTRlength = simulator.getTotalUtrLength(bed12File)

    if(readNumber == 0):
        readNumber = (totalUTRlength / readLenght) *  readCoverage
        readNumber = int(readNumber * (random.uniform(-0.2, 0.2) + 1))

    #message("Simulating " + str(readNumber) + " reads with sequencing error of " + str(seqError))
    simulator.simulateReads(bed12File, bed12FastaFile, explvFile, bedReads, faReads, readLenght, readNumber, seqError)

    bamReadsWithTC = os.path.join(outputDirectory, sampleName + "_reads.bam")
    utrSummary = os.path.join(outputDirectory, sampleName + "_utrsummary.tsv")

    simulator.addTcConversions(bed, faReads, bamReadsWithTC, pulseTimePoint, chaseTimePoint, utrSummary, conversionRate, readNumber, sampleInfo, labledTranscripots)

    os.unlink(faReads)
    os.unlink(bedReads)


def run():
    ########################################################################
    # Argument parsing
    ########################################################################

    # TODO: parameter for simulating expression levels
    # TODO: more realistic simulation of half lifes

    # Info
    usage = "SLAMdunk software for simulating SLAM-seq data"

    # Main Parsers
    parser = ArgumentParser(description=usage, formatter_class=RawDescriptionHelpFormatter)
    parser.add_argument('--version', action='version', version='%(prog)s ' + __version__)

    # Initialize Subparsers
    subparsers = parser.add_subparsers(help="", dest="command")

    allparse = subparsers.add_parser('all', help='Simulated full SlamSeq samples')
    allparse.add_argument("-r", "--reference", type=str, required=True, dest="referenceFile", help="Reference fasta file")
    allparse.add_argument("-b", "--bed", type=str, required=True, dest="bed", help="BED file")
    allparse.add_argument("-l", "--read-length", type=int, required=True, dest="readLength", help="All UTRs short than the read length are removed.")
    allparse.add_argument("-o", "--outputDir", type=str, required=False, dest="outputDir", default=".", help="Output directory for mapped BAM files.")
    allparse.add_argument("-s", "--snp-rate", type=float, required=False, default=0.001, dest="snpRate", help="SNP rate in UTRs")
    allparse.add_argument("-cov", "--read-coverage", type=int, required=False, default=20, dest="readCoverage", help="Read coverage (if read number is not specified)")
    allparse.add_argument("-e", "--sequencing-error", type=float, required=False, default=0.05, dest="seqError", help="Sequencing error")
    allparse.add_argument("-p", "--pulse", type=str, required=False, dest="pulse", help="Pulse in minutes")
    allparse.add_argument("-ra", "--rates", type=str, required=False, default=None, dest="rates", help="List of rates")
    allparse.add_argument("-c", "--chase", type=str, required=False, default="", dest="chase", help="Chase in minutes")
    allparse.add_argument("-tc", "--tc-rate", type=float, required=False, dest="conversionRate", default=0.024, help="T->C conversion rate")
    allparse.add_argument("-minhl", "--min-halflife", type=int, required=False, default=30, dest="minHalfLife", help="Lower bound for the simulated half lifes in minutes")
    allparse.add_argument("-maxhl", "--max-halflife", type=int, required=False, default=720, dest="maxHalfLife", help="Upper bound for the simulated half lifes in minutes")
    allparse.add_argument("-t", "--threads", type=int, required=False, default=1, dest="threads", help="Thread number")
    allparse.add_argument("-rep", "--replicates", type=int, required=False, default=1, dest="replicates", help="Number of replicates")
    allparse.add_argument('-st', "--skip-turnover", required=False, dest="skipTurnover", action='store_true', help="Take half-life from score filed of input BED file")

    preparebedparse = subparsers.add_parser('preparebed', help='Prepares a UTR BED file for SlamSim')
    preparebedparse.add_argument("-b", "--bed", type=str, required=True, dest="bed", help="BED file")
    preparebedparse.add_argument("-l", "--read-length", type=int, required=True, dest="readLength", help="All UTRs short than the read length are removed.")
    preparebedparse.add_argument("-o", "--outputDir", type=str, required=False, dest="outputDir", default=".", help="Output directory for mapped BAM files.")

    turnoverparse = subparsers.add_parser('turnover', help='Simulate utrs and turnover rate')
    turnoverparse.add_argument("-b", "--bed", type=str, required=True, dest="bed", help="BED file")
    turnoverparse.add_argument("-minhl", "--min-halflife", type=int, required=False, default=30, dest="minHalfLife", help="Lower bound for the simulated half lifes in minutes")
    turnoverparse.add_argument("-maxhl", "--max-halflife", type=int, required=False, default=720, dest="maxHalfLife", help="Upper bound for the simulated half lifes in minutes")
    turnoverparse.add_argument("-o", "--outputDir", type=str, required=False, dest="outputDir", default=".", help="Output directory for mapped BAM files.")

    utrsparse = subparsers.add_parser('utrs', help='Simulate utrs and turnover rate')
    utrsparse.add_argument("-r", "--reference", type=str, required=True, dest="referenceFile", help="Reference fasta file")
    utrsparse.add_argument("-b", "--bed", type=str, required=True, dest="bed", help="BED file")
    utrsparse.add_argument("-l", "--read-length", type=int, required=True, dest="readLength", help="Read length")
    utrsparse.add_argument("-o", "--outputDir", type=str, required=False, dest="outputDir", default=".", help="Output directory for mapped BAM files.")
    utrsparse.add_argument("-s", "--snp-rate", type=float, required=False, default=0.001, dest="snpRate", help="SNP rate in UTRs")


    simulateparse = subparsers.add_parser('reads', help='Simulate SLAM-seq read data')
    simulateparse.add_argument("-o", "--outputDir", type=str, required=False, dest="outputDir", default=".", help="Output directory for mapped BAM files.")
    simulateparse.add_argument("--sample-name", type=str, required=True, dest="sampleName", help="Name of sample")
    simulateparse.add_argument("-b", "--bed", type=str, required=True, dest="bed", help="BED file")
    simulateparse.add_argument("-l", "--read-length", type=int, required=True, dest="readLength", help="Read length")
    simulateparse.add_argument("-n", "--read-number", type=int, required=False, default=0, dest="readNumber", help="Number of reads to simulate")
    simulateparse.add_argument("-cov", "--read-coverage", type=int, required=False, default=20, dest="readCoverage", help="Read coverage (if read number is not specified)")
    simulateparse.add_argument("-e", "--sequencing-error", type=float, required=False, default=0.05, dest="seqError", help="Sequencing error")
    simulateparse.add_argument("-p", "--pulse", type=int, required=True, dest="pulse", help="Pulse in minutes")
    simulateparse.add_argument("-c", "--chase", type=int, required=False, default=0, dest="chase", help="Chase in minutes")
    simulateparse.add_argument("-tc", "--tc-rate", type=float, required=False, dest="conversionRate", default=0.024, help="T->C conversion rate")

    evalparser = subparsers.add_parser('eval-counts', help='Evaluate count files')
    evalparser.add_argument("-s", "--simulated", type=str, required=True, dest="simulated", help="")
    evalparser.add_argument("-d", "--slamdun", type=str, required=True, dest="slamdunk", help="")
    evalparser.add_argument("-o", "--outputFile", type=str, required=True, dest="outputFile", help="")

    evalreadsparser = subparsers.add_parser('eval-reads', help='Evaluate read files')
    evalreadsparser.add_argument("-o", "--outputFile", type=str, required=True, dest="outputFile", help="")
    evalreadsparser.add_argument("-b", "--bed", type=str, required=True, dest="bed", help="BED file")
    evalreadsparser.add_argument("-r", "--reference", type=str, required=True, dest="referenceFile", help="Reference fasta file")
    evalreadsparser.add_argument('bam', action='store', help='Bam file(s)' , nargs="+")

    evalconversionplotparse = subparsers.add_parser('plot.conversions', help='Plots differences in simulated and found conversion rates')
    evalconversionplotparse.add_argument("-sim", "--simDir", type=str, required=True, dest="simDir", help="")
    evalconversionplotparse.add_argument("-slam", "--slamdunkDir", type=str, required=True, dest="slamDir", help="")
    evalconversionplotparse.add_argument("-o", "--outputFile", type=str, required=True, dest="outputFile", help="")
    evalconversionplotparse.add_argument("-tc", "--tc-rate", type=float, required=False, dest="conversionRate", default=0.03, help="T->C conversion rate")

    evalhalflifeparse = subparsers.add_parser('plot.halflifes', help='Plots half lifes')
    evalhalflifeparse.add_argument("-sim", "--simulated-hl", type=str, required=True, dest="simHL", help="Simulated half-lifes")
    evalhalflifeparse.add_argument("-pred", "--predicted-hl", type=str, required=True, dest="predHL", help="Predicted half-lifes")
    evalhalflifeparse.add_argument("-true", "--true-hl", type=str, required=True, dest="trueHL", help="Predicted half-lifes")
    evalhalflifeparse.add_argument("-o", "--outputFile", type=str, required=True, dest="outputFile", help="")
    evalhalflifeparse.add_argument("-e", "--erroroutputFile", type=str, required=True, dest="erroutputFile", help="")

    evalhalflifeplotparse = subparsers.add_parser('plot.halflifespergene', help='Plots half lifes')
    evalhalflifeplotparse.add_argument("-sim", "--simDir", type=str, required=True, dest="simDir", help="")
    evalhalflifeplotparse.add_argument("-slam", "--slamdunkDir", type=str, required=True, dest="slamDir", help="")
    evalhalflifeplotparse.add_argument("-t", "--timepoints", type=str, required=True, dest="timepoints", help="")
    evalhalflifeplotparse.add_argument("-o", "--outputFile", type=str, required=True, dest="outputFile", help="")
    evalhalflifeplotparse.add_argument("-tc", "--tc-rate", type=float, required=False, dest="conversionRate", default=0.03, help="T->C conversion rate")
    evalhalflifeplotparse.add_argument("-b", "--bed", type=str, required=True, dest="bed", help="BED file")

    utilcrateparse = subparsers.add_parser('util.conversionrate', help='Get conversion rate from mapped BAM files')
    utilcrateparse.add_argument('bam', action='store', help='Bam file(s)' , nargs="+")
    utilcrateparse.add_argument("-r", "--reference", type=str, required=True, dest="referenceFile", help="Reference fasta file")
    utilcrateparse.add_argument("-region", "--region", type=str, required=True, dest="region", help="")
    utilcrateparse.add_argument('-rev',required=False, dest="reverse", action='store_true')

    args = parser.parse_args()

    ########################################################################
    # Routine selection
    ########################################################################

    def prepareBed(outputDirectory, bed, readLength):

        createDir(outputDirectory)
        slamSimBed = os.path.join(outputDirectory, replaceExtension(basename(bed), ".bed", "_original"))
        simulator.prepareBED(bed, slamSimBed, readLength)

    def turnOver(outputDirectory, bed, minHalflife, maxHalfLife, skipTurnover=False):
        message("Simulating turnover")
        createDir(outputDirectory)
        trunoverBed = os.path.join(outputDirectory, replaceExtension(basename(bed), ".bed", "_utrs"))
        if not skipTurnover:
            simulator.simulateTurnOver(bed, trunoverBed, minHalflife, maxHalfLife)
        else:
            copyfile(bed, trunoverBed)

    def Utrs(outputDirectory, bed, referenceFasta, readLength, polyALength, snpRate):
        message("Simulating UTRs")
        createDir(outputDirectory)
        bed12 = os.path.join(outputDirectory, replaceExtension(basename(bed), ".bed12", "_utrs"))
        bed12Fasta = os.path.join(outputDirectory, replaceExtension(basename(bed), ".fa", "_utrs"))
        explv = os.path.join(outputDirectory, replaceExtension(basename(bed), ".eplv", "_utrs"))
        vcfFile = os.path.join(outputDirectory, replaceExtension(basename(bed), ".vcf", "_utrs"))

        totalUTRlength = simulator.prepareUTRs(bed, bed12, bed12Fasta, referenceFasta, readLength, polyALength, explv, snpRate, vcfFile)

    command = args.command
    if (command == "preparebed") :
        prepareBed(args.outputDir, args.bed, args.readLength)

    elif (command == "turnover"):
        turnOver(args.outputDir, args.bed, args.minHalfLife, args.maxHalfLife)

    elif (command == "utrs") :
        polyALength = 0
        Utrs(args.outputDir, args.bed, args.referenceFile, args.readLength, polyALength, args.snpRate)

    elif (command == "reads") :
        createDir(args.outputDir)
        reads(args.outputDir, args.bed, args.sampleName, args.readLength, args.readNumber, args.readCoverage, args.seqError, args.pulse, args.chase, args.conversionRate)
    elif (command == "eval-counts") :
        outputPath = os.path.dirname(args.outputFile)
        createDir(outputPath)
        simulator.evaluate(args.simulated, args.slamdunk, args.outputFile, mainOutput)
    elif (command == "eval-reads") :
        outputPath = os.path.dirname(args.outputFile)
        createDir(outputPath)
        for bam in args.bam:
            simulator.evaluateReads(bam, args.referenceFile, args.bed, args.outputFile, mainOutput)
    elif (command == "plot.conversions") :

        simDir = args.simDir
        slamDir = args.slamDir
        outputPDF = args.outputFile
        conversionRate = args.conversionRate

        outputPath = os.path.dirname(outputPDF)
        createDir(outputPath)

        simulator.plotconversiondifferences(simDir, slamDir, conversionRate, outputPDF)

    elif (command == "plot.halflifespergene") :

        bed = args.bed
        simDir = args.simDir
        slamDir = args.slamDir
        outputPDF = args.outputFile
        conversionRate = args.conversionRate
        timePoints = args.timepoints
        outputPath = os.path.dirname(outputPDF)
        createDir(outputPath)

        simulator.plotHalfLifes(bed, simDir, slamDir, timePoints, conversionRate, outputPDF)

    elif (command == "plot.halflifes") :

        trueHLFile = args.trueHL
        simHLFile = args.simHL
        predHLFile = args.predHL


        outputPDF = args.outputFile
        erroutputCSV = args.erroutputFile

        simulator.evalHalfLifes(trueHLFile, simHLFile, predHLFile, outputPDF, erroutputCSV)

    elif (command == "util.conversionrate") :

        ref = args.referenceFile
        bams = args.bam
        region = args.region
        region = region.replace(",", "")
        chromosome = region.split(":")[0]
        start = int(region.split(":")[1].split("-")[0])
        end = int(region.split(":")[1].split("-")[1])
        strand = "+"
        if(args.reverse):
            strand = "-"
        for bam in bams:
            simulator.getConversionRateFromBam(bam, ref, chromosome, start, end, strand)

    elif (command == "all") :

        #args.outputDir, args.bed, args.sampleName, args.readLength, args.readNumber, args.readCoverage, args.seqError, args.pulse, args.chase, args.conversionRate

        referenceFile = args.referenceFile

        baseFolder = args.outputDir
        annotationFile = args.bed

        readLength = args.readLength
        readCoverage = args.readCoverage
        sequencingError = args.seqError
        polyALength = 0

        #timePoints = [0, 15, 30, 60, 180, 360, 720, 1440]
        if not args.pulse == None:
            timePoints = args.pulse.split(",")
            chaseTimePoints = []
        if len(args.chase) > 0:
            chaseTimePoints = args.chase.split(",")

        labledTranscripots = None
        if not args.rates == None:
            labledTranscripots = args.rates.split(",")

        replicates = args.replicates

        n = args.threads

        createDir(baseFolder)

        annotationPrefix = removeExtension(basename(annotationFile))
        simulatedAnnotationPref = os.path.join(baseFolder, annotationPrefix)

        prepareBed(baseFolder, annotationFile, readLength)

        # TODO parameter to skip this
        turnOver(baseFolder, simulatedAnnotationPref + "_original.bed", args.minHalfLife, args.maxHalfLife, args.skipTurnover)

        Utrs(baseFolder, simulatedAnnotationPref + "_original.bed", referenceFile, readLength, polyALength, args.snpRate)

        sampleFile = open(os.path.join(baseFolder, "samples.tsv"), "w")

        sampleNumber = 1
        jobs = []

        if(labledTranscripots == None):
            for timePoint in timePoints:
                for replicate in range(1, replicates + 1):
                    sampleName =  "sample_" + str(sampleNumber) + "_pulse_" + str(timePoint) + "min_rep" + str(replicate)
                    sampleInfo = SampleInfo(ID = sampleNumber, Name = sampleName, Type = "pulse", Time = str(timePoint))

                    jobs.append(delayed(reads)(baseFolder,
                                simulatedAnnotationPref + "_original_utrs.bed",
                                sampleName,
                                readLength, 0, readCoverage, sequencingError,
                                int(timePoint), 0, args.conversionRate, sampleInfo))

                    sampleNumber += 1
                    print(os.path.join(baseFolder, sampleName + "_reads.bam"), sampleName, "pulse", timePoint, sep="\t", file=sampleFile)

            for timePoint in chaseTimePoints:
                for replicate in range(1, replicates + 1):
                    sampleName =  "sample_" + str(sampleNumber) + "_chase_" + str(timePoint) + "min_rep" + str(replicate)
                    sampleInfo = SampleInfo(ID = sampleNumber, Name = sampleName, Type = "chase", Time = str(timePoint))

                    jobs.append(delayed(reads)(baseFolder,
                                simulatedAnnotationPref + "_original_utrs.bed",
                                sampleName,
                                readLength, 0, readCoverage, sequencingError,
                                int(timePoints[-1]), int(timePoint), args.conversionRate, sampleInfo))

                    sampleNumber += 1
                    print(os.path.join(baseFolder, sampleName + "_reads.bam"), sampleName, "chase", timePoint, sep="\t", file=sampleFile)
        else:
            for rate in labledTranscripots:
                for replicate in range(1, replicates + 1):
                    sampleName =  "sample_" + str(sampleNumber) + "_rate_" + str(rate) + "_rep" + str(replicate)
                    sampleInfo = SampleInfo(ID = sampleNumber, Name = sampleName, Type = "rate", Time = str(rate))

                    jobs.append(delayed(reads)(baseFolder,
                                simulatedAnnotationPref + "_original_utrs.bed",
                                sampleName,
                                readLength, 0, readCoverage, sequencingError,
                                0, 0, args.conversionRate, sampleInfo, float(rate)))

                    sampleNumber += 1
                    print(os.path.join(baseFolder, sampleName + "_reads.bam"), sampleName, "rate", rate, sep="\t", file=sampleFile)


        sampleFile.close()

        results = Parallel(n_jobs=n, verbose=False)(jobs)

    else:
        parser.error("Too few arguments.")

if __name__ == '__main__':
    run()
