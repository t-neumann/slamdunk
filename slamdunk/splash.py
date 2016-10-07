#!/usr/bin/env python

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

from slamdunk.dunks import simulator
from slamdunk.utils.misc import replaceExtension

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

def run():
    ########################################################################
    # Argument parsing
    ########################################################################
    
    # TODO: parameter for simulating expression levels
    # TODO: more realistic simulation of half lifes
    
    # Info
    usage = "SLAMdunk software for simulating SLAM-seq data"
    version = "1.0"
    
    # Main Parsers
    parser = ArgumentParser(description=usage, formatter_class=RawDescriptionHelpFormatter, version=version)
    
    # Initialize Subparsers
    subparsers = parser.add_subparsers(help="", dest="command")
    
    preparebedparse = subparsers.add_parser('preparebed', help='Prepares a UTR BED file for SlamSim')
    preparebedparse.add_argument("-b", "--bed", type=str, required=True, dest="bed", help="BED file")
    preparebedparse.add_argument("-l", "--read-length", type=int, required=True, dest="readLength", help="All UTRs short than the read length are removed.")
    preparebedparse.add_argument("-o", "--outputDir", type=str, required=False, dest="outputDir", default=".", help="Output directory for mapped BAM files.")
    
    turnoverparse = subparsers.add_parser('utrs', help='Simulate utrs and turnover rate')
    turnoverparse.add_argument("-r", "--reference", type=str, required=True, dest="referenceFile", help="Reference fasta file")
    turnoverparse.add_argument("-b", "--bed", type=str, required=True, dest="bed", help="BED file")
    turnoverparse.add_argument("-minhl", "--min-halflife", type=int, required=False, default=30, dest="minHalfLife", help="Lower bound for the simulated half lifes in minutes")
    turnoverparse.add_argument("-maxhl", "--max-halflife", type=int, required=False, default=720, dest="maxHalfLife", help="Upper bound for the simulated half lifes in minutes")
    turnoverparse.add_argument("-o", "--outputDir", type=str, required=False, dest="outputDir", default=".", help="Output directory for mapped BAM files.")
    turnoverparse.add_argument("-s", "--snp-rate", type=float, required=False, default=0.001, dest="snpRate", help="SNP rate in UTRs")
    
    simulateparse = subparsers.add_parser('reads', help='Simulate SLAM-seq read data')
    simulateparse.add_argument("-o", "--outputDir", type=str, required=False, dest="outputDir", default=".", help="Output directory for mapped BAM files.")
    simulateparse.add_argument("--sample-name", type=str, required=True, dest="sampleName", help="Name of sample")
    simulateparse.add_argument("-b", "--bed", type=str, required=True, dest="bed", help="BED file")
    simulateparse.add_argument("-l", "--read-length", type=int, required=True, dest="readLength", help="Read length")
    simulateparse.add_argument("-n", "--read-number", type=int, required=False, default=0, dest="readNumber", help="Number of reads to simulate")
    simulateparse.add_argument("-cov", "--read-coverage", type=int, required=False, default=20, dest="readCoverage", help="Read coverage (if read number is not specified)")
    simulateparse.add_argument("-e", "--sequencing-error", type=float, required=False, default=0.05, dest="seqError", help="Sequencing error")
    simulateparse.add_argument("-t", "--timepoint", type=int, required=True, dest="timePoint", help="Timepoint in minutes")
    simulateparse.add_argument("-tc", "--tc-rate", type=float, required=False, dest="conversionRate", default=0.03, help="T->C conversion rate")
        
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
    
    command = args.command
    if (command == "preparebed") :
        outputDirectory = args.outputDir
        createDir(outputDirectory)
        bed = args.bed
        readLength = args.readLength
        slamSimBed = os.path.join(outputDirectory, replaceExtension(basename(bed), ".bed", "_original"))
        simulator.prepareBED(bed, slamSimBed, readLength)
        
    elif (command == "utrs") :
        outputDirectory = args.outputDir
        createDir(outputDirectory)
        
        bed = args.bed
        trunoverBed = os.path.join(outputDirectory, replaceExtension(basename(bed), ".bed", "_utrs"))
        minHalflife = args.minHalfLife
        maxHalfLife = args.maxHalfLife
        snpRate = args.snpRate
        
        simulator.simulateTurnOver(bed, trunoverBed, minHalflife, maxHalfLife)
        
        referenceFasta = args.referenceFile
        bed12 = os.path.join(outputDirectory, replaceExtension(basename(bed), ".bed12", "_utrs"))
        bed12Fasta = os.path.join(outputDirectory, replaceExtension(basename(bed), ".fa", "_utrs"))
        explv = os.path.join(outputDirectory, replaceExtension(basename(bed), ".eplv", "_utrs"))
        vcfFile = os.path.join(outputDirectory, replaceExtension(basename(bed), ".vcf", "_utrs"))
        
        totalUTRlength = simulator.prepareUTRs(bed, bed12, bed12Fasta, referenceFasta, 50, explv, snpRate, vcfFile)
        
    elif (command == "reads") :
        message("Start simulation")
        
        outputDirectory = args.outputDir
        createDir(outputDirectory)
        
        sampleName = args.sampleName
        readLenght = args.readLength
        
        bed = args.bed
        bed12File = replaceExtension(bed, ".bed12")
        bed12FastaFile = replaceExtension(bed, ".fa")
        explvFile = replaceExtension(bed, ".eplv")
        
        timePoint = args.timePoint
        
        bedReads = os.path.join(outputDirectory, sampleName + "_reads_tmp.bed")
        faReads = os.path.join(outputDirectory, sampleName + "_reads_tmp.fa")
        
        totalUTRlength = simulator.getTotalUtrLength(bed12File)
        readCoverage = args.readCoverage
        readNumber = args.readNumber
        if(readNumber == 0):
            readNumber = (totalUTRlength / readLenght) *  readCoverage
            readNumber = int(readNumber * (random.uniform(-0.2, 0.2) + 1))
        seqError = args.seqError
        
        message("Simulating " + str(readNumber) + " reads with sequencing error of " + str(seqError))
        simulator.simulateReads(bed12File, bed12FastaFile, explvFile, bedReads, faReads, readLenght, readNumber, seqError)
        
        bamReadsWithTC = os.path.join(outputDirectory, sampleName + "_reads.bam")
        utrSummary = os.path.join(outputDirectory, sampleName + "_utrsummary.csv")
        conversionRate = args.conversionRate
        simulator.addTcConversions(bed, faReads, bamReadsWithTC, timePoint, utrSummary, conversionRate, readNumber)
        
        os.unlink(faReads)
        os.unlink(bedReads)    
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
    
if __name__ == '__main__':
    run()

