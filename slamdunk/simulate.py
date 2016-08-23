#!/usr/bin/env python

#########################################################################
# Main routine for the SLAMdunk simulation
#########################################################################
# Imports
#########################################################################

from __future__ import print_function
import sys, os, random

from argparse import ArgumentParser, RawDescriptionHelpFormatter
    
from os.path import basename

from dunks import simulator
from joblib import Parallel, delayed
from utils.misc import replaceExtension, readSampleNames, checkStep

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
    if not os.path.exists(directory):
        message("Creating output directory: " + directory)
        os.makedirs(directory)



def run():
    ########################################################################
    # Argument parsing
    ########################################################################
    
    # Info
    usage = "SLAMdunk software for simulating SLAM-seq data"
    version = "1.0"
    
    # Main Parsers
    parser = ArgumentParser(description=usage, formatter_class=RawDescriptionHelpFormatter, version=version)
    
    # Initialize Subparsers
    subparsers = parser.add_subparsers(help="", dest="command")
    
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
    simulateparse.add_argument("-tc", "--tc-rate", type=float, required=False, dest="conversionRate", default=0.2, help="T->C conversion rate")
        
    args = parser.parse_args()
    
    ########################################################################
    # Routine selection
    ########################################################################
    
    command = args.command
    
    if (command == "utrs") :
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
        
        
    
if __name__ == '__main__':
    run()
