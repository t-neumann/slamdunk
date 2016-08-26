#!/usr/bin/env python

# Date located in: -
from __future__ import print_function
#import pysam, random, os, sys
import random
import pysam
import math
import os
import glob
import sys

from utils import SNPtools
from utils.BedReader import BedIterator
from utils.misc import shell, run
from slamseq.SlamSeqFile import SlamSeqBamFile
from Bio import SeqIO

projectPath = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
rNASeqReadSimulatorPath = os.path.join(projectPath, "bin", "RNASeqReadSimulator-master/")
pathEvalHalfLifes = os.path.join(projectPath, "plot", "eval_halflife_per_gene_plots.R")
pathEvalConversionrates = os.path.join(projectPath, "plot", "eval_conversion_rate_plots.R")


def getRndBaseWithoutDup(base):
    rndBase = getRndBase()
    while(rndBase == base):
        rndBase = getRndBase()
    return rndBase
        
def getRndBase():
    return ['A', 'T', 'G', 'C'][random.randrange(0, 3, 1)]

def getCmpBase(base):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'} 
    return complement[base]

def printUTR(utr, outBed, minLength):
    if utr.getLength() > minLength:
        print(utr.chromosome, utr.start, utr.stop, utr.name, utr.score, utr.strand, sep="\t", file=outBed)

def prepareBED(bed, slamSimBed, minLength):
    utrs = []
    for utr in BedIterator(bed):
        utrs.append(utr)
        
    utrs.sort(key=lambda x: (x.name, -x.getLength()))
    
    outBed = open(slamSimBed, "w")
    
    partList = []
    lastUtr = None
    for utr in utrs:
        if utr.hasStrand() and utr.hasNonEmptyName():
            currentUtr = utr.name
            if currentUtr == lastUtr:
                partList.append(utr)
            else:
                if(not lastUtr is None):
                    printUTR(partList[0], outBed, minLength)
                partList = [utr]
            lastUtr = currentUtr
        else:
            print("Warning: Invalid BED entry found: " + str(utr))
        
    if(not lastUtr is None):
        printUTR(partList[0], outBed, minLength)
    
    outBed.close()

def simulateUTR(utrSeq, utr, polyALenght, snpRate, vcfFile):
    utrSeq = list(utrSeq)
    snpCount = 0
    for i in xrange(0, len(utrSeq)):
        if(random.uniform(0, 1) < snpRate):
            rndBase = getRndBaseWithoutDup(utrSeq[i])
            print("Introducing SNP " + utrSeq[i] + " -> " + rndBase + " in UTR " + utr.name + " at position " + utr.chromosome + ":" + str(utr.start + i))
            snpPosition = 0
            snpCount += 1
            if(utr.strand == "+"):
                snpPosition = utr.start + i + 1
                snpRef = utrSeq[i]
                snpAlt = rndBase
            else:
                snpPosition = utr.stop - i
                snpRef = getCmpBase(utrSeq[i])
                snpAlt = getCmpBase(rndBase)
                
            print(utr.chromosome, snpPosition, utr.name + "_" + str(snpCount), snpRef, snpAlt, ".", "PASS", ".", sep="\t", file=vcfFile)
            utrSeq[i] = rndBase
    return "".join(utrSeq) + (polyALenght * 'A')

def prepareUTRs(bed, bed12, bed12Fasta, referenceFasta, readLength, explv, snpRate, vcfFile):
    
    # Read utrs from BED file
    utrs = parseUtrBedFile(bed)
    
    vcf = open(vcfFile, "w")
    print("##fileformat=VCFv4.1", file=vcf)
    print("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO", file=vcf)
    
    bed12FastaTmp = bed12Fasta + "_tmp.fa"
    utrFasta = shell("bedtools getfasta -name -s -fi " + referenceFasta + " -bed " + bed + " -fo " + bed12FastaTmp)
    print(utrFasta)
    
    bed12FastaTmpFile = open(bed12FastaTmp, "r")
    utrFasta = bed12FastaTmpFile.read()
    bed12FastaTmpFile.close()
    os.unlink(bed12FastaTmp)
    
    bed12FastaFile = open(bed12Fasta, "w")
    utrName = None
    for line in utrFasta.splitlines():
        if(line[0] == ">"):
            print(line, file=bed12FastaFile)
            utrName = line[1:] 
        else:
            print(simulateUTR(line, utrs[utrName], readLength, snpRate, vcf), file=bed12FastaFile)
    bed12FastaFile.close()
    vcf.close()
    
    bed12File = open(bed12, "w")
    
    totalLength = 0
    
    minFragmentLength = 150
    maxFragmentLength = 250
    for utr in BedIterator(bed):
        
        fragmentLength = random.randrange(minFragmentLength, maxFragmentLength, 1)
        
        start = max(0, utr.getLength() - fragmentLength)
        end = utr.getLength() + readLength

        totalLength += (end - start)
        print(utr.name, start, end, utr.name, utr.score, "+", start, end, "255,0,0", "1", min(utr.getLength() + readLength / 4, fragmentLength + readLength / 4), 0, sep="\t", file=bed12File)
        
    bed12File.close()    
    
    output = shell(rNASeqReadSimulatorPath + "src/genexplvprofile.py --geometric 0.8 " + bed12 + " > " + explv)
    print(output)
        
    return totalLength
    
def simulateReads(bed12, bed12Fasta, explv, bedReads, faReads, readLength, readCount, seqError):
    
    output = shell(rNASeqReadSimulatorPath + "src/gensimreads.py -l " + str(readLength) + " -e " + explv + " -n " + str(readCount) + " -b " + rNASeqReadSimulatorPath + "demo/input/sampleposbias.txt --stranded " + bed12 + " > " + bedReads)
    print(output)
    output = shell(rNASeqReadSimulatorPath + "src/getseqfrombed.py -r " + str(seqError) + " -l " + str(readLength) + " " + bedReads + " " + bed12Fasta + " > " + faReads)
    print(output)
    
def getRndHalfLife(minHalfLife, maxHalfLife):
    return random.randrange(minHalfLife, maxHalfLife, 1)

def simulateTurnOver(bed, turnoverBed, minHalfLife, maxHalfLife):
    turnoverFile = open(turnoverBed, "w")
    for utr in BedIterator(bed):
        print(utr.chromosome, utr.start, utr.stop, utr.name, getRndHalfLife(minHalfLife, maxHalfLife), utr.strand, sep='\t', file=turnoverFile)
    turnoverFile.close()

def printFastaEntry(sequence, name, index, conversions, readOutSAM):
    #a = pysam.AlignedSegment()
    print(name + "_" + str(index) + "_" + str(conversions),
          "4",
          "*",
          "0",
          "0",
          "*",
          "*",
          "0",
          "0",
          sequence,
          "<" * len(sequence),
          "TC:i:" + str(conversions),
          "ID:i:" + str(index),
           file=readOutSAM, sep="\t")
    
def convertRead(read, name, index, conversionRate, readOutSAM):
    
    tCount = 0
    TcCount = 0
    # TODO: Uncomment to go back to pysam
    #seq = list(read.sequence)
    seq = list(read.seq.tostring())
    for i in xrange(0, len(seq)):
        if seq[i] == 'T':
            tCount += 1
            if random.uniform(0, 1) <= conversionRate:
                seq[i] = 'C'
                TcCount += 1
    
    printFastaEntry("".join(seq), name, index, TcCount, readOutSAM)
    
    return tCount, TcCount

def getLambdaFromHalfLife(halfLife):
    return math.log(2) / float(halfLife)

def addTcConversionsToReads(utr, reads, timePoint, readOutSAM, conversionRate):    
    print(utr.name + " reads found: " + str(len(reads)))
    
    varLambda = getLambdaFromHalfLife(utr.score)
    readsToConvert = int(len(reads) * (1 - math.exp(-varLambda * timePoint)))
    print("Converting " + str(readsToConvert) + " reads (lambda = " + str(varLambda) + ")")
    
    totalTCount = 0
    totalTcCount = 0
    readSample = sorted(random.sample(xrange(0, len(reads)), readsToConvert))
    #print(readSample)
    #print("Reads selected: " + str(len(readSample)))
    sampledReads = 0
    notSampledReads = 0
    for i in xrange(0, len(reads)):
        read = reads[i]
        if(i in readSample):
            tCount, TcCount = convertRead(read, utr.name, i, conversionRate, readOutSAM)
            sampledReads += 1
        else:
            tCount, TcCount = convertRead(read, utr.name, i, 0, readOutSAM)
            notSampledReads += 1
        totalTcCount += TcCount
        totalTCount += tCount
    
    #print(sampledReads, notSampledReads, len(reads), sampledReads * 1.0 / len(reads),  totalTcCount, totalTCount, totalTcCount * 1.0 / totalTCount )
    return readsToConvert, totalTCount, totalTcCount

def getUtrName(readName):
    return readName.split("_")[0]

def printUtrSummary(utr, totalReadCount, readsToConvert, totalTCount, totalTcCount, utrSummary, readsCPM):
    conversionRate = 0
    if totalTCount > 0:
        conversionRate = totalTcCount * 1.0 / totalTCount
    #print(utr.chromosome, utr.start, utr.stop, utr.name, utr.score, utr.strand, totalReadCount, readsToConvert, totalTCount, totalTcCount, conversionRate, sep="\t", file=utrSummary)
    print(utr.chromosome, 
          utr.start, 
          utr.stop, 
          utr.name, #utr.score, 
          utr.strand,
          conversionRate,
          readsCPM,
          totalTCount,
          totalTcCount,
          totalReadCount, 
          readsToConvert,
          "-1" , sep="\t", file=utrSummary)

def parseUtrBedFile(bed):
    utrs = {}
    for utr in BedIterator(bed):
        utrs[utr.name] = utr
    return utrs

def addTcConversions(bed, readInFile, readOutFile, timePoint, utrSummaryFile, conversionRate, librarySize):
    
    # Read utrs from BED file
    utrs = parseUtrBedFile(bed)
    
    readOutTemp = readOutFile + "_tmp.sam"
    #bamheader = { 'HD': {'VN': '1.0'} }
    #readOutBAM = pysam.AlignmentFile(readOutTemp, "wb", header=bamheader, add_sq_text=False)
    readOutSAM = open(readOutTemp, "w")
    utrSummary = open(utrSummaryFile, "w")
    
    reads = []
    lastUtrName = None
    utrName = None
    
    fasta_sequences = SeqIO.parse(open(readInFile),'fasta')

    for entry in fasta_sequences:
        
    # TODO: Uncomment to go back to pysam
    #with pysam.FastxFile(readInFile) as fh:
        #for entry in fh:
            #utrName = getUtrName(entry.name)
        utrName = getUtrName(entry.id)
        if(utrName == lastUtrName):
            reads.append(entry)
        elif(lastUtrName == None):
            reads.append(entry)
        else:
            readsCPM = len(reads)  * 1000000.0 / librarySize;
            readsToConvert, totalTCount, totalTcCount = addTcConversionsToReads(utrs[lastUtrName], reads, timePoint, readOutSAM, conversionRate)
            printUtrSummary(utrs[lastUtrName], len(reads), readsToConvert, totalTCount, totalTcCount, utrSummary, readsCPM)
            reads = []
        lastUtrName = utrName
    readsCPM = len(reads) * 1000000.0 / librarySize;
    readsToConvert, totalTCount, totalTcCount = addTcConversionsToReads(utrs[lastUtrName], reads, timePoint, readOutSAM, conversionRate)
    printUtrSummary(utrs[lastUtrName], len(reads), readsToConvert, totalTCount, totalTcCount, utrSummary, readsCPM)
        
            
    readOutSAM.close()       
    utrSummary.close()  
    
    
    readOutTempBAM = readOutFile + "_tmp.bam"
    run("samtools view -Sb " + readOutTemp + " > " + readOutTempBAM)
    # Sort reads by query name (doesn't matter for mapping, but makes evaluation easier
    run("samtools sort -o " + readOutFile + " " + readOutTempBAM)
    os.unlink(readOutTemp)
    os.unlink(readOutTempBAM)
            
def getTotalUtrLength(bed12File):
    totalUtrLength = 0
    for utr in BedIterator(bed12File):
        totalUtrLength += utr.getLength()
    return totalUtrLength
        
        
def plotconversiondifferences(simDir, slamDir, conversionRate, outputPDF):
    
    simFiles = sorted(glob.glob(simDir + "*_utrsummary.csv"))
    slamdunkFiles = sorted(glob.glob(slamDir + "*_reads_slamdunk_mapped_filtered_tcount.csv"))
    
    if(len(simFiles) == len(slamdunkFiles)):
        run("Rscript " + pathEvalConversionrates + " -c " + str(conversionRate) + " -s " + ",".join(simFiles) + " -f " + ",".join(slamdunkFiles) + " -o " + outputPDF, sys.stderr, dry=False, verbose=False)
    else:
        raise RuntimeError("Couldn't match files with timepoints")
    

def plotHalfLifes(bed, simDir, slamDir, timePointsStr, conversionRate, outputPDF):
    
    simFiles = sorted(glob.glob(simDir + "*_utrsummary.csv"))
    slamdunkFiles = sorted(glob.glob(slamDir + "*_reads_slamdunk_mapped_filtered_tcount.csv"))
    timePoints = timePointsStr.split(",")
    
    if(len(simFiles) == len(slamdunkFiles) and len(slamdunkFiles) == len(timePoints)):
        run("Rscript " + pathEvalHalfLifes + " -b " + bed + " -c " + str(conversionRate) + " -s " + ",".join(simFiles) + " -f " + ",".join(slamdunkFiles) + " -t " + ",".join(timePoints) + " -o " + outputPDF, sys.stderr, dry=False, verbose=False)
    else:
        raise RuntimeError("Couldn't match files with timepoints")
    
def getConversionRateFromBam(bam, ref, chromosome, start, end, strand):
    
    testFile = SlamSeqBamFile(bam, ref, SNPtools.SNPDictionary(None))
    
    sumConversionRate = 0
    readCount = 0
    #for chromosome in testFile.getChromosomes():
    #    readIterator = testFile.readsInChromosome(chromosome)
    #readIterator = testFile.readInRegion("chr7", 3217778, 3221036, "+", 55)
    readIterator = testFile.readInRegion(chromosome, start, end, strand, 100)
        
    for read in readIterator:
        conversionRate = 0
        if(read.tCount > 0):
            conversionRate = read.tcCount * 1.0 / read.tCount
        
        #if(read.tcCount > 0):
        sumConversionRate += conversionRate
        readCount += 1
        
        #if(readCount % 1000 == 0 and readCount > 0):
        #    print(str(readCount) + ": " + str(sumConversionRate) + " / " + str(readCount) + " = " + str(sumConversionRate / readCount))
        #if(readCount >= 10000):
        #    break
    
    print("Read count: " + str(readCount))
    print("Avg. conversion rate: " + str(sumConversionRate / readCount))