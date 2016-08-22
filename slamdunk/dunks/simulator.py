#!/usr/bin/env python

# Date located in: -
from __future__ import print_function
#import pysam, random, os, sys
import random
import pysam
import math

from utils.BedReader import BedIterator
from utils.misc import shell

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
    
    utrFasta = shell("bedtools getfasta -name -s -fi " + referenceFasta + " -bed " + bed + " -fo - ")
    print(utrFasta)
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
    
    output = shell("genexplvprofile.py --geometric 0.8 " + bed12 + " > " + explv)
    print(output)
        
    return totalLength
    
def simulateReads(bed12, bed12Fasta, explv, bedReads, faReads, readLength, readCount, seqError):
    
    output = shell("gensimreads.py -l " + str(readLength) + " -e " + explv + " -n " + str(readCount) + " -b /project/ngs/philipp/slamdunk-analysis/simulation/tools/RNASeqReadSimulator/demo/input/sampleposbias.txt --stranded " + bed12 + " > " + bedReads)
    print(output)
    output = shell("getseqfrombed.py -r " + str(seqError) + " -l " + str(readLength) + " " + bedReads + " " + bed12Fasta + " > " + faReads)
    print(output)
    
def simulateTurnOver(bed, turnoverBed, minTurnover, maxTurnover):
    turnoverFile = open(turnoverBed, "w")
    for utr in BedIterator(bed):
        print(utr.chromosome, utr.start, utr.stop, utr.name, random.uniform(minTurnover, maxTurnover), utr.strand, sep='\t', file=turnoverFile)
    turnoverFile.close()

def printFastaEntry(sequence, name, index, conversions, outputBAM):
    a = pysam.AlignedSegment()
    a.query_name = name + "_" + str(index) + "_" + str(conversions)
    a.flag = 4
    a.query_sequence = sequence
    a.query_qualities = pysam.qualitystring_to_array("<" * len(sequence))
    a.tags = (("TC", conversions),
              ("ID", index))
    outputBAM.write(a)
    
def convertRead(read, name, index, conversionRate, outputBAM):
    
    tCount = 0
    TcCount = 0
    seq = list(read.sequence)
    for i in xrange(0, len(seq)):
        if seq[i] == 'T':
            tCount += 1
            if random.random() <= conversionRate:
                seq[i] = 'C'
                TcCount += 1
    
    printFastaEntry("".join(seq), name, index, TcCount, outputBAM)
    
    return tCount, TcCount
    

def addTcConversionsToReads(utr, reads, timePoint, outputBAM):    
    conversionRate = 0.01
    print(utr.name + " reads found: " + str(len(reads)))
    
    varLambda = float(utr.score)
    readsToConvert = int(len(reads) * (1 - math.exp(-varLambda * timePoint)))
    print("Converting " + str(readsToConvert) + " reads (lambda = " + str(varLambda) + ")")
    
    totalTCount = 0
    totalTcCount = 0
    readSample = random.sample(range(0, len(reads)), readsToConvert)
    for i in xrange(0, len(reads)):
        read = reads[i]
        if(i in readSample):
            tCount, TcCount = convertRead(read, utr.name, i, conversionRate, outputBAM)
        else:
            tCount, TcCount = convertRead(read, utr.name, i, 0, outputBAM)
        totalTcCount += TcCount
        totalTCount += tCount
    
    return readsToConvert, totalTCount, totalTcCount

def getUtrName(readName):
    return readName.split("_")[0]

def printUtrSummary(utr, totalReadCount, readsToConvert, totalTCount, totalTcCount, utrSummary):
    conversionRate = 0
    if totalTCount > 0:
        conversionRate = totalTcCount * 1.0 / totalTCount
    print(utr.chromosome, utr.start, utr.stop, utr.name, utr.score, utr.strand, totalReadCount, readsToConvert, totalTCount, totalTcCount, conversionRate, sep="\t", file=utrSummary)

def parseUtrBedFile(bed):
    utrs = {}
    for utr in BedIterator(bed):
        utrs[utr.name] = utr
    return utrs

def addTcConversions(bed, readInFile, readOutFile, timePoint, utrSummaryFile):
    
    # Read utrs from BED file
    utrs = parseUtrBedFile(bed)
    
    bamheader = { 'HD': {'VN': '1.0'} }
    readOutBAM = pysam.AlignmentFile(readOutFile, "wb", header=bamheader)
    utrSummary = open(utrSummaryFile, "w")
    
    reads = []
    lastUtrName = None
    utrName = None
    with pysam.FastxFile(readInFile) as fh:
        for entry in fh:
            utrName = getUtrName(entry.name)
            if(utrName == lastUtrName):
                reads.append(entry)
            elif(lastUtrName == None):
                reads.append(entry)
            else:
                readsToConvert, totalTCount, totalTcCount = addTcConversionsToReads(utrs[lastUtrName], reads, timePoint, readOutBAM)
                printUtrSummary(utrs[lastUtrName], len(reads), readsToConvert, totalTCount, totalTcCount, utrSummary)
                reads = []
            lastUtrName = utrName
        readsToConvert, totalTCount, totalTcCount = addTcConversionsToReads(utrs[lastUtrName], reads, timePoint, readOutBAM)
        printUtrSummary(utrs[lastUtrName], len(reads), readsToConvert, totalTCount, totalTcCount, utrSummary)
        
            
    readOutBAM.close()       
    utrSummary.close()    
            
def getTotalUtrLength(bed12File):
    totalUtrLength = 0
    for utr in BedIterator(bed12File):
        totalUtrLength += utr.getLength()
    return totalUtrLength
        