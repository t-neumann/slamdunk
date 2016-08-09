#!/usr/bin/env python

# Date located in: -
from __future__ import print_function

import csv
import sys
import itertools as IT

from os.path import basename
from utils.misc import getReadCount, getSampleName
from utils.BedReader import BedIterator

from utils import SNPtools
from slamseq.SlamSeqFile import SlamSeqBamFile, ReadDirection, SlamSeqInterval


def collapse(expandedCSV, collapsedCSV, readNumber, log):
    
    tcDict = {}
    
    outCSV = open(collapsedCSV, 'w')
    
    with open(expandedCSV, 'r') as f:
        
        # Skip header
        next(f)

        for line in f:
            fields = line.rstrip().split('\t')
            if (len(fields) == 12) :
            
                gene = fields[3]
                nontc = fields[4]
                tc = fields[6]
                fwdReads = fields[9]
                revReads = fields[10]
                snps = fields[11]
                
                if (gene in tcDict.keys()) :
                    tcDict[gene]['nontc'] += int(nontc)
                    tcDict[gene]['tc'] += int(tc)
                    tcDict[gene]['fwdReads'] += int(fwdReads)
                    tcDict[gene]['revReads'] += int(revReads)
                    tcDict[gene]['snps'] += int(snps)
                else :
                    tcDict[gene] = {}
                    tcDict[gene]['nontc'] = int(nontc)
                    tcDict[gene]['tc'] = int(tc)
                    tcDict[gene]['fwdReads'] = int(fwdReads)
                    tcDict[gene]['revReads'] = int(revReads)
                    tcDict[gene]['snps'] = int(snps)
            
        else :
            print("Error in TC file format - unexpected number of fields (" + str(len(fields)) + ") in the following line:\n" + line, file=log)
            
    print("gene_name", "non_tc_read_count", "non_tc_norm_read_count", "tc_read_count", "tc_norm_read_count", "tc_read_perc", "avg_conversion_rate", "fwd_reads", "rev_reads", "snp_In_UTR", sep='\t', file=outCSV)

    for gene in sorted(tcDict.keys()) :
        
        print(gene,end="\t",file=outCSV)
        print(tcDict[gene]['nontc'],end="\t",file=outCSV)
        print(float(tcDict[gene]['nontc']) / float(readNumber) * float(1000000),end="\t",file=outCSV)
        print(tcDict[gene]['tc'],end="\t",file=outCSV)
        print(float(tcDict[gene]['tc']) / float(readNumber) * float(1000000),end="\t",file=outCSV)
        
        percent = 0
        if (tcDict[gene]['nontc'] + tcDict[gene]['tc'] > 0) :
            percent = float(tcDict[gene]['tc']) / float((tcDict[gene]['nontc'] + tcDict[gene]['tc']))
        elif (tcDict[gene]['tc'] > 0) :
            percent = 1
        print(percent,end="\t",file=outCSV)
        
        print(tcDict[gene]['fwdReads'],end="\t",file=outCSV)
        print(tcDict[gene]['revReads'],end="\t",file=outCSV)
        print(tcDict[gene]['snps'],file=outCSV)
        
        #sys.stdin.readline()
        
    outCSV.close()    
    

def getMean(values, skipZeros=True):
    count = 0.0
    totalSum = 0.0
    for value in values:
        if not skipZeros or value > 0:
            totalSum = totalSum + value
            count += 1
    if(count > 0):
        return totalSum / count
    else:
        return 0.0

def computeTconversions(ref, bed, snpsFile, bam, maxReadLength, minQual, outputCSV, outputBedgraphPlus, outputBedgraphMinus, strictTCs, log):
    
    flagstat = getReadCount(bam)
    readNumber = flagstat.MappedReads

    fileCSV = open(outputCSV,'w')
    
    snps = SNPtools.SNPDictionary(snpsFile)
    snps.read()
    
    #Go through one chr after the other
    testFile = SlamSeqBamFile(bam, ref, snps)
    
    conversionBedGraph = {}
                         
    progress = 0
    for utr in BedIterator(bed):

        if(not utr.hasStrand()):
            raise RuntimeError("Input BED file does not contain stranded intervals.")

        readIterator = testFile.readInRegion(utr.chromosome, utr.start, utr.stop, utr.strand, maxReadLength)
      
        tcCountUtr = [0] * utr.getLength()
        coverageUtr = [0] * utr.getLength()

        countFwd = 0
        tcCountFwd = 0
        countRev = 0
        tCountRev = 0
        
        multiMapFwd = 0
        multiMapRev = 0
        
        for read in readIterator:
            
            # Overwrite any conversions for non-TC reads (reads with < 2 TC conversions)
            if (not read.isTcRead and strictTCs) :
                read.tcCount = 0
                read.mismatches = []
                read.conversionRates = 0.0
                read.tcRate = 0.0
            
            if(read.direction == ReadDirection.Reverse):
                countRev += 1
                if read.tcCount > 0:
                    tCountRev += 1
                if read.isMultimapper:
                    multiMapRev += 1
            else:
                countFwd += 1
                if read.tcCount > 0:
                    tcCountFwd += 1
                if read.isMultimapper:
                    multiMapFwd += 1
            
            for mismatch in read.mismatches:
                if(mismatch.isTCMismatch(read.direction == ReadDirection.Reverse) and mismatch.referencePosition >= 0 and mismatch.referencePosition < utr.getLength()):
                    tcCountUtr[mismatch.referencePosition] += 1
            
            for i in xrange(read.startRefPos, read.endRefPos):
                if(i >= 0 and i < utr.getLength()):
                    coverageUtr[i] += 1
            

        if((utr.strand == "+" and countFwd > 0) or (utr.strand == "-" and countRev > 0)):        
            tcRateUtr = [ x * 100.0 / y if y > 0 else 0 for x, y in zip(tcCountUtr, coverageUtr)]
            
            readCount = countFwd
            tcReadCount = tcCountFwd
            multiMapCount = multiMapFwd
            
            if(utr.strand == "-"):
                readCount = countRev
                tcReadCount = tCountRev
                multiMapCount = multiMapRev
            
            if((utr.strand == "-" and countFwd > countRev) or (utr.strand == "+" and countRev > countFwd)):
                print("Warning: " + utr.name + " is located on the " + utr.strand + " strand but read counts are higher for the opposite strand (fwd: " + countFwd + ", rev: " + countRev + ")", file=sys.stderr)
                
            
            refSeq = readIterator.getRefSeq()
    
            # Get number of covered Ts/As in the UTR and compute average conversion rate for all covered Ts/As
            coveredTcount = 0
            avgConversationRate = 0
            coveredPositions = 0
            # Get number of reads on T positions and number of reads with T->C conversions on T positions
            coverageOnTs = 0
            conversionsOnTs = 0
            
            for position in xrange(0, len(coverageUtr)):
#                 if position >= len(refSeq) :
#                     print(utr)
#                     print(position,file=sys.stderr)
#                     print(len(refSeq),file=sys.stderr)
#                     print(refSeq,file=sys.stderr)
#                     sys.stdin.readline()
                if(coverageUtr[position] > 0 and ((utr.strand == "+" and refSeq[position] == "T") or (utr.strand == "-" and refSeq[position] == "A"))):
                    coveredTcount += 1
                    avgConversationRate += tcRateUtr[position]
                    
                    coverageOnTs += coverageUtr[position]
                    conversionsOnTs += tcCountUtr[position]
                    #print(refSeq[position] + "\t", utr.start + position - 1, refSeq[position], coverageUtr[position], tcRateUtr[position], file=sys.stderr)
                    # print conversion rates for all covered Ts/As
                    #print(utr.chromosome, utr.start + position - 1, utr.start + position, tcRateUtr[position], sep="\t")
                    conversionBedGraph[utr.chromosome + ":" + str(utr.start + position - 1) + ":" + str(utr.strand)] = tcRateUtr[position] 
                if(coverageUtr[position] > 0):
                    coveredPositions += 1
            
            if(coveredTcount > 0):
                avgConversationRate = avgConversationRate / coveredTcount
            else:
                avgConversationRate = 0
                
            # reads per million mapped to the UTR
            readsCPM = 0
            if(readNumber > 0):
                readsCPM = readCount  * 1000000.0 / readNumber
            
         
            # Convert to SlamSeqInterval and print
            #slamSeqUtr = SlamSeqInterval(utr.chromosome, utr.start, utr.stop, utr.strand, utr.name, readsCPM, avgConversationRate, coveredTcount, coveredPositions, readCount, tcReadCount)
            conversionRate = 0
            if (coverageOnTs > 0) :
                conversionRate = float(conversionsOnTs) / float(coverageOnTs)
            slamSeqUtr = SlamSeqInterval(utr.chromosome, utr.start, utr.stop, utr.strand, utr.name, readsCPM, coverageOnTs, conversionsOnTs, conversionRate, readCount, tcReadCount, multiMapCount)
        else:
            slamSeqUtr = SlamSeqInterval(utr.chromosome, utr.start, utr.stop, utr.strand, utr.name, 0, -1, 0, 0, 0, 0, 0)
            #slamSeqUtr = SlamSeqInterval(utr.chromosome, utr.start, utr.stop, utr.strand, utr.name, 0, 0, 0, 0, 0, 0)
        print(slamSeqUtr, file=fileCSV)
        #print(slamSeqUtr)
        
        #if progress % 10 == 0:
            #print("Progress: " + str(progress) + "\r", file=sys.stderr, end="")
        progress += 1
        
    fileCSV.close()
    
    fileBedgraphPlus = open(outputBedgraphPlus,'w')
    fileBedgraphMinus = open(outputBedgraphMinus,'w')
    
    for position in conversionBedGraph:
        positionData = position.split(":")
        if(positionData[2] == "+"):
            print(positionData[0], positionData[1], int(positionData[1]) + 1, conversionBedGraph[position], file=fileBedgraphPlus)
        else:
            print(positionData[0], positionData[1], int(positionData[1]) + 1, conversionBedGraph[position], file=fileBedgraphMinus)
            
    fileBedgraphPlus.close()
    fileBedgraphMinus.close()

# @deprecated 
def count(ref, bed, snpsFile, bam, maxReadLength, minQual, outputCSV, log):
    flagstat = getReadCount(bam)
    readNumber = flagstat.MappedReads

    fileCSV = open(outputCSV,'w')
    
    snps = SNPtools.SNPDictionary(snpsFile)
    snps.read()

    #Go through one chr after the other
    testFile = SlamSeqBamFile(bam, ref, snps)
                      
    #chr    start    stop    reads with T->C    read without T->C    Percentage of read with T->C    Number of forward reads mapping to region    Average conversion rate for all reads mapping to the utr    Number of reverse reads mapping to region    T->C SNPs found in region
    print("chr", "start", "end", "gene_name", "non_tc_read_count", "non_tc_norm_read_count", "tc_read_count", "tc_norm_read_count", "tc_read_perc", "avg_conversion_rate", "fwd_reads", "rev_reads", "snp_In_UTR", sep='\t', file=fileCSV)
       
    for utr in BedIterator(bed):
        tcCount = [0] * maxReadLength
    
        readIterator = testFile.readInRegion(utr.chromosome, utr.start, utr.stop, utr.strand, maxReadLength)
        
        readCount = 0
        countFwd = 0
        countRev = 0
        avgConversionRate = 0.0
        for read in readIterator:
        
            if(read.direction == ReadDirection.Reverse):
                countRev += 1
            else:
                countFwd += 1
                
            readCount += 1
            tcCount[read.tcCount] += 1
            avgConversionRate += read.tcRate
        
        snpInUTR = 0
        if(countRev > countFwd):
            snpInUTR = snps.getAGSNPsInUTR(utr.chromosome, utr.start, utr.stop, 2)
        else:
            snpInUTR = snps.getTCSNPsInUTR(utr.chromosome, utr.start, utr.stop, 1)
        
        percTC = 0
        if(readCount > 0):
            percTC = ((readCount - tcCount[0]) * 1.0 / readCount)
        #chr    start    stop    read without T->C    reads with T->C    Percentage of read with T->C    Number of forward reads mapping to region    Number of reverse reads mapping to region    T->C SNPs found in region
        if(readCount > 0):
            avgConversionRate = avgConversionRate * 1.0 / readCount
        else:
            avgConversionRate = 0.0
        print(utr.chromosome, utr.start, utr.stop, utr.name, tcCount[0], tcCount[0] * 1000000.0 / readNumber, (readCount - tcCount[0]), (readCount - tcCount[0]) * 1000000.0 / readNumber, percTC, avgConversionRate, countFwd, countRev, snpInUTR, sep='\t', file=fileCSV)
    
    fileCSV.close()


def summary(bams, samples, outputFile, colNumber):
    filenames = bams
    handles = [open(filename, 'rb') for filename in filenames]    
    readers = [csv.reader(f, delimiter='\t') for f in handles]
    names = ["chr", "start", "end", "gene"] + [getSampleName(basename(f), samples) for f in filenames]

    with  open(outputFile, 'wb') as h:
        writer = csv.writer(h, delimiter=';', lineterminator='\n', )
        writer.writerow(names)
        header = True
        for rows in IT.izip_longest(*readers, fillvalue=['']*2):
            if(not header):
                combined_row = []
                firstRow = True
                for row in rows:
                    if(firstRow):
                        combined_row.extend(row[0:4])
                        combined_row.append(row[colNumber])
                        
                        firstRow = False
                    else:
                        combined_row.append(row[colNumber])
            
                writer.writerow(combined_row)
            else:
                header = False
    
    for f in handles:
        f.close()

    
