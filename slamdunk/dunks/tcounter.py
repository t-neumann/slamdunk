#!/usr/bin/env python

# Date located in: -
from __future__ import print_function

import csv
import itertools as IT

from os.path import basename
from utils.misc import getReadCount, getSampleName
from utils.BedReader import BedIterator

from utils import SNPtools
from slamseq.SlamSeqFile import SlamSeqFile, ReadDirection


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
            
    print("gene_name", "non_tc_read_count", "non_tc_norm_read_count", "tc_read_count", "tc_norm_read_count", "tc_read_perc", "fwd_reads", "rev_reads", "snp_In_UTR", sep='\t', file=outCSV)

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

def count(ref, bed, snpsFile, bam, maxReadLength, minQual, outputCSV, log):
    flagstat = getReadCount(bam)
    readNumber = flagstat.MappedReads

    fileCSV = open(outputCSV,'w')
    
    snps = SNPtools.SNPDictionary(snpsFile)

    #Go through one chr after the other
    testFile = SlamSeqFile(bam, ref, snps)
                      
    #chr    start    stop    reads with T->C    read without T->C    Percentage of read with T->C    Number of forward reads mapping to region    Number of reverse reads mapping to region    T->C SNPs found in region
    print("chr", "start", "end", "gene_name", "non_tc_read_count", "non_tc_norm_read_count", "tc_read_count", "tc_norm_read_count", "tc_read_perc", "fwd_reads", "rev_reads", "snp_In_UTR", sep='\t', file=fileCSV)
       
    for utr in BedIterator(bed):
        tcCount = [0] * maxReadLength
                     
        readIterator = testFile.readInRegion(utr.chromosome, utr.start, utr.stop, maxReadLength)
        
        readCount = 0
        countFwd = 0
        countRev = 0
        for read in readIterator:
        
            if(read.direction == ReadDirection.Reverse):
                countRev += 1
            else:
                countFwd += 1
                
            readCount += 1
            tcCount[read.tcCount] += 1
        
        snpInUTR = 0
        if(countRev > countFwd):
            snpInUTR = snps.getAGSNPsInUTR(utr.chromosome, utr.start, utr.stop, 2)
        else:
            snpInUTR = snps.getTCSNPsInUTR(utr.chromosome, utr.start, utr.stop, 1)
        
        percTC = 0
        if(readCount > 0):
            percTC = ((readCount - tcCount[0]) * 1.0 / readCount)
        #chr    start    stop    read without T->C    reads with T->C    Percentage of read with T->C    Number of forward reads mapping to region    Number of reverse reads mapping to region    T->C SNPs found in region
        print(utr.chromosome, utr.start, utr.stop, utr.name, tcCount[0], tcCount[0] * 1000000.0 / readNumber, (readCount - tcCount[0]), (readCount - tcCount[0]) * 1000000.0 / readNumber, percTC, countFwd, countRev, snpInUTR, sep='\t', file=fileCSV)
    
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

    
