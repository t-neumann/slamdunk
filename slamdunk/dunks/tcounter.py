#!/usr/bin/env python

# Date located in: -
from __future__ import print_function
import sys

import pysam
import csv
import itertools as IT

from utils import removeExtension
from os.path import basename
from dunks.utils import files_exist

snpDict = {}

def getSampleName(fileName, samples):
    if samples == None:
        return removeExtension(fileName)
    else:
        for key in samples:
            if(key in fileName):
                return samples[key]
        
    return

def isAGSnp(chromosome, position):
    key = chromosome + str(int(position) + 1)
    return key in snpDict and snpDict[key] == 2

def getSNPsInUTR(chromosome, start, stop, snpType):
    count = 0
    for i in range(start, stop):
        if(snpType == 1):
            if(isTCSnp(chromosome, i)):
                count += 1
        else:
            if(isAGSnp(chromosome, i)):
                count += 1
    return count

def isTCSnp(chromosome, position):
    key = chromosome + str(int(position) + 1)
    return key in snpDict and snpDict[key] == 1

def getTCNgm(read):
    return int(read.get_tag("TC"))

def getTC(read, refSeq, chromosome, refRegion, regionStart, maxReadLength, minQual):
#     print(read, file=sys.stderr)
    tcCount = 0;
    agCount = 0;
    for pair in read.get_aligned_pairs(matches_only=True):
        readPos = pair[0]
        refPos = pair[1] - int(regionStart) + maxReadLength + 1
        refBase = refSeq[refPos]
        readBase = read.query_sequence[readPos]
        readQlty = read.query_qualities[readPos]
        if readBase != refBase:
#             print(pair, readPos, refPos, readBase, refBase, readQlty, file=sys.stderr)
            if(readQlty >= minQual):
                if(read.is_reverse):
                    if(refBase == 'A' and readBase == 'G'):
                        if(not isAGSnp(chromosome, int(pair[1]))):
                            agCount += 1
#                         print(pair, readPos, refPos, readBase, refBase, readQlty, file=sys.stderr)
                else:
                    if(refBase == 'T' and readBase == 'C'):
                        if(not isTCSnp(chromosome, int(pair[1]))):
                            tcCount += 1
#                         print(pair, readPos, refPos, readBase, refBase, readQlty, file=sys.stderr)
                    
    return [tcCount, agCount]

def countReads(bam):
    # TODO
    return 0

def getMappedReadCount(bam):
    flagStat = bam + ".flagstat"
    readNumber = 0
    if(files_exist(flagStat)):
        with open(flagStat, 'r') as content_file:
            content = content_file.read()
            if content.count("\n") > 10:
                readNumber = int(content.split("\n")[4].split(" ")[0])
    else:
        readNumber = countReads(bam)
    return readNumber

def count(ref, bed, snps, bam, maxReadLength, minQual, outputCSV, log):
    
    readNumber = getMappedReadCount(bam)
#     print(readNumber)
#     return
    fileCSV = open(outputCSV,'w')
    
    reffile = pysam.FastaFile(ref)
    samfile = pysam.AlignmentFile(bam, "rb")
    
    if not snps is None:
        with open(snps, "r") as f:
            for line in f:
                #print(line)
                cols = line.rstrip().split("\t")
                #print(cols[2], cols[3])
                if(cols[2] == "A" and cols[3] == "G"):
                    key = cols[0] + cols[1]
                    snpDict[key] = 2
                if(cols[2] == 'T' and cols[3] == 'C'):
                    key = cols[0] + cols[1]
                    snpDict[key] = 1
    
    
    lineCount = 0
    
    #chr    start    stop    reads with T->C    read without T->C    Percentage of read with T->C    Number of forward reads mapping to region    Number of reverse reads mapping to region    T->C SNPs found in region
    print("chr", "start", "end", "gene_name", "non_tc_read_count", "non_tc_norm_read_count", "tc_read_count", "tc_norm_read_count", "tc_read_perc", "fwd_reads", "rev_reads", "snp_In_UTR", sep='\t', file=fileCSV)
    
    errorCount = 0
    with open(bed, "r") as f:
        for line in f:
            cols = line.rstrip().split("\t")
            region = cols[0] + ":" + cols[1] + "-" + cols[2]
            geneName = cols[3]
            tcCount = [0] * maxReadLength
            readCount = 0
            
            
            try:
                refRegion = cols[0] + ":" + str(int(cols[1]) - maxReadLength) + "-" + str(int(cols[2]) + maxReadLength)
                refSeq = reffile.fetch(region=refRegion)
                chromosome = cols[0]
                countFwd = 0
                countRev = 0
                for read in samfile.fetch(region=region):
                    if(read.is_reverse):
                        countRev += 1
                    else:
                        countFwd += 1
                    #chr = samfile.getrname(read.reference_id)
                    tcOnly,agOnly = getTC(read, refSeq, chromosome, refRegion, cols[1], maxReadLength, minQual)
                    #Make strand specific in the future???
                    tc = tcOnly + agOnly
                    tcCount[tc] += 1
                    readCount += 1
                
                snpInUTR = 0
                finalCount = 0;
                if(countRev > countFwd):
                    snpInUTR = getSNPsInUTR(chromosome, int(cols[1]), int(cols[2]), 2)
                else:
                    snpInUTR = getSNPsInUTR(chromosome, int(cols[1]), int(cols[2]), 1)
                
                percTC = 0
                if(readCount > 0):
                    percTC = ((readCount - tcCount[0]) * 1.0 / readCount)
                #chr    start    stop    read without T->C    reads with T->C    Percentage of read with T->C    Number of forward reads mapping to region    Number of reverse reads mapping to region    T->C SNPs found in region
                print(cols[0], cols[1], cols[2], geneName, tcCount[0], tcCount[0] / readNumber * 1000000, (readCount - tcCount[0]), (readCount - tcCount[0]) / readNumber * 1000000, percTC, countFwd, countRev, snpInUTR, sep='\t', file=fileCSV)
            except ValueError:
    #             print("", file=sys.stderr)
                errorCount += 1
                print("Error counting TC reads for region " + region, file=log)
                print("Msg: " + str(sys.exc_info()[0]), file=log)
            lineCount += 1
#             if(lineCount % 10 == 0):
#                 sys.stderr.write("\r%d" % lineCount)
#                 sys.stderr.flush()
    
    samfile.close()
    reffile.close()
    
    fileCSV.close()
    
    print("Coudln't count T->C for " + str(errorCount) + " regions. ", file=log)


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

    