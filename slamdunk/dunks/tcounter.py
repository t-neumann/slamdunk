#!/usr/bin/env python

# Date located in: -
from __future__ import print_function
import sys, os

import pysam

# from argparse import ArgumentParser, RawDescriptionHelpFormatter
# 
# usage = ""
# parser = ArgumentParser(description=usage, formatter_class=RawDescriptionHelpFormatter)
# parser.add_argument("-i", "--input", type=str, required=True, dest="bam", help="BAM file" )
# parser.add_argument("-r", "--reference", type=str, required=True, dest="ref", help="Reference fasta file" )
# parser.add_argument("-b", "--bed", type=str, required=True, dest="bed", help="BED file" )
# parser.add_argument("-s", "--snps", type=str, required=False, dest="snp", help="SNP file" )
# parser.add_argument("-l", "--max-read-length", type=int, required=True, dest="maxLength", help="Max read length in BAM file" )
# parser.add_argument("-q", "--min-base-qual", type=int, default=0, required=False, dest="minQual", help="Min base quality for T -> C conversions" )
# args = parser.parse_args()

snpDict = {}

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


def count(ref, bed, snps, bam, maxReadLength, minQual, outputCSV):
    
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
    print("chr", "start", "end", "gene_name", "non_tc_read_count", "tc_read_count", "tc_read_perc", "fwd_reads", "rev_reads", "snp_In_UTR", sep='\t', file=fileCSV)
    
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
                print(cols[0], cols[1], cols[2], geneName, tcCount[0], (readCount - tcCount[0]), percTC, countFwd, countRev, snpInUTR, sep='\t', file=fileCSV)
            except ValueError:
    #             print("", file=sys.stderr)
                errorCount += 1
                #print("Error counting TC reads for region " + region, file=sys.stderr)
                #print("Msg: " + str(sys.exc_info()[0]), file=sys.stderr)
            lineCount += 1
#             if(lineCount % 10 == 0):
#                 sys.stderr.write("\r%d" % lineCount)
#                 sys.stderr.flush()
    
    samfile.close()
    reffile.close()
    
    fileCSV.close()
    
    print("Coudln't count T->C for " + str(errorCount) + " regions. ", file=sys.stderr)
