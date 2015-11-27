#!/usr/bin/env python

# Date located in: -
from __future__ import print_function
import sys, os

import pysam

from argparse import ArgumentParser, RawDescriptionHelpFormatter
    
from os.path import basename

usage = ""
parser = ArgumentParser(description=usage, formatter_class=RawDescriptionHelpFormatter)
parser.add_argument("-i", "--input", type=str, required=True, dest="bam", help="BAM file" )
args = parser.parse_args()

samfile = pysam.AlignmentFile(args.bam, "rb")
outfile = pysam.AlignmentFile("-", "wb", template=samfile)

prevChr = ""
prevStart = ""

duplicateBuffer = {}

for read in samfile:
    
    flag = read.cigarstring
    chr = read.reference_id
    start = read.reference_start
    seq = read.query_sequence
    #qual = read.query_qualities.tostring()
        
    if (chr != prevChr or start != prevStart) :
                    
        if (prevChr != "") :
            for curSeq in duplicateBuffer :
                for curFlag in duplicateBuffer[curSeq]:
                    outfile.write(duplicateBuffer[curSeq][curFlag])
                        #i = 1
            #print("Buffer cleared")
            duplicateBuffer.clear()
#             
#     if (start == 122606245):
#             print(read)
                  
    #duplicateBuffer[qual] = {}
    #duplicateBuffer[qual][flag] = {}
    #duplicateBuffer[qual][flag][seq] = read
    
    duplicateBuffer[seq] = {}
    duplicateBuffer[seq][flag] = read
    
#     if (start == 122606245):
#         print(duplicateBuffer)
#         sys.stdin.readline()
     
    prevChr = chr
    prevStart = start
    
for seq in duplicateBuffer:
    for flag in duplicateBuffer[seq] :
        outfile.write(duplicateBuffer[seq][flag])
duplicateBuffer.clear()