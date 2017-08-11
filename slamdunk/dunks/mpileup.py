#!/usr/bin/env python

from __future__ import print_function
import sys, os

from argparse import ArgumentParser, RawDescriptionHelpFormatter
    
import re, os, sys
import pysam, random
from collections import defaultdict

usage = "Custom SNP filtering"
version = "0.1.0"
# Main Parsers
parser = ArgumentParser(description=usage, formatter_class=RawDescriptionHelpFormatter, version=version)

parser.add_argument("-p", "--pileup", type=str, required=True, dest="pileup", help="samtools mpileup file")

args = parser.parse_args()

# This is the pileup format
# chr pos ref coverage read bases baseQ
# chr10   3105486 G       1       ^$,     F
# . is match fwd , is match rev ACGTN match fwd, acgtn match rev, +[0-9][AGCTN] ins, -[0-9][AGCTN] del 

# int minCoverage = 8;
#         int minReads2 = 2;
#         int minAvgQual = 15;
#         double minVarFreq = 0.01;
#         double minFreqForHom = 0.75;
#         double pValueThreshold = 0.99;
#         double strandPvalueThreshold = 0.01;
#         boolean variantsOnly = false;
#         boolean snpsOnly = false;
#         boolean indelsOnly = false;
#         boolean strandFilter = true;
#         String sampleList = "";
#                    snpsOnly = true;

#varscanCmd = "java -jar " + getBinary("VarScan.v2.4.1.jar") + " mpileup2snp  --strand-filter 0 --output-vcf --min-var-freq " + str(minVarFreq) + " --min-coverage " + str(minCov) + " --variants 1"

with open(args.pileup) as f:
    for line in f:
        
        fields = line.rstrip().split("\t")
        

#         pile = fields[4]
#          
#         # Replace starts
#         pile = re.sub("\^.","",pile)
#         # Replace matches
#         pile = re.sub("\.","",pile)
#         pile = re.sub(",","",pile)
#         # Replace dels
#         pile = re.sub("\*","",pile)
#         # Replace ends
#         pile = re.sub("\$","",pile)
#          
#         indel = False
#         type = ""
#         length = ""
#         bases = ""
#         parsed = 0
#          
#         if len(pile) == 0:
#             if not "ref" in mutDict[chr][pos]:
#                 mutDict[chr][pos]["ref"] = 0
#          
#         for char in pile:
#             if (char == "+") :
#                 indel = True
#                 type = "ins"
#             elif (char == "-") :
#                 indel = True
#                 type = "del"
#             elif (char.isdigit()) :
#                 length = length + char
#             else :
#                 if indel:
#                     if bases == "":
#                         length = int(length)
#                     bases = bases + char
#                     parsed += 1
#  
#                     if (parsed == length) :
#                          
#                         indelType = type + bases
#                          
#                         indel = False
#                         length = ""
#                         parsed = 0
#                         bases = ""
#                         type = ""
#                          
#                         if not indelType in mutDict[chr][pos]:
#                             mutDict[chr][pos][indelType] = 0
#                         mutDict[chr][pos][indelType] += 1
#                 else :
#                     conversion = ref.upper() + ">" + char.upper()
#                     if not conversion in mutDict[chr][pos]:
#                         mutDict[chr][pos][conversion] = 0
#                     mutDict[chr][pos][conversion] += 1