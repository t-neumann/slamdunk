#!/usr/bin/env python

from __future__ import print_function
import sys, os, re

import pysam
import tempfile

from os.path import basename
from utils import run
from dunks.utils import removeExtension

projectPath = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
pathComputeOverallRates = os.path.join(projectPath, "plot", "compute_overall_rates.R")

baseNumber = 5
toBase = [ 'A', 'C', 'G', 'T', 'N' ]

def atoi(text):
    return int(text) if text.isdigit() else text

def natural_keys(text):
    '''
    alist.sort(key=natural_keys) sorts in human order
    http://nedbatchelder.com/blog/200712/human_sorting.html
    (See Toothy's implementation in the comments)
    '''
    return [ atoi(c) for c in re.split('(\d+)', text) ]

def encodeBase(base):
    if(base.upper() == 'A'):
        return 0
    if(base.upper() == 'C'):
        return 1
    if(base.upper() == 'G'):
        return 2
    if(base.upper() == 'T'):
        return 3
    
    return 4

def incRate(rates, refBase, readBase):
    rates[baseNumber * encodeBase(refBase) + encodeBase(readBase)] += 1
    return rates

def getRate(rates, refBase, readBase):
    return rates[baseNumber * encodeBase(refBase) + encodeBase(readBase)]
    
def getTCNgm(read):
    return int(read.get_tag("TC"))

def computeRatesForRead(read, refSeq, minQual):
    rates = [0] * 25

    for pair in read.get_aligned_pairs(matches_only=True):
        readPos = pair[0]
        refPos = pair[1]
        refBase = refSeq[refPos]
        readBase = read.query_sequence[readPos]
        readQlty = read.query_qualities[readPos]
        
        if(readQlty >= minQual):
            rates = incRate(rates, refBase, readBase)
                    
    return rates

def sumLists(a, b):
    return [int(x) + int(y) for x, y in zip(a, b)]

#Check if two rates arrays are equal
def compareLists(a, b):
    if(len(a) != len(b)):
        return False
    for x,y in zip(a, b):
        if(a != b):
            return False
    
    return True

#Print rates in correct format for plotting
def printRates(ratesFwd, ratesRev, f):
    print("\t", end='', file=f)
    for i in range(0, 5):
        print(toBase[i] + "\t" + toBase[i].lower() + "\t", end='', file=f)
    print(file=f)
    for i in range(0, 5):
        print(toBase[i] + "\t", end='', file=f)
        for j in range(0,5):
            print(str(ratesFwd[i * 5 + j]) + "\t" + str(ratesRev[i * 5 + j]) + "\t", end='', file=f)
        print(file=f)

#Get the T -> C count from rates array.         
def getTCCount(rev, rates):
    if(not rev):
        return getRate(rates, 'T', 'C')
    else:
        return getRate(rates, 'A', 'G')
            
# ref = args.ref
# # bed = args.bed
# bam = args.bam
# #maxReadLength = args.maxLength
# minQual = args.minQual
# rateFile = args.rateFile
# tcFile = args.tcFile

def statsComputeOverallRates(ref, bam, minQual, outputCSV, outputPDF, printOnly=False, verbose=True, force=True):
    reffile = pysam.FastaFile(ref)
    samfile = pysam.AlignmentFile(bam, "rb")
    
    #Init
    totalRatesFwd = [0] * 25
    totalRatesRev = [0] * 25
    tcCount = [0] * 100
    refCount = 1
    totalRefCount = len(reffile.references)
    
    refs = list(reffile.references)
    refs.sort(key=natural_keys)
    #Go through one chr after the other
    for ref in refs:
        readCount = 0
        
        #Read chr into memory
        refSeq = reffile.fetch(ref)
        #Get total number of reads mapping to chr
        totalCount = samfile.count(ref)
        #Get all reads on chr
        for read in samfile.fetch(ref):
            
            try:
                #Compute rates for current read
                rates = computeRatesForRead(read, refSeq, minQual)
                #Get T -> C conversions for current read
                tc = getTCCount(read.is_reverse, rates)
                tcCount[tc] += 1
                
#                 #If mapped with NGM-slamseq check if results are the same
#                 if(read.has_tag("RA")):
#                     ratesNgm = map(int, read.get_tag("RA").split(","))
#                     tcNgm = getTCCount(read.is_reverse, ratesNgm)
#                     if(not compareLists(rates, ratesNgm) or tc != tcNgm):
#                         print("Difference found:")
#                         print(read)
#                         print(ratesNgm)
#                         print(rates)
#                         print("TC (ngm): " + str(tcNgm))
#                         print("TC (pys): " + str(tc))
#                         #sys.stdin.read(1)
        
                #Add rates from read to total rates
                if(read.is_reverse):
                    totalRatesRev = sumLists(totalRatesRev, rates)
                else:
                    totalRatesFwd = sumLists(totalRatesFwd, rates)
            except:
                print("")
                print("Error computing rates for read " + read.query_name)
                #print("Msg: " + sys.exc_info()[0])
                print(read)
            #Progress info
            readCount += 1
            if(readCount % 1000 == 0):
                sys.stderr.write("\rProcessing %s (%d/%d): %d%%                        " % (ref, refCount, totalRefCount, int(readCount * 100.0 / totalCount)))
                sys.stderr.flush()
                
        refCount += 1
    
    #Cleanup
    sys.stderr.write("\n")
    samfile.close()
    reffile.close()
    
    
#     #Writing T -> C counts to file
#     i = 0;
#     foTC = open(tcFile, "w")
#     for x in tcCount:
#         print(i, x, sep='\t', file=foTC)
#         i += 1
#     foTC.close()
    
    #Print rates in correct format for plotting
    fo = open(outputCSV, "w")
    printRates(totalRatesFwd, totalRatesRev, fo)
    fo.close()
    
    f = tempfile.NamedTemporaryFile(delete=False)
    print(removeExtension(basename(bam)), outputCSV, sep='\t', file=f)
    f.close()
        
    run(pathComputeOverallRates + " -f " + f.name + " -O " + outputPDF, dry=printOnly, verbose=verbose)
        
