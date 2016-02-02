#!/usr/bin/env python

from __future__ import print_function
import os
import tempfile

from os.path import basename
from utils.misc import run
from utils.misc import removeExtension, checkStep, getReadCount, matchFile
from slamseq.SlamSeqFile import SlamSeqFile, ReadDirection


projectPath = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
pathComputeOverallRates = os.path.join(projectPath, "plot", "compute_overall_rates.R")

baseNumber = 5
toBase = [ 'A', 'C', 'G', 'T', 'N' ]


def sumLists(a, b):
    return [int(x) + int(y) for x, y in zip(a, b)]

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

def statsComputeOverallRates(referenceFile, bam, minQual, outputCSV, outputPDF, log, printOnly=False, verbose=True, force=False):
    
    if(not checkStep([bam, referenceFile], [outputCSV], force)):
        print("Skipped computing overall rates for file " + bam, file=log)
    else:
        #Init
        totalRatesFwd = [0] * 25
        totalRatesRev = [0] * 25
        tcCount = [0] * 100
        
        #Go through one chr after the other
        testFile = SlamSeqFile(bam, referenceFile, None)
        
        chromosomes = testFile.getChromosomes()
        
        for chromosome in chromosomes:
            readIterator = testFile.readsInChromosome(chromosome)
                
            for read in readIterator:
                
                #Compute rates for current read
                rates = read.conversionRates
                #Get T -> C conversions for current read
                tc = read.tcCount
                tcCount[tc] += 1
                
                #Add rates from read to total rates
                if(read.direction == ReadDirection.Reverse):
                    totalRatesRev = sumLists(totalRatesRev, rates)
                else:
                    totalRatesFwd = sumLists(totalRatesFwd, rates)
        
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
    
    if(not checkStep([bam, referenceFile], [outputPDF], force)):
        print("Skipped computing overall rate pdfs for file " + bam, file=log)
    else:
        f = tempfile.NamedTemporaryFile(delete=False)
        print(removeExtension(basename(bam)), outputCSV, sep='\t', file=f)
        f.close()
            
        run(pathComputeOverallRates + " -f " + f.name + " -O " + outputPDF, log, dry=printOnly, verbose=verbose)
            
    
def readSummary(mappedFiles, filteredFiles, snpsFiles, samples, ouputCSV, log, printOnly=False, verbose=True, force=False):
    if(len(mappedFiles) == len(filteredFiles) and len(filteredFiles) == len(snpsFiles)):
        outputFile = open(ouputCSV, "w")
#         print("Filename", "Name",  "Sequenced reads", "Mapped reads", "Filtered reads", "SNP count", "T->C SNP count", sep=";", file=outputFile)
        print("Filename", "Name",  "Sequenced reads", "Mapped reads", "Filtered reads", sep=";", file=outputFile)
        for sample in samples:
            name = samples[sample]
            mappedFile = matchFile(sample, mappedFiles)
            filteredFile = matchFile(sample, filteredFiles)
            snpFile = matchFile(sample, snpsFiles)
            if(snpFile == None or filteredFile == None or mappedFile == None):
                raise RuntimeError("Couldn't match all files.")
            else:
                mappedStats = getReadCount(mappedFile)
                filteredStats = getReadCount(filteredFile)
#                 snpCount, tcSnpCount = snps.countSNPsInFile(snpFile)
#                 print(sample, name, mappedStats.TotalReads, mappedStats.MappedReads, filteredStats.MappedReads, snpCount, tcSnpCount, file=outputFile, sep=";")
                print(sample, name, mappedStats.TotalReads, mappedStats.MappedReads, filteredStats.MappedReads, file=outputFile, sep=";")
        outputFile.close()
    else:
        print("Files missing", file=log)
        

    
