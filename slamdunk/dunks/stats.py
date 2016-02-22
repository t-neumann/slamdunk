#!/usr/bin/env python

from __future__ import print_function
import os
import tempfile

from os.path import basename
from utils.misc import run
from utils.misc import removeExtension, checkStep, getReadCount, matchFile
from slamseq.SlamSeqFile import SlamSeqBamFile, ReadDirection
from utils import SNPtools
from utils.BedReader import BedIterator, BedEntry

projectPath = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
pathComputeOverallRates = os.path.join(projectPath, "plot", "compute_overall_rates.R")
pathConversionPerReadPos = os.path.join(projectPath, "plot", "conversion_per_read_position.R")

utrNormFactor = 100
baseNumber = 5
toBase = [ 'A', 'C', 'G', 'T', 'N' ]

def sumLists(a, b):
    return [int(x) + int(y) for x, y in zip(a, b)]

def maxLists(a, b):
    return [max(int(x),int(y)) for x, y in zip(a, b)]

def normalizePos(pos, length, factor):
    return int(round(float(pos) / float(length) * factor))


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
        
# def perUtrTemplate(referenceFile, utrBed, snpsFile, bam, maxReadLength, minQual, outputCSV, log):
# 
#     fileCSV = open(outputCSV,'w')
#     
#     snps = SNPtools.SNPDictionary(snpsFile)
# 
#     #Go through one chr after the other
#     testFile = SlamSeqBamFile(bam, referenceFile, snps)
#                       
#        
#     for utr in BedIterator(utrBed):
#                      
#         readIterator = testFile.readInRegion(utr.chromosome, utr.start, utr.stop, maxReadLength)
#         
#         for read in readIterator:        
#             print(read)
#         
#     fileCSV.close()
  
# def perReadTemplate(referenceFile, bam, minQual, outputCSV, log, printOnly=False, verbose=True, force=False):
#     
#     if(not checkStep([bam, referenceFile], [outputCSV], force)):
#         print("Skipped computing xy for file " + bam, file=log)
#     else:
#         #Go through one chr after the other
#         testFile = SlamSeqBamFile(bam, referenceFile, None)
#         
#         chromosomes = testFile.getChromosomes()
#         
#         for chromosome in chromosomes:
#             readIterator = testFile.readsInChromosome(chromosome)
#                 
#             for read in readIterator:
#                 print(read.n)            

def statsComputeOverallRates(referenceFile, bam, minQual, outputCSV, outputPDF, log, printOnly=False, verbose=True, force=False):
    
    if(not checkStep([bam, referenceFile], [outputCSV], force)):
        print("Skipped computing overall rates for file " + bam, file=log)
    else:
        #Init
        totalRatesFwd = [0] * 25
        totalRatesRev = [0] * 25
        tcCount = [0] * 100
        
        #Go through one chr after the other
        testFile = SlamSeqBamFile(bam, referenceFile, None)
        
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
            
    
def readSummary(mappedFiles, filteredFiles, dedupFiles, snpsFiles, samples, ouputCSV, log, printOnly=False, verbose=True, force=False):
    if(len(mappedFiles) == len(snpsFiles)):
        outputFile = open(ouputCSV, "w")
                
#         print("Filename", "Name",  "Sequenced reads", "Mapped reads", "Filtered reads", "SNP count", "T->C SNP count", sep=";", file=outputFile)
        header = ";".join(["Filename", "Name",  "Sequenced reads", "Mapped reads"])
        if (filteredFiles != None) :
            header = header + ";Filtered reads"
        if (dedupFiles != None) :
            header = header + ";Dedup reads"
        print(header, file=outputFile)
        
        for sample in samples:
            name = samples[sample]
            mappedFile = matchFile(sample, mappedFiles)
            filteredFile = None
            mappedFile == None
            
            if (filteredFiles != None) :
                filteredFile = matchFile(sample, filteredFiles)
            if (dedupFiles != None) :
                dedupFile = matchFile(sample, dedupFiles)

            snpFile = matchFile(sample, snpsFiles)
            if((filteredFiles != None and  filteredFile == None) or (dedupFiles != None and dedupFile == None) or snpFile == None or mappedFile == None):
                raise RuntimeError("Couldn't match all files.")
            else:
                mappedStats = getReadCount(mappedFile)
                print(sample, name, mappedStats.TotalReads, mappedStats.MappedReads, file=outputFile, sep=";",end="")
                
                if (filteredFiles != None) :
                    filteredStats = getReadCount(filteredFile)
                    print(";" + str(filteredStats.MappedReads), file=outputFile, end="")
                
                if (dedupFiles != None) :
                    dedupStats = getReadCount(dedupFile)
                    print(";" + str(dedupStats.MappedReads), file=outputFile, end="")
                    
                print("",file=outputFile)
                
#                 snpCount, tcSnpCount = snps.countSNPsInFile(snpFile)
#                 print(sample, name, mappedStats.TotalReads, mappedStats.MappedReads, filteredStats.MappedReads, snpCount, tcSnpCount, file=outputFile, sep=";")
                #print(sample, name, mappedStats.TotalReads, mappedStats.MappedReads, filteredStats.MappedReads, dedupStats.MappedReads, file=outputFile, sep=";")
        outputFile.close()
    else:
        print("Files missing", file=log)
        
    
def tcPerReadPos(referenceFile, bam, minQual, maxReadLength, outputCSV, outputPDF, snpsFile, log, printOnly=False, verbose=True, force=False):
    
    if(not checkStep([bam, referenceFile], [outputCSV], force)):
        print("Skipped computing T->C per reads position for file " + bam, file=log)
    else:
        
        totalReadCountFwd = [0] * maxReadLength
        totalReadCountRev = [0] * maxReadLength
        
        tcPerPosRev = [0] * maxReadLength
        tcPerPosFwd = [0] * maxReadLength
        
        allPerPosRev = [0] * maxReadLength
        allPerPosFwd = [0] * maxReadLength

        
        snps = SNPtools.SNPDictionary(snpsFile)
        
        #Go through one chr after the other
        testFile = SlamSeqBamFile(bam, referenceFile, snps)
        
        chromosomes = testFile.getChromosomes()
        
        for chromosome in chromosomes:
            readIterator = testFile.readsInChromosome(chromosome)
                
            for read in readIterator:
                
                tcCounts = [0] * maxReadLength
                mutCounts = [0] * maxReadLength
                
                for mismatch in read.mismatches:
                    mutCounts[mismatch.readPosition] += 1   
                    if(mismatch.isTCMismatch(read.direction == ReadDirection.Reverse)):
                        tcCounts[mismatch.readPosition] += 1
                
                query_length = len(read.sequence)
                if(read.direction == ReadDirection.Reverse):
                    tcPerPosRev = sumLists(tcPerPosRev, tcCounts)
                    allPerPosRev = sumLists(allPerPosRev, mutCounts)
                    
                    for i in range(0, query_length):
                        totalReadCountRev[i] += 1
                else:
                    tcPerPosFwd = sumLists(tcPerPosFwd, tcCounts)
                    allPerPosFwd = sumLists(allPerPosFwd, mutCounts)
                    
                    for i in range(0, query_length):
                        totalReadCountFwd[i] += 1
                        

        foTC = open(outputCSV, "w")
        for i in range(0, maxReadLength):
            print(allPerPosFwd[i], allPerPosRev[i], tcPerPosFwd[i], tcPerPosRev[i], totalReadCountFwd[i], totalReadCountRev[i],sep='\t', file=foTC)
        foTC.close()
       
    if(not checkStep([outputCSV], [outputPDF], force)):
        print("Skipped computing T->C per reads position plot for file " + bam, file=log)
    else: 
        run(pathConversionPerReadPos + " -i " + outputCSV + " -o " + outputPDF, log, dry=printOnly, verbose=verbose)
        
def tcPerUtr(referenceFile, utrBed, bam, minQual, maxReadLength, outputCSV, outputPDF, snpsFile, log, printOnly=False, verbose=True, force=False):
        
    if(not checkStep([bam, referenceFile], [outputCSV], force)):
        print("Skipped computing T->C per UTR position for file " + bam, file=log)
    else:
    
        counter = 0
            
        totalUtrCount = [0] * utrNormFactor
        
        tcPerPosRev = [0] * utrNormFactor
        tcPerPosFwd = [0] * utrNormFactor
         
        allPerPosRev = [0] * utrNormFactor
        allPerPosFwd = [0] * utrNormFactor

#         tcPerPos = [0] * utrNormFactor
#         allPerPos = [0] * utrNormFactor
        
        snps = SNPtools.SNPDictionary(snpsFile)
        
        #Go through one utr after the other
        testFile = SlamSeqBamFile(bam, referenceFile, snps)
        
        for utr in BedIterator(utrBed):
                                         
            readIterator = testFile.readInRegion(utr.chromosome, utr.start, utr.stop, maxReadLength)
            
            tcForwardCounts = [0] * utrNormFactor
            mutForwardCounts = [0] * utrNormFactor
            tcReverseCounts = [0] * utrNormFactor
            mutReverseCounts = [0] * utrNormFactor

#             tcForwardCounts =  [0] * utrNormFactor
#             mutForwardCounts = [0] * utrNormFactor
            
            for read in readIterator:
                
                tcCounts = [0] * utrNormFactor
                mutCounts = [0] * utrNormFactor
                
                for mismatch in read.mismatches:
                    if (mismatch.referencePosition >= 0 and mismatch.referencePosition < utr.getLength()) :
                        normPos = normalizePos(mismatch.referencePosition, utr.getLength(), utrNormFactor - 1)
                        mutCounts[normPos] = 1
                    if(mismatch.isTCMismatch(read.direction == ReadDirection.Reverse)):
                        tcCounts[normPos] = 1
                
                if(read.direction == ReadDirection.Reverse):
                    tcReverseCounts = maxLists(tcReverseCounts, tcCounts)
                    mutReverseCounts = maxLists(mutReverseCounts, mutCounts)
                     
                else:
                    tcForwardCounts = maxLists(tcForwardCounts, tcCounts)
                    mutForwardCounts = maxLists(mutForwardCounts, mutCounts)

#                 tcForwardCounts = maxLists(tcForwardCounts, tcCounts)
#                 mutForwardCounts = maxLists(mutForwardCounts, mutCounts)
                    
                        
            tcPerPosFwd = sumLists(tcPerPosFwd, tcForwardCounts)
            allPerPosFwd = sumLists(allPerPosFwd, mutForwardCounts)
             
            tcPerPosRev = sumLists(tcPerPosRev, tcReverseCounts)
            allPerPosRev = sumLists(allPerPosRev, mutReverseCounts)

#             tcPerPos = sumLists(tcPerPos, tcForwardCounts)
#             allPerPos = sumLists(allPerPos, mutForwardCounts)
            
            for i in range(0, utrNormFactor):
                totalUtrCount[i] += 1
                
            counter += 1
            
            if (verbose and counter % 10000 == 0) :
                print("Handled " + str(counter) + " UTRs.",file=log)
    
        foTC = open(outputCSV, "w")
        
        for i in range(0, utrNormFactor):
            print(allPerPosFwd[i], allPerPosRev[i], tcPerPosFwd[i], tcPerPosRev[i], totalUtrCount[i] ,sep='\t', file=foTC)
#             print(allPerPos[i], tcPerPos[i],totalUtrCount[i] ,sep='\t', file=foTC)
        foTC.close()
       
    if(not checkStep([outputCSV], [outputPDF], force)):
        print("Skipped computing T->C per UTR position plot for file " + bam, file=log)
    else: 
        run(pathConversionPerReadPos + " -u -i " + outputCSV + " -o " + outputPDF, log, dry=printOnly, verbose=verbose)
