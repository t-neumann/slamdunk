#!/usr/bin/env python

from __future__ import print_function
import os
import tempfile
import math
import pysam

from os.path import basename
from slamdunk.utils.misc import run, removeExtension, checkStep, getReadCount, matchFile, complement , getPlotter, callR
from slamdunk.slamseq.SlamSeqFile import SlamSeqBamFile, ReadDirection
from slamdunk.utils import SNPtools
from slamdunk.utils.BedReader import BedIterator

# projectPath = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
# pathComputeOverallRates = os.path.join(projectPath, "plot", "compute_overall_rates.R")
# pathComputeGlobalRates = os.path.join(projectPath, "plot", "globalRatePlotter.R")
# pathComputeTCContext = os.path.join(projectPath, "plot", "compute_context_TC_rates.R")
# pathConversionPerReadPos = os.path.join(projectPath, "plot", "conversion_per_read_position.R")
# pathSampleComparison = os.path.join(projectPath, "plot", "compute_sample_comparison_statistics.R")
# pathComputeHalfLifes = os.path.join(projectPath, "plot", "compute_halflifes.R")

utrNormFactor = 200
baseNumber = 5
toBase = [ 'A', 'C', 'G', 'T', 'N' ]

def sumLists(a, b):
    return [int(x) + int(y) for x, y in zip(a, b)]

def maxLists(a, b):
    return [max(int(x), int(y)) for x, y in zip(a, b)]

def normalizePos(pos, length, factor):
    # return int(round(float(pos) / float(length) * factor))
    return int(math.floor(float(pos) / float(length) * factor))

# Print rates in correct format for plotting
def printRates(ratesFwd, ratesRev, f):
    print("\t", end='', file=f)
    for i in range(0, 5):
        print(toBase[i] + "\t" + toBase[i].lower() + "\t", end='', file=f)
    print(file=f)
    for i in range(0, 5):
        print(toBase[i] + "\t", end='', file=f)
        for j in range(0, 5):
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
        # Init
        totalRatesFwd = [0] * 25
        totalRatesRev = [0] * 25
        tcCount = [0] * 100
         
        # Go through one chr after the other
        testFile = SlamSeqBamFile(bam, referenceFile, None)
         
        chromosomes = testFile.getChromosomes()
         
        for chromosome in chromosomes:
            readIterator = testFile.readsInChromosome(chromosome)
                 
            for read in readIterator:
                 
                # Compute rates for current read
                rates = read.conversionRates
                # Get T -> C conversions for current read
                tc = read.tcCount
                tcCount[tc] += 1
                 
                # Add rates from read to total rates
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
         
        # Print rates in correct format for plotting
        fo = open(outputCSV, "w")
        printRates(totalRatesFwd, totalRatesRev, fo)
        fo.close()
     
    if(not checkStep([bam, referenceFile], [outputPDF], force)):
        print("Skipped computing overall rate pdfs for file " + bam, file=log)
    else:
        f = tempfile.NamedTemporaryFile(delete=False)
        print(removeExtension(basename(bam)), outputCSV, sep='\t', file=f)
        f.close()
             
        #run(pathComputeOverallRates + " -f " + f.name + " -O " + outputPDF, log, dry=printOnly, verbose=verbose)
        callR(getPlotter("compute_overall_rates") + " -f " + f.name + " -O " + outputPDF, log, dry=printOnly, verbose=verbose)

# This is for TC conversions
    
# def statsComputeTCContext(referenceFile, bam, minQual, outputCSV, outputPDF, log, printOnly=False, verbose=True, force=False):
#     
#     if(not checkStep([bam, referenceFile], [outputCSV], force)):
#         print("Skipped computing overall rates for file " + bam, file=log)
#     else:
#         #Init
#         #combinations = ["AT","CT","GT","TT","NT","AA","CA","GA","TA","NA"]
#         frontCombinations = ["AT","CT","GT","TT","NT"]
#         backCombinations = ["TA","TC","TG","TT","TN"]
#         
#         counts = {}
#         counts['5prime']= {}
#         counts['3prime']= {}
#         counts['5prime']['fwd'] = {}
#         counts['5prime']['rev'] = {}
#         counts['3prime']['fwd'] = {}
#         counts['3prime']['rev'] = {}
#         
#         for combination in frontCombinations :
#             counts['5prime']['fwd'][combination] = 0
#             counts['5prime']['rev'][combination] = 0
#             
#         for combination in backCombinations:
#             counts['3prime']['fwd'][combination] = 0
#             counts['3prime']['rev'][combination] = 0
#         
#         #Go through one chr after the other
#         testFile = SlamSeqBamFile(bam, referenceFile, None)
#         
#         chromosomes = testFile.getChromosomes()
#         
#         for chromosome in chromosomes:
#             readIterator = testFile.readsInChromosome(chromosome)
#                 
#             for read in readIterator:
#                 
#                 # Iterator over mismatches for T -> C conversion contexts
#                 # Get mismatches rates for current read
#                 for mismatch in read.mismatches:
#                     if(mismatch.isTCMismatch(read.direction == ReadDirection.Reverse)):
#                         frontContext = mismatch.referenceBase5PrimeContext + mismatch.referenceBase
#                         backContext = mismatch.referenceBase + mismatch.referenceBase3PrimeContext
#                         
#                         if (read.direction == ReadDirection.Reverse) :
#                             frontContext = complement(frontContext)
#                             backContext = complement(backContext)
#                             counts['5prime']['rev'][frontContext] += 1
#                             counts['3prime']['rev'][backContext] += 1
#                         else :
#                             counts['5prime']['fwd'][frontContext] += 1
#                             counts['3prime']['fwd'][backContext] += 1
#         
#         #Print rates in correct format for plotting
#         fo = open(outputCSV, "w")
#         
#         print("\t".join(frontCombinations),file=fo)
#         
#         frontFwdLine = ""
#         frontRevLine = ""
#         backFwdLine = ""
#         backRevLine = ""
#         
#         for combination in frontCombinations :
#             frontFwdLine += str(counts['5prime']['fwd'][combination]) + "\t"
#             frontRevLine += str(counts['5prime']['rev'][combination]) + "\t"
#         
#         print(frontFwdLine.rstrip(),file=fo)
#         print(frontRevLine.rstrip(),file=fo)
#         
#         print("\t".join(backCombinations),file=fo)
# 
#         for combination in backCombinations :
#             backFwdLine += str(counts['3prime']['fwd'][combination]) + "\t"
#             backRevLine += str(counts['3prime']['rev'][combination]) + "\t"
# 
#         print(backFwdLine.rstrip(),file=fo)
#         print(backRevLine.rstrip(),file=fo)
#         
#         fo.close()
#     
#     if(not checkStep([bam, referenceFile], [outputPDF], force)):
#         print("Skipped computing overall rate pdfs for file " + bam, file=log)
#     else:
#         f = tempfile.NamedTemporaryFile(delete=False)
#         print(removeExtension(basename(bam)), outputCSV, sep='\t', file=f)
#         f.close()
#         
#         run(pathComputeTCContext + " -f " + f.name + " -O " + outputPDF, log, dry=printOnly, verbose=verbose)
 
# This is for Ts in reads

def statsComputeTCContext(referenceFile, bam, minQual, outputCSV, outputPDF, log, printOnly=False, verbose=True, force=False):
     
    if(not checkStep([bam, referenceFile], [outputCSV], force)):
        print("Skipped computing overall rates for file " + bam, file=log)
    else:
        # Init
        # combinations = ["AT","CT","GT","TT","NT","AA","CA","GA","TA","NA"]
        frontCombinations = ["AT", "CT", "GT", "TT", "NT"]
        backCombinations = ["TA", "TC", "TG", "TT", "TN"]
         
        counts = {}
        counts['5prime'] = {}
        counts['3prime'] = {}
        counts['5prime']['fwd'] = {}
        counts['5prime']['rev'] = {}
        counts['3prime']['fwd'] = {}
        counts['3prime']['rev'] = {}
         
        for combination in frontCombinations :
            counts['5prime']['fwd'][combination] = 0
            counts['5prime']['rev'][combination] = 0
             
        for combination in backCombinations:
            counts['3prime']['fwd'][combination] = 0
            counts['3prime']['rev'][combination] = 0
             
        bamFile = pysam.AlignmentFile(bam, "rb")
         
        # Go through one chr after the other
        testFile = SlamSeqBamFile(bam, referenceFile, None)
         
        chromosomes = testFile.getChromosomes()
         
        for chromosome in chromosomes:
                 
            for read in bamFile.fetch(region=chromosome):
                 
                i = 0
                while i < len(read.query_sequence):
                    if(read.query_sequence[i] == "T" and not read.is_reverse) :
                        frontContext = None
                        backContext = None
                        if (i > 0) :
                            frontContext = read.query_sequence[i - 1]
                        if (i < (len(read.query_sequence) - 1)) :
                            backContext  = read.query_sequence[i + 1]
                         
                        if (frontContext != None) :
                            counts['5prime']['fwd'][frontContext + "T"] += 1
                        if (backContext != None) :
                            counts['3prime']['fwd']["T" + backContext] += 1
                             
                    if(read.query_sequence[i] == "A" and read.is_reverse) :
                        frontContext = None
                        backContext = None
                        if (i > 0) :
                            backContext = read.query_sequence[i - 1]
                        if (i < (len(read.query_sequence) - 1)) :
                            frontContext  = read.query_sequence[i + 1]
                         
                        if (frontContext != None) :
                            counts['5prime']['rev'][complement(frontContext + "A")] += 1
                        if (backContext != None) :
                            counts['3prime']['rev'][complement("A" + backContext)] += 1
                     
                    i += 1
         
        # Print rates in correct format for plotting
        fo = open(outputCSV, "w")
         
        print("\t".join(frontCombinations), file=fo)
         
        frontFwdLine = ""
        frontRevLine = ""
        backFwdLine = ""
        backRevLine = ""
         
        for combination in frontCombinations :
            frontFwdLine += str(counts['5prime']['fwd'][combination]) + "\t"
            frontRevLine += str(counts['5prime']['rev'][combination]) + "\t"
         
        print(frontFwdLine.rstrip(), file=fo)
        print(frontRevLine.rstrip(), file=fo)
         
        print("\t".join(backCombinations), file=fo)
 
        for combination in backCombinations :
            backFwdLine += str(counts['3prime']['fwd'][combination]) + "\t"
            backRevLine += str(counts['3prime']['rev'][combination]) + "\t"
 
        print(backFwdLine.rstrip(), file=fo)
        print(backRevLine.rstrip(), file=fo)
         
        fo.close()
     
    if(not checkStep([bam, referenceFile], [outputPDF], force)):
        print("Skipped computing overall rate pdfs for file " + bam, file=log)
    else:
        f = tempfile.NamedTemporaryFile(delete=False)
        print(removeExtension(basename(bam)), outputCSV, sep='\t', file=f)
        f.close()
         
        #run(pathComputeTCContext + " -f " + f.name + " -O " + outputPDF, log, dry=printOnly, verbose=verbose)
        callR(getPlotter("compute_context_TC_rates") + " -f " + f.name + " -O " + outputPDF, log, dry=printOnly, verbose=verbose)
        
# def statsComputeTCContext(referenceFile, bam, bed, minQual, outputCSV, outputPDF, log, printOnly=False, verbose=True, force=False):
#     
#     #if(not checkStep([bam, referenceFile], [outputCSV], force)):
#     #    print("Skipped computing overall rates for file " + bam, file=log)
#     #else:
#         # Init
#         # combinations = ["AT","CT","GT","TT","NT","AA","CA","GA","TA","NA"]
#         #frontCombinations = ["AT", "CT", "GT", "TT", "NT"]
#         #backCombinations = ["TA", "TC", "TG", "TT", "TN"]
#         
#         counts = {}
#         counts['5prime'] = {}
#         counts['3prime'] = {}
#         
# #         for combination in frontCombinations :
# #             counts['5prime'][combination] = 0
# #             
# #         for combination in backCombinations:
# #             counts['3prime'][combination] = 0
#         
#         ref = pysam.FastaFile(referenceFile)
#         
#         #chromosomes = list(ref.references)
#          
#         #for chromosome in chromosomes:
#          
#         for utr in BedIterator(bed):
#             
#             chromosome = utr.chromosome + ":" + str(utr.start) + "-" + str(utr.stop)
#             seq = ref.fetch(region=chromosome).upper()
#             
#             #seq = ref.fetch(chromosome,0).upper()
#             #print("Handling chromosome " + chromosome)
#             
# #             i = 0
# #             while i < len(seq):
# #                 if(seq[i] == "T" and i > 0 and i < (len(seq) - 1)) :
# #                         frontContext = seq[i - 1]
# #                         backContext  = seq[i + 1]
# #                     
# #                         counts['5prime'][frontContext + "T"] += 1
# #                         counts['3prime']["T" + backContext] += 1
# #                     
# #                 i += 1
#             for (front, middle, back) in zip(seq[0::1], seq[1::1], seq[2::1]):
#                 triplet = front + middle  + back
#                 #if (middle == "T") :
#                 if (not counts['5prime'].has_key(triplet[0:2])) :
#                     counts['5prime'][triplet[0:2]] = 0
#                 counts['5prime'][triplet[0:2]] += 1
#                 if (not counts['3prime'].has_key(triplet[1:3])) :
#                     counts['3prime'][triplet[1:3]] = 0
#                 counts['3prime'][triplet[1:3]] += 1
#         
#         # Print rates in correct format for plotting
#         fo = open(outputCSV, "w")
#         
#         for baseA in toBase :
#             headerFront = ""
#             headerBack = ""
#             
#             for baseB in toBase :
#                 headerFront += baseB + baseA + "\t"
#                 headerBack += baseA + baseB + "\t"
#                 
#             valuesFront = ""
#             valuesBack = ""
#             
#             for baseB in toBase :
#                 if (counts['5prime'].has_key(baseB + baseA)) :
#                     valuesFront += str(counts['5prime'][baseB + baseA]) + "\t"
#                 else :
#                     valuesFront += "0\t"
#                 if (counts['3prime'].has_key(baseA + baseB)):
#                     valuesBack += str(counts['3prime'][baseA + baseB]) + "\t"
#                 else:
#                     valuesBack += "0\t"
#                 
#             print(headerFront.rstrip(),file=fo)
#             print(valuesFront.rstrip(),file=fo)
#             print(headerBack.rstrip(),file=fo)
#             print(valuesBack.rstrip(),file=fo)
#             
# #         print("\t".join(frontCombinations), file=fo)
# #         
# #         frontFwdLine = ""
# #         backFwdLine = ""
# #         
# #         for combination in frontCombinations :
# #             frontFwdLine += str(counts['5prime'][combination]) + "\t"
# #         
# #         print(frontFwdLine.rstrip(), file=fo)
# #         
# #         print("\t".join(backCombinations), file=fo)
# # 
# #         for combination in backCombinations :
# #             backFwdLine += str(counts['3prime'][combination]) + "\t"
# # 
# #         print(backFwdLine.rstrip(), file=fo)
#         
#         fo.close()
#     
#     #if(not checkStep([bam, referenceFile], [outputPDF], force)):
#     #    print("Skipped computing overall rate pdfs for file " + bam, file=log)
#     #else:
# #         f = tempfile.NamedTemporaryFile(delete=False)
# #         print(removeExtension(basename(bam)), outputCSV, sep='\t', file=f)
# #         f.close()
# #         
# #         run(pathComputeTCContext + " -f " + f.name + " -O " + outputPDF, log, dry=printOnly, verbose=verbose)

def statsComputeOverallRatesPerUTR(referenceFile, bam, minQual, outputCSV, outputPDF, utrBed, maxReadLength, log, printOnly=False, verbose=True, force=False):
    
    if(not checkStep([bam, referenceFile], [outputCSV], force)):
        print("Skipped computing overall rates for file " + bam, file=log)
    else:
        # Go through one chr after the other
        testFile = SlamSeqBamFile(bam, referenceFile, None)
        
        fo = open(outputCSV, "w")
        print("Name", "Chr", "Start", "End", "Strand", "ReadCount", sep="\t", end="\t", file=fo)
        for i in range(0, 5):
            for j in range(0, 5):
                print(toBase[i].upper() + "_" + toBase[j].upper(), end="", file=fo)
                if(i != 4 or j != 4):
                    print("\t", end="", file=fo)
        print(file=fo)
                        
        for utr in BedIterator(utrBed):
                                         
            readIterator = testFile.readInRegion(utr.chromosome, utr.start, utr.stop, utr.strand, maxReadLength)
            
            # Init
            totalRates = [0] * 25
            
            readCount = 0
            for read in readIterator:
                
                # Compute rates for current read
                rates = read.conversionRates
                
                # Add rates from read to total rates
                totalRates = sumLists(totalRates, rates)
                readCount += 1
                    
            print(utr.name, utr.chromosome, utr.start, utr.stop, utr.strand, readCount, "\t".join(str(x) for x in totalRates), sep="\t", file=fo)
        fo.close()
        
    if(not checkStep([bam, referenceFile], [outputPDF], force)):
        print("Skipped computing global rate pdfs for file " + bam, file=log)
    else:
        f = tempfile.NamedTemporaryFile(delete=False)
        print(removeExtension(basename(bam)), outputCSV, sep='\t', file=f)
        f.close()
              
        #run(pathComputeGlobalRates + " -f " + f.name + " -O " + outputPDF, log, dry=printOnly, verbose=verbose)
        callR(getPlotter("globalRatePlotter") + " -f " + f.name + " -O " + outputPDF, log, dry=printOnly, verbose=verbose)
           
    
def readSummary(mappedFiles, filteredFiles, dedupFiles, snpsFiles, samples, outputPrefix, log, printOnly=False, verbose=True, force=False):
    if(len(mappedFiles) == len(snpsFiles)):
        outputFile = open(outputPrefix + "_summary.txt", "w")
                
#         print("Filename", "Name",  "Sequenced reads", "Mapped reads", "Filtered reads", "SNP count", "T->C SNP count", sep=";", file=outputFile)
        header = ";".join(["Filename", "Name", "Sequenced reads", "Mapped reads"])
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
                print(sample, name, mappedStats.TotalReads, mappedStats.MappedReads, file=outputFile, sep=";", end="")
                
                if (filteredFiles != None) :
                    filteredStats = getReadCount(filteredFile)
                    print(";" + str(filteredStats.MappedReads), file=outputFile, end="")
                
                if (dedupFiles != None) :
                    dedupStats = getReadCount(dedupFile)
                    print(";" + str(dedupStats.MappedReads), file=outputFile, end="")
                    
                print("", file=outputFile)
                
#                 snpCount, tcSnpCount = snps.countSNPsInFile(snpFile)
#                 print(sample, name, mappedStats.TotalReads, mappedStats.MappedReads, filteredStats.MappedReads, snpCount, tcSnpCount, file=outputFile, sep=";")
                # print(sample, name, mappedStats.TotalReads, mappedStats.MappedReads, filteredStats.MappedReads, dedupStats.MappedReads, file=outputFile, sep=";")
        outputFile.close()
    else:
        print("Files missing", file=log)
        
def sampleSummary(readCounts, outputPrefix, log, printOnly=False, verbose=True, force=False):
    if(not checkStep([readCounts], [], force)):
        print("Skipped computing pairwise correlation plots for file " + readCounts, file=log)
    else: 
        run(pathSampleComparison + " -i " + readCounts + " -o " + outputPrefix, log, dry=printOnly, verbose=verbose)        
    
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
        snps.read()
        
        # Go through one chr after the other
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
            print(allPerPosFwd[i], allPerPosRev[i], tcPerPosFwd[i], tcPerPosRev[i], totalReadCountFwd[i], totalReadCountRev[i], sep='\t', file=foTC)
        foTC.close()
       
    if(not checkStep([outputCSV], [outputPDF], force)):
        print("Skipped computing T->C per reads position plot for file " + bam, file=log)
    else: 
        #run(pathConversionPerReadPos + " -i " + outputCSV + " -o " + outputPDF, log, dry=printOnly, verbose=verbose)
        callR(getPlotter("conversion_per_read_position") + " -i " + outputCSV + " -o " + outputPDF, log, dry=printOnly, verbose=verbose)
        
def tcPerUtr(referenceFile, utrBed, bam, minQual, maxReadLength, outputCSV, outputPDF, snpsFile, log, printOnly=False, verbose=True, force=False):
        
    if(not checkStep([bam, referenceFile], [outputCSV], force)):
        print("Skipped computing T->C per UTR position for file " + bam, file=log)
    else:
    
        counter = 0
            
        totalUtrCountFwd = [0] * utrNormFactor
        totalUtrCountRev = [0] * utrNormFactor
        
        tcPerPosRev = [0] * utrNormFactor
        tcPerPosFwd = [0] * utrNormFactor
         
        allPerPosRev = [0] * utrNormFactor
        allPerPosFwd = [0] * utrNormFactor
        
        snps = SNPtools.SNPDictionary(snpsFile)
        snps.read()
        
        # Go through one utr after the other
        testFile = SlamSeqBamFile(bam, referenceFile, snps)
        
        for utr in BedIterator(utrBed):
                                         
            readIterator = testFile.readInRegion(utr.chromosome, utr.start, utr.stop, utr.strand, maxReadLength)
            
            tcForwardCounts = [0] * utrNormFactor
            mutForwardCounts = [0] * utrNormFactor
            tcReverseCounts = [0] * utrNormFactor
            mutReverseCounts = [0] * utrNormFactor
            
            for read in readIterator:
                
                tcCounts = [0] * utrNormFactor
                mutCounts = [0] * utrNormFactor
                
                for mismatch in read.mismatches:
                             
                    mismatchPos = mismatch.referencePosition

                    # mismatchPos = read.startRefPos
                        
                    if (utr.strand == "+") :
                                                
                        # New try for UTRs (remove + 1
                        if (mismatchPos >= (utr.getLength() - utrNormFactor) and mismatchPos < utr.getLength()) :
                        # if (mismatchPos >= (utr.getLength() - utrNormFactor) and mismatchPos < utr.getLength() + 1) :
                            mismatchPos = utrNormFactor - (utr.getLength() - mismatchPos)
                            
                            if(mismatch.isTCMismatch(read.direction == ReadDirection.Reverse)):
                                tcCounts[mismatchPos] += 1
                            else :
                                mutCounts[mismatchPos] += 1                    
                    else :
                        
#                         mismatchPos = read.endRefPos
                    
                        if (mismatchPos >= 0 and mismatchPos < min(utr.getLength(), utrNormFactor)) :
                            if(mismatch.isTCMismatch(read.direction == ReadDirection.Reverse)):
                                tcCounts[mismatchPos] += 1
                            else :
                                mutCounts[mismatchPos] += 1
                            
                if(read.direction == ReadDirection.Reverse):
                    
                    tcReverseCounts = sumLists(tcReverseCounts, tcCounts)
                    mutReverseCounts = sumLists(mutReverseCounts, mutCounts)
                    
                    start = max(0, min(min(utr.getLength(), utrNormFactor), read.startRefPos))
                    end = max(0, min(min(utr.getLength(), utrNormFactor), read.endRefPos))
                    
                    for i in range(start, end):
                        
                        totalUtrCountRev[i] += 1
                     
                else:
                            
                    tcForwardCounts = sumLists(tcForwardCounts, tcCounts)
                    mutForwardCounts = sumLists(mutForwardCounts, mutCounts)
                    
                    start = min(utr.getLength(), max(utr.getLength() - utrNormFactor, read.startRefPos))
                    end = min(utr.getLength(), max(utr.getLength() - utrNormFactor, read.endRefPos))
                
                    for i in range(start, end):
                        normPos = utrNormFactor - (utr.getLength() - i)
                        totalUtrCountFwd[normPos] += 1                 
                        
            tcPerPosFwd = sumLists(tcPerPosFwd, tcForwardCounts)
            allPerPosFwd = sumLists(allPerPosFwd, mutForwardCounts)
             
            tcPerPosRev = sumLists(tcPerPosRev, tcReverseCounts)
            allPerPosRev = sumLists(allPerPosRev, mutReverseCounts)
            
            counter += 1
            
            if (verbose and counter % 10000 == 0) :
                print("Handled " + str(counter) + " UTRs.", file=log)
    
        foTC = open(outputCSV, "w")
        
        reverseAllPerPosRev = allPerPosRev[::-1]
#         reverseAllPerPosRev = reverseAllPerPosRev[0:utrNormFactor - 1]
        reverseTcPerPosRev = tcPerPosRev[::-1]
#         reverseTcPerPosRev = reverseTcPerPosRev[0:utrNormFactor - 1]
        reverseTotalUtrCountRev = totalUtrCountRev[::-1]
#         reverseTotalUtrCountRev = reverseTotalUtrCountRev[0:utrNormFactor - 1]
        
        #print(allPerPosFwd)
        #print(tcPerPosFwd)
        #print(totalUtrCountFwd)

#         allPerPosFwd = allPerPosFwd[1:utrNormFactor]
#         tcPerPosFwd = tcPerPosFwd[1:utrNormFactor]
#         totalUtrCountFwd = totalUtrCountFwd[1:utrNormFactor]
        
        #for i in range(0, utrNormFactor - 1):
        for i in range(0, utrNormFactor):
            #print(allPerPosFwd[i], reverseAllPerPosRev[i], tcPerPosFwd[i], reverseTcPerPosRev[i], totalUtrCountFwd[i], reverseTotalUtrCountRev[i], sep='\t', file=foTC)
            print(allPerPosFwd[i], reverseAllPerPosRev[i], tcPerPosFwd[i], reverseTcPerPosRev[i], totalUtrCountFwd[i], reverseTotalUtrCountRev[i], sep='\t', file=foTC)
        foTC.close()
       
    if(not checkStep([outputCSV], [outputPDF], force)):
        print("Skipped computing T->C per UTR position plot for file " + bam, file=log)
    else: 
        #run(pathConversionPerReadPos + " -u -i " + outputCSV + " -o " + outputPDF, log, dry=printOnly, verbose=verbose)
        callR(getPlotter("conversion_per_read_position") + " -u -i " + outputCSV + " -o " + outputPDF, log, dry=printOnly, verbose=verbose)


def halflifes(bams, outputCSV, timepoints, log, printOnly=False, verbose=True, force=False):
    
    #run("Rscript " + pathComputeHalfLifes + " -f " + bams + " -t " + timepoints + " -o " + outputCSV, log, dry=printOnly, verbose=verbose)
    callR(getPlotter("compute_halflifes") + " -f " + bams + " -t " + timepoints + " -o " + outputCSV, log, dry=printOnly, verbose=verbose)

