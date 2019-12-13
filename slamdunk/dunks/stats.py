#!/usr/bin/env python

# Copyright (c) 2015 Tobias Neumann, Philipp Rescheneder.
#
# This file is part of Slamdunk.
#
# Slamdunk is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as
# published by the Free Software Foundation, either version 3 of the
# License, or (at your option) any later version.
#
# Slamdunk is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Affero General Public License for more details.
#
# You should have received a copy of the GNU Affero General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

from __future__ import print_function
import tempfile
import math
import pysam
import os
import glob
import numpy as np

from slamdunk.utils.misc import removeExtension, replaceExtension, checkStep, getSampleInfo, complement , getPlotter, callR, SlamSeqInfo  # @UnresolvedImport
from slamdunk.slamseq.SlamSeqFile import SlamSeqBamFile, ReadDirection  # @UnresolvedImport
from slamdunk.utils import SNPtools  # @UnresolvedImport
from slamdunk.utils.BedReader import BedIterator  # @UnresolvedImport

from slamdunk.version import __version__  # @UnresolvedImport

utrNormFactor = 250
baseNumber = 5
toBase = [ 'A', 'C', 'G', 'T', 'N' ]

def sumLists(a, b):
    return [int(x) + int(y) for x, y in zip(a, b)]

def sumCounts(countFile, column="ReadCount"):

    columnId = -1

    sum = 0

    with open(countFile) as f:
        for line in f:
            if (not line.startswith("#")):

                columns = line.rstrip().split("\t")

                # Find column
                if line.startswith("Chromosome") :

                    id = 0
                    for col in columns:
                        if col == column:
                            columnId = id
                        id += 1

                    # Column not found
                    if (columnId < 0) :
                        return 0

                else :

                    sum += int(columns[columnId])

    return sum

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

def statsComputeOverallRates(referenceFile, bam, minBaseQual, outputCSV, outputPDF, log, printOnly=False, verbose=True, force=False):

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
            readIterator = testFile.readsInChromosome(chromosome, minBaseQual)

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

        # Print rates in correct format for plotting
        fo = open(outputCSV, "w")
        print("# slamdunk rates v" + __version__, file=fo)
        printRates(totalRatesFwd, totalRatesRev, fo)
        fo.close()

    if(not checkStep([bam, referenceFile], [outputPDF], force)):
        print("Skipped computing overall rate pdfs for file " + bam, file=log)
    else:

        #f = tempfile.NamedTemporaryFile(delete=False)
        #print(removeExtension(basename(bam)), outputCSV, sep='\t', file=f)
        #f.close()

        callR(getPlotter("compute_overall_rates") + " -f " + outputCSV + " -n " + removeExtension(os.path.basename(bam)) + " -O " + outputPDF, log, dry=printOnly, verbose=verbose)

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

def statsComputeTCContext(referenceFile, bam, minBaseQual, outputCSV, outputPDF, log, printOnly=False, verbose=True, force=False):

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
        f = tempfile.NamedTemporaryFile(mode='w',delete=False)
        print(removeExtension(os.path.basename(bam)), outputCSV, sep='\t', file=f)
        f.close()

        callR(getPlotter("compute_context_TC_rates") + " -f " + f.name + " -O " + outputPDF, log, dry=printOnly, verbose=verbose)


def statsComputeOverallRatesPerUTR(referenceFile, bam, minBaseQual, strictTCs, outputCSV, outputPDF, utrBed, maxReadLength, log, printOnly=False, verbose=True, force=False):

    sampleInfo = getSampleInfo(bam)

    slamseqInfo = SlamSeqInfo(bam)

    if(not checkStep([bam, referenceFile], [outputCSV], force)):
        print("Skipped computing overall rates for file " + bam, file=log)
    else:

        # Go through one chr after the other
        testFile = SlamSeqBamFile(bam, referenceFile, None)

        # UTR stats for MultiQC
        utrStats = dict()

        plotConversions = ['A>T', 'A>G', 'A>C',
                           'C>A', 'C>G', 'C>T',
                           'G>A', 'G>C', 'G>T',
                           'T>A', 'T>G', 'T>C',
        ]

        for conversion in plotConversions:
            utrStats[conversion] = list()

        f = tempfile.NamedTemporaryFile(mode='w', delete=False)

        for utr in BedIterator(utrBed):

            readIterator = testFile.readInRegion(utr.chromosome, utr.start, utr.stop, utr.strand, maxReadLength, minBaseQual)

            # Init
            totalRates = [0] * 25

            readCount = 0
            for read in readIterator:

                if (not read.isTcRead and strictTCs and read.tcCount > 0) :
                    pass
                else :

                    # Compute rates for current read
                    rates = read.conversionRates

                    # Add rates from read to total rates
                    totalRates = sumLists(totalRates, rates)
                    readCount += 1

            print(utr.name, utr.chromosome, utr.start, utr.stop, utr.strand, readCount, "\t".join(str(x) for x in totalRates), sep="\t", file=f)

            # Process rates for MultiQC
            # Copied directly, too lazy to do it properly now

            utrDict = {}

            conversionSum = 0

            A_A = totalRates[0]
            conversionSum =+ A_A
            A_C = totalRates[1]
            conversionSum =+ A_C
            A_G = totalRates[2]
            conversionSum =+ A_G
            A_T = totalRates[3]
            conversionSum =+ A_T

            C_A = totalRates[5]
            conversionSum =+ C_A
            C_C = totalRates[6]
            conversionSum =+ C_C
            C_G = totalRates[7]
            conversionSum =+ C_G
            C_T = totalRates[8]
            conversionSum =+ C_T

            G_A = totalRates[10]
            conversionSum =+ G_A
            G_C = totalRates[11]
            conversionSum =+ G_C
            G_G = totalRates[12]
            conversionSum =+ G_G
            G_T = totalRates[13]
            conversionSum =+ G_T

            T_A = totalRates[15]
            conversionSum =+ T_A
            T_C = totalRates[16]
            conversionSum =+ T_C
            T_G = totalRates[17]
            conversionSum =+ T_G
            T_T = totalRates[18]
            conversionSum =+ T_T

            if utr.strand == "-":

                A_A, T_T = T_T,A_A
                G_G, C_C = C_C,G_G
                A_C, T_G = T_G, A_C
                A_G, T_C = T_C, A_G
                A_T, T_A = T_A, A_T
                C_A, G_T = G_T, C_A
                C_G, G_C = G_C, C_G
                C_T, G_A = G_A, C_T

            if conversionSum > 0:

                Asum = A_A + A_C + A_G + A_T
                Csum = C_A + C_C + C_G + C_T
                Gsum = G_A + G_C + G_G + G_T
                Tsum = T_A + T_C + T_G + T_T

                if Asum > 0 :
                    A_T = A_T / float(Asum) * 100
                    A_G = A_G / float(Asum) * 100
                    A_C = A_C / float(Asum) * 100
                else :
                    A_T = 0
                    A_G = 0
                    A_C = 0
                if Csum > 0:
                    C_A = C_A / float(Csum) * 100
                    C_G = C_G / float(Csum) * 100
                    C_T = C_T / float(Csum) * 100
                else :
                    C_A = 0
                    C_G = 0
                    C_T = 0
                if Gsum > 0:
                    G_A = G_A / float(Gsum) * 100
                    G_C = G_C / float(Gsum) * 100
                    G_T = G_T / float(Gsum) * 100
                else :
                    G_A = 0
                    G_C = 0
                    G_T = 0
                if Tsum > 0:
                    T_A = T_A / float(Tsum) * 100
                    T_G = T_G / float(Tsum) * 100
                    T_C = T_C / float(Tsum) * 100
                else :
                    T_A = 0
                    T_G = 0
                    T_C = 0

                utrStats['A>T'].append(A_T)
                utrStats['A>G'].append(A_G)
                utrStats['A>C'].append(A_C)

                utrStats['C>A'].append(C_A)
                utrStats['C>G'].append(C_G)
                utrStats['C>T'].append(C_T)

                utrStats['G>A'].append(G_A)
                utrStats['G>T'].append(G_T)
                utrStats['G>C'].append(G_C)

                utrStats['T>A'].append(T_A)
                utrStats['T>G'].append(T_G)
                utrStats['T>C'].append(T_C)

        f.close()

        fo = open(outputCSV, "w")

        print("# slamdunk utrrates v" + __version__, file=fo)

        print("# Median-Conversions=",end="",file=fo)

        first = True
        for conversion in plotConversions:
            if (not first) :
                print(',',file=fo, end="")
            else :
                first = False
            print(conversion + ":" + str(np.median(utrStats[conversion])),file=fo, end="")
        print(file=fo)

        print("Name", "Chr", "Start", "End", "Strand", "ReadCount", sep="\t", end="\t", file=fo)
        for i in range(0, 5):
            for j in range(0, 5):
                print(toBase[i].upper() + "_" + toBase[j].upper(), end="", file=fo)
                if(i != 4 or j != 4):
                    print("\t", end="", file=fo)
        print(file=fo)

        with open(f.name, "r") as valueFile:
            fo.write(valueFile.read())

        fo.close()

    if(not checkStep([bam, referenceFile], [outputPDF], force)):
        print("Skipped computing global rate pdfs for file " + bam, file=log)
    else:
        f = tempfile.NamedTemporaryFile(mode='w',delete=False)
        print(sampleInfo.Name, outputCSV, sep='\t', file=f)
        f.close()

        callR(getPlotter("globalRatePlotter") + " -f " + f.name + " -O " + outputPDF, log, dry=printOnly, verbose=verbose)


def readSummary(filteredFiles, countDirectory, outputFile, log, printOnly=False, verbose=True, force=False):

    # Print sort by ID
    contentDict = {}

    tsvFile = open(outputFile, "w")

    print("# slamdunk summary v" + __version__, file=tsvFile)

    if (countDirectory != None) :
        f = tempfile.NamedTemporaryFile(mode='w', delete=False)

    for bam in filteredFiles:
        slamseqInfo = SlamSeqInfo(bam)
        sampleInfo = getSampleInfo(bam)

        if (countDirectory != None) :

            countedReads = 0

            countFile = os.path.join(countDirectory, replaceExtension(os.path.basename(bam), ".tsv", "_tcount"))
            if not os.path.exists(countFile):
                print("TCount directory does not seem to contain tcount file for:\t" + countFile)
            else :
                print(sampleInfo.Name, countFile, sep='\t', file=f)
                countedReads = sumCounts(countFile)

            if(int(sampleInfo.ID) in contentDict):
                ID = len(contentDict) + 1
            else:
                ID = sampleInfo.ID

            contentDict[int(ID)] = "\t".join([bam, sampleInfo.Name, sampleInfo.Type, sampleInfo.Time, str(slamseqInfo.SequencedReads), str(slamseqInfo.MappedReads), str(slamseqInfo.DedupReads), str(slamseqInfo.MQFilteredReads), str(slamseqInfo.IdFilteredReads), str(slamseqInfo.NmFilteredReads), str(slamseqInfo.MultimapperReads), str(slamseqInfo.FilteredReads), str(countedReads), slamseqInfo.AnnotationName])

        else :

            if(int(sampleInfo.ID) in contentDict):
                ID = len(contentDict) + 1
            else:
                ID = sampleInfo.ID

            contentDict[int(ID)] = "\t".join([bam, sampleInfo.Name, sampleInfo.Type, sampleInfo.Time, str(slamseqInfo.SequencedReads), str(slamseqInfo.MappedReads), str(slamseqInfo.DedupReads), str(slamseqInfo.MQFilteredReads), str(slamseqInfo.IdFilteredReads), str(slamseqInfo.NmFilteredReads), str(slamseqInfo.MultimapperReads), str(slamseqInfo.FilteredReads), slamseqInfo.AnnotationName])

    if (countDirectory != None) :

        f.close()

        callR(getPlotter("PCAPlotter") + " -f " + f.name + " -O " + replaceExtension(outputFile, ".pdf", "_PCA") + " -P " + replaceExtension(outputFile, ".txt", "_PCA"), log, dry=printOnly, verbose=verbose)

        print("FileName", "SampleName", "SampleType", "SampleTime", "Sequenced", "Mapped", "Deduplicated", "MQ-Filtered", "Identity-Filtered", "NM-Filtered", "Multimap-Filtered", "Retained", "Counted", "Annotation", sep="\t", file=tsvFile)


    else :
        print("FileName", "SampleName", "SampleType", "SampleTime", "Sequenced", "Mapped", "Deduplicated", "MQ-Filtered", "Identity-Filtered", "NM-Filtered", "Multimap-Filtered", "Retained", "Annotation", sep="\t", file=tsvFile)

    for key in sorted(contentDict):
        print(contentDict[key], file=tsvFile)

    tsvFile.close()

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
            readIterator = testFile.readsInChromosome(chromosome, minQual)

            for read in readIterator:

                tcCounts = [0] * maxReadLength
                mutCounts = [0] * maxReadLength

                for mismatch in read.mismatches:
                    if(mismatch.isTCMismatch(read.direction == ReadDirection.Reverse)):
                        tcCounts[mismatch.readPosition] += 1
                    else :
                        mutCounts[mismatch.readPosition] += 1


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

        print("# slamdunk tcperreadpos v" + __version__, file=foTC)

        for i in range(0, maxReadLength):
            print(allPerPosFwd[i], allPerPosRev[i], tcPerPosFwd[i], tcPerPosRev[i], totalReadCountFwd[i], totalReadCountRev[i], sep='\t', file=foTC)
        foTC.close()

    if(not checkStep([outputCSV], [outputPDF], force)):
        print("Skipped computing T->C per reads position plot for file " + bam, file=log)
    else:
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

            readIterator = testFile.readInRegion(utr.chromosome, utr.start, utr.stop, utr.strand, maxReadLength, minQual)

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

        print("# slamdunk tcperutr v" + __version__, file=foTC)

        reverseAllPerPosRev = allPerPosRev[::-1]
        reverseTcPerPosRev = tcPerPosRev[::-1]
        reverseTotalUtrCountRev = totalUtrCountRev[::-1]

        for i in range(0, utrNormFactor):
            print(allPerPosFwd[i], reverseAllPerPosRev[i], tcPerPosFwd[i], reverseTcPerPosRev[i], totalUtrCountFwd[i], reverseTotalUtrCountRev[i], sep='\t', file=foTC)
        foTC.close()

    if(not checkStep([outputCSV], [outputPDF], force)):
        print("Skipped computing T->C per UTR position plot for file " + bam, file=log)
    else:
        callR(getPlotter("conversion_per_read_position") + " -u -i " + outputCSV + " -o " + outputPDF, log, dry=printOnly, verbose=verbose)

def computeSNPMaskedRates (ref, bed, snpsFile, bam, maxReadLength, minQual, coverageCutoff, variantFraction, outputCSV, outputPDF, strictTCs, log, printOnly=False, verbose=True, force=False):

    if(not checkStep([bam, ref], [outputCSV], force)):
        print("Skipped computing T->C per UTR with SNP masking for file " + bam, file=log)
    else:
        fileCSV = open(outputCSV,'w')

        snps = SNPtools.SNPDictionary(snpsFile)
        snps.read()

        #Go through one chr after the other
        testFile = SlamSeqBamFile(bam, ref, snps)

        progress = 0
        for utr in BedIterator(bed):

            if(not utr.hasStrand()):
                raise RuntimeError("Input BED file does not contain stranded intervals.")

            if utr.start < 0:
                raise RuntimeError("Negativ start coordinate found. Please check the following entry in your BED file: " + utr)

            readIterator = testFile.readInRegion(utr.chromosome, utr.start, utr.stop, utr.strand, maxReadLength, minQual)

            unmaskedTCCount = 0
            maskedTCCount = 0
            readCount = 0

            for read in readIterator:

                # Overwrite any conversions for non-TC reads (reads with < 2 TC conversions)
                if (not read.isTcRead and strictTCs) :
                    read.tcCount = 0
                    read.mismatches = []
                    read.conversionRates = 0.0
                    read.tcRate = 0.0

                isTC = False
                isTrueTC = False

                for mismatch in read.mismatches:
                    if(mismatch.isTCMismatch(read.direction == ReadDirection.Reverse) and mismatch.referencePosition >= 0 and mismatch.referencePosition < utr.getLength()):
                        isTrueTC = True

                    unmasked = False
                    if (read.direction == ReadDirection.Reverse and mismatch.referenceBase == "A" and mismatch.readBase == "G"):
                        unmasked = True
                    elif (read.direction != ReadDirection.Reverse and mismatch.referenceBase == "T" and mismatch.readBase == "C") :
                        unmasked = True

                    if (unmasked and mismatch.referencePosition >= 0 and mismatch.referencePosition < utr.getLength()) :
                        isTC = True

                readCount += 1

                if (isTC) :
                    unmaskedTCCount += 1

                if (isTrueTC) :
                    maskedTCCount += 1

            containsSNP = 0

            if (unmaskedTCCount != maskedTCCount) :
                containsSNP = 1

            print(utr.name + "\t" + str(readCount) + "\t" + str(unmaskedTCCount) + "\t" + str(maskedTCCount) + "\t" + str(containsSNP), file=fileCSV)

            progress += 1

        fileCSV.close()

    if(not checkStep([outputCSV], [outputPDF], force)):
        print("Skipped computing T->C per UTR position plot for file " + bam, file=log)
    else:
        callR(getPlotter("SNPeval") + " -i " + outputCSV + " -c " + str(coverageCutoff) + " -v " + str(variantFraction) + " -o " + outputPDF, log, dry=printOnly, verbose=verbose)

def halflifes(bams, outputCSV, timepoints, log, printOnly=False, verbose=True, force=False):
    callR(getPlotter("compute_halflifes") + " -f " + bams + " -t " + timepoints + " -o " + outputCSV, log, dry=printOnly, verbose=verbose)

def mergeRates(bams, outputCSV, column, columnName, log, printOnly=False, verbose=True, force=False):
    cmd = getPlotter("merge_rate_files") + " -f " + bams + " -o " + outputCSV
    cmd = cmd + " -c \"" + column + "\""
    cmd = cmd + " -n " + str(columnName)
    callR(cmd, log, dry=printOnly, verbose=verbose)
