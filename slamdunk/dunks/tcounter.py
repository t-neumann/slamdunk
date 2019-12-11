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

import sys
import pysam
import os
import re

from os.path import basename

from slamdunk.utils.misc import replaceExtension, getSampleInfo, SlamSeqInfo, md5, callR, getPlotter  # @UnresolvedImport
from slamdunk.utils.BedReader import BedIterator  # @UnresolvedImport

from slamdunk.utils import SNPtools  # @UnresolvedImport
from slamdunk.slamseq.SlamSeqFile import SlamSeqBamFile, ReadDirection, SlamSeqInterval  # @UnresolvedImport

from slamdunk.version import __version__, __bam_version__, __count_version__  # @UnresolvedImport

def pysamIndex(outputBam):
    pysam.index(outputBam)  # @UndefinedVariable

def collapse(expandedCSV, collapsedCSV, log):

    tcDict = {}

    outCSV = open(collapsedCSV, 'w')

    readNumber = 0

    with open(expandedCSV, 'r') as f:

        # Skip header
#         next(f)

        for line in f:

            if (not line.startswith('#') and not line.startswith("Chromosome")) :
                fields = line.rstrip().split('\t')

                # For now, ignore everything after column 14
                if (len(fields) >= 14) :

                    gene = fields[3]
                    length = fields[4]
                    Tcontent = fields[8]
                    coverageOnTs = fields[9]
                    conversionsOnTs = fields[10]
                    readCount = fields[11]
                    tcReadCount = fields[12]
                    multimapCount = fields[13]

                    if (gene in tcDict.keys()) :
                        tcDict[gene]['length'] += max(int(length),0)
                        tcDict[gene]['Tcontent'] += max(int(Tcontent),0)
                        tcDict[gene]['coverageOnTs'] += max(int(coverageOnTs),0)
                        tcDict[gene]['conversionsOnTs'] += max(int(conversionsOnTs),0)
                        tcDict[gene]['readCount'] += max(int(readCount),0)
                        tcDict[gene]['tcReadCount'] += max(int(tcReadCount),0)
                        tcDict[gene]['multimapCount'] += max(int(multimapCount),0)
                    else :
                        tcDict[gene] = {}
                        tcDict[gene]['length'] = max(int(length),0)
                        tcDict[gene]['Tcontent'] = max(int(Tcontent),0)
                        tcDict[gene]['coverageOnTs'] = max(int(coverageOnTs),0)
                        tcDict[gene]['conversionsOnTs'] = max(int(conversionsOnTs),0)
                        tcDict[gene]['readCount'] = max(int(readCount),0)
                        tcDict[gene]['tcReadCount'] = max(int(tcReadCount),0)
                        tcDict[gene]['multimapCount'] = max(int(multimapCount),0)

                    readNumber += int(readCount)

                else :
                    print("Error in TC file format - unexpected number of fields (" + str(len(fields)) + ") in the following line:\n" + line, file=log)

    print("gene_name", "length", "readsCPM", "conversionRate", "Tcontent", "coverageOnTs", "conversionsOnTs", "readCount", "tcReadCount", "multimapCount", sep='\t', file=outCSV)

    for gene in sorted(tcDict.keys()) :

        print(gene,end="\t",file=outCSV)
        print(tcDict[gene]['length'],end="\t",file=outCSV)
        print(float(tcDict[gene]['readCount']) / float(readNumber) * 1000000,end="\t",file=outCSV)
        conversionRate = 0
        if (tcDict[gene]['coverageOnTs'] > 0) :
            conversionRate = float(tcDict[gene]['conversionsOnTs']) / tcDict[gene]['coverageOnTs']
        print(conversionRate,end="\t",file=outCSV)
        print(tcDict[gene]['Tcontent'],end="\t",file=outCSV)
        print(tcDict[gene]['coverageOnTs'],end="\t",file=outCSV)
        print(tcDict[gene]['conversionsOnTs'],end="\t",file=outCSV)
        print(tcDict[gene]['readCount'],end="\t",file=outCSV)
        print(tcDict[gene]['tcReadCount'],end="\t",file=outCSV)
        print(tcDict[gene]['multimapCount'],file=outCSV)

    outCSV.close()

def getMean(values, skipZeros=True):
    count = 0.0
    totalSum = 0.0
    for value in values:
        if not skipZeros or value > 0:
            totalSum = totalSum + value
            count += 1
    if(count > 0):
        return totalSum / count
    else:
        return 0.0

def computeTconversions(ref, bed, snpsFile, bam, maxReadLength, minQual, outputCSV, outputBedgraphPlus, outputBedgraphMinus, conversionThreshold, log, mle = False):

    referenceFile = pysam.FastaFile(ref)

    sampleInfo = getSampleInfo(bam)

    slamseqInfo = SlamSeqInfo(bam)
    #readNumber = slamseqInfo.MappedReads
    readNumber = slamseqInfo.FilteredReads

    bedMD5 = md5(bed)

    if(mle):
        fileNameTest = replaceExtension(outputCSV, ".tsv", "_perread")
        fileTest = open(fileNameTest,'w')
        print("#slamdunk v" + __version__, __count_version__, "sample info:", sampleInfo.Name, sampleInfo.ID, sampleInfo.Type, sampleInfo.Time, sep="\t", file=fileTest)
        print("#annotation:", os.path.basename(bed), bedMD5, sep="\t", file=fileTest)
        #print("utr", "n", "k", file=fileTest)
        print(SlamSeqInterval.Header, file=fileTest)


    fileCSV = open(outputCSV,'w')
    print("#slamdunk v" + __version__, __count_version__, "sample info:", sampleInfo.Name, sampleInfo.ID, sampleInfo.Type, sampleInfo.Time, sep="\t", file=fileCSV)
    print("#annotation:", os.path.basename(bed), bedMD5, sep="\t", file=fileCSV)
    print(SlamSeqInterval.Header, file=fileCSV)

    snps = SNPtools.SNPDictionary(snpsFile)
    snps.read()

    #Go through one chr after the other
    testFile = SlamSeqBamFile(bam, ref, snps)
    if not testFile.bamVersion == __bam_version__:
        raise RuntimeError("Wrong filtered BAM file version detected (" + testFile.bamVersion + "). Expected version " + __bam_version__ + ". Please rerun slamdunk filter.")

    bedMD5 = md5(bed)
    if slamseqInfo.AnnotationMD5 != bedMD5:
        print("Warning: MD5 checksum of annotation (" + bedMD5 + ") does not matched MD5 in filtered BAM files (" + slamseqInfo.AnnotationMD5 + "). Most probably the annotation filed changed after the filtered BAM files were created.", file=log)

    conversionBedGraph = {}

    for utr in BedIterator(bed):
        Tcontent = 0
        slamSeqUtr = SlamSeqInterval(utr.chromosome, utr.start, utr.stop, utr.strand, utr.name, Tcontent, 0, 0, 0, 0, 0, 0, 0)
        slamSeqUtrMLE = SlamSeqInterval(utr.chromosome, utr.start, utr.stop, utr.strand, utr.name, Tcontent, 0, 0, 0, 0, 0, 0, 0)
        if(not utr.hasStrand()):
            raise RuntimeError("Input BED file does not contain stranded intervals.")

        if utr.start < 0:
            raise RuntimeError("Negativ start coordinate found. Please check the following entry in your BED file: " + utr)
        # Retreive reference sequence
        region = utr.chromosome + ":" + str(utr.start + 1) + "-" + str(utr.stop)

        if(utr.chromosome in list(referenceFile.references)):
            #print(refRegion,file=sys.stderr)
            # pysam-0.15.0.1
            #refSeq = referenceFile.fetch(region=region).upper()
            refSeq = referenceFile.fetch(reference=utr.chromosome, start=utr.start, end=utr.stop).upper()
            if (utr.strand == "-") :
                #refSeq = complement(refSeq[::-1])
                Tcontent = refSeq.count("A")
            else :
                Tcontent = refSeq.count("T")


            slamSeqUtr._Tcontent = Tcontent

        readIterator = testFile.readInRegion(utr.chromosome, utr.start, utr.stop, utr.strand, maxReadLength, minQual, conversionThreshold)

        tcCountUtr = [0] * utr.getLength()
        coverageUtr = [0] * utr.getLength()

        tInReads = []
        tcInRead = []

        countFwd = 0
        tcCountFwd = 0
        countRev = 0
        tCountRev = 0

        multiMapFwd = 0
        multiMapRev = 0

        for read in readIterator:

            # Overwrite any conversions for non-TC reads (reads with < 2 TC conversions)
            if (not read.isTcRead) :
                read.tcCount = 0
                read.mismatches = []
                read.conversionRates = 0.0
                read.tcRate = 0.0

            if(read.direction == ReadDirection.Reverse):
                countRev += 1
                if read.tcCount > 0:
                    tCountRev += 1
                if read.isMultimapper:
                    multiMapRev += 1
            else:
                countFwd += 1
                if read.tcCount > 0:
                    tcCountFwd += 1
                if read.isMultimapper:
                    multiMapFwd += 1

            for mismatch in read.mismatches:
                if(mismatch.isTCMismatch(read.direction == ReadDirection.Reverse) and mismatch.referencePosition >= 0 and mismatch.referencePosition < utr.getLength()):
                    tcCountUtr[mismatch.referencePosition] += 1

            testN = read.getTcount()
            testk = 0
            for mismatch in read.mismatches:
                if(mismatch.referencePosition >= 0 and mismatch.referencePosition < utr.getLength()):
                    if(mismatch.isT(read.direction == ReadDirection.Reverse)):
                        testN += 1
                    if(mismatch.isTCMismatch(read.direction == ReadDirection.Reverse)):
                        testk += 1
            #print(utr.name, read.name, read.direction, testN, testk, read.sequence, sep="\t")
            tInReads.append(testN)
            tcInRead.append(testk)
            #print(utr.name, testN, testk, sep="\t", file=fileTest)

            for i in range(read.startRefPos, read.endRefPos):
                if(i >= 0 and i < utr.getLength()):
                    coverageUtr[i] += 1


        if((utr.strand == "+" and countFwd > 0) or (utr.strand == "-" and countRev > 0)):
            tcRateUtr = [ x * 100.0 / y if y > 0 else 0 for x, y in zip(tcCountUtr, coverageUtr)]

            readCount = countFwd
            tcReadCount = tcCountFwd
            multiMapCount = multiMapFwd

            if(utr.strand == "-"):
                readCount = countRev
                tcReadCount = tCountRev
                multiMapCount = multiMapRev

            if((utr.strand == "-" and countFwd > countRev) or (utr.strand == "+" and countRev > countFwd)):
                print("Warning: " + utr.name + " is located on the " + utr.strand + " strand but read counts are higher for the opposite strand (fwd: " + countFwd + ", rev: " + countRev + ")", file=sys.stderr)


            refSeq = readIterator.getRefSeq()

            # Get number of covered Ts/As in the UTR and compute average conversion rate for all covered Ts/As
            coveredTcount = 0
            avgConversationRate = 0
            coveredPositions = 0
            # Get number of reads on T positions and number of reads with T->C conversions on T positions
            coverageOnTs = 0
            conversionsOnTs = 0

            for position in range(0, len(coverageUtr)):

                if(coverageUtr[position] > 0 and ((utr.strand == "+" and refSeq[position] == "T") or (utr.strand == "-" and refSeq[position] == "A"))):
                    coveredTcount += 1
                    avgConversationRate += tcRateUtr[position]

                    coverageOnTs += coverageUtr[position]
                    conversionsOnTs += tcCountUtr[position]
                    conversionBedGraph[utr.chromosome + ":" + str(utr.start + position) + ":" + str(utr.strand)] = tcRateUtr[position]
                if(coverageUtr[position] > 0):
                    coveredPositions += 1

            if(coveredTcount > 0):
                avgConversationRate = avgConversationRate / coveredTcount
            else:
                avgConversationRate = 0

            # reads per million mapped to the UTR
            readsCPM = 0
            if(readNumber > 0):
                readsCPM = readCount  * 1000000.0 / readNumber


            # Convert to SlamSeqInterval and print
            conversionRate = 0
            if (coverageOnTs > 0) :
                conversionRate = float(conversionsOnTs) / float(coverageOnTs)
            slamSeqUtr = SlamSeqInterval(utr.chromosome, utr.start, utr.stop, utr.strand, utr.name, Tcontent, readsCPM, coverageOnTs, conversionsOnTs, conversionRate, readCount, tcReadCount, multiMapCount)
            slamSeqUtrMLE = SlamSeqInterval(utr.chromosome, utr.start, utr.stop, utr.strand, utr.name, Tcontent, readsCPM, coverageOnTs, conversionsOnTs, conversionRate, ",".join(str(x) for x in tInReads), ",".join(str(x) for x in tcInRead), multiMapCount)

        print(slamSeqUtr, file=fileCSV)
        if(mle):
            print(slamSeqUtrMLE, file=fileTest)

    fileCSV.close()
    if(mle):
        fileTest.close()

    fileBedgraphPlus = open(outputBedgraphPlus,'w')
    fileBedgraphMinus = open(outputBedgraphMinus,'w')

    for position in conversionBedGraph:
        positionData = position.split(":")
        if(positionData[2] == "+"):
            print(positionData[0], positionData[1], int(positionData[1]) + 1, conversionBedGraph[position], file=fileBedgraphPlus)
        else:
            print(positionData[0], positionData[1], int(positionData[1]) + 1, conversionBedGraph[position], file=fileBedgraphMinus)

    fileBedgraphPlus.close()
    fileBedgraphMinus.close()

    if(mle):
        fileNameMLE = replaceExtension(outputCSV, ".tsv", "_mle")
        callR(getPlotter("compute_conversion_rate_mle") +  " -f " + fileNameTest + " -r " + "0.024" + " -o " + fileNameMLE + " &> /dev/null")

def genomewideConversionRates(referenceFile, snpsFile, bam, minBaseQual, outputBedGraphPrefix, conversionThreshold, coverageCutoff, log):

    ref = pysam.FastaFile(referenceFile)

    snps = SNPtools.SNPDictionary(snpsFile)
    snps.read()

    # Go through one chr after the other
    testFile = SlamSeqBamFile(bam, referenceFile, snps)

    chromosomes = testFile.getChromosomes()

    bedGraphInfo = re.sub("_slamdunk_mapped.*","",basename(outputBedGraphPrefix))
    print(bedGraphInfo)

    fileBedGraphRatesPlus = open(outputBedGraphPrefix + "_TC_rates_genomewide.bedGraph", 'w')
    fileBedGraphRatesMinus = open(outputBedGraphPrefix + "_AG_rates_genomewide.bedGraph", 'w')
    fileBedGraphCoveragePlus = open(outputBedGraphPrefix + "_coverage_plus_genomewide.bedGraph", 'w')
    fileBedGraphCoverageMinus = open(outputBedGraphPrefix + "_coverage_minus_genomewide.bedGraph", 'w')
    fileBedGraphTCConversions = open(outputBedGraphPrefix + "_TC_conversions_genomewide.bedGraph", 'w')
    fileBedGraphAGConversions = open(outputBedGraphPrefix + "_AG_conversions_genomewide.bedGraph", 'w')
    fileBedGraphT = open(outputBedGraphPrefix + "_coverage_T_genomewide.bedGraph", 'w')
    fileBedGraphA = open(outputBedGraphPrefix + "_coverage_A_genomewide.bedGraph", 'w')

    print("track type=bedGraph name=\"" + bedGraphInfo + " tc-conversions\" description=\"# T->C conversions / # reads on T per position genome-wide\"", file=fileBedGraphRatesPlus)
    print("track type=bedGraph name=\"" + bedGraphInfo + " ag-conversions\" description=\"# A->G conversions / # reads on A per position genome-wide\"", file=fileBedGraphRatesMinus)
    print("track type=bedGraph name=\"" + bedGraphInfo + " plus-strand coverage\" description=\"# Reads on plus strand genome-wide\"", file=fileBedGraphCoveragePlus)
    print("track type=bedGraph name=\"" + bedGraphInfo + " minus-strand coverage\" description=\"# Reads on minus strand genome-wide\"", file=fileBedGraphCoverageMinus)
    print("track type=bedGraph name=\"" + bedGraphInfo + " T->C conversions\" description=\"# T->C conversions on plus strand genome-wide\"", file=fileBedGraphTCConversions)
    print("track type=bedGraph name=\"" + bedGraphInfo + " A->G conversions\" description=\"# A->G conversions on minus strand genome-wide\"", file=fileBedGraphAGConversions)
    print("track type=bedGraph name=\"" + bedGraphInfo + " T-coverage\" description=\"# Plus-strand reads on Ts genome-wide\"", file=fileBedGraphT)
    print("track type=bedGraph name=\"" + bedGraphInfo + " A-coverage\" description=\"# Minus-strand reads on As genome-wide\"", file=fileBedGraphA)

    for chromosome in chromosomes:

        chrLength = testFile.getChromosomeLength(chromosome)

        tcCount = [0] * chrLength
        agCount = [0] * chrLength

        coveragePlus = [0] * chrLength
        coverageMinus = [0] * chrLength

        tCoverage = [0] * chrLength
        aCoverage = [0] * chrLength

        readIterator = testFile.readsInChromosome(chromosome, minBaseQual, conversionThreshold)

        for read in readIterator:
            if (not read.isTcRead) :
                read.tcCount = 0
                read.mismatches = []
                read.conversionRates = 0.0
                read.tcRate = 0.0

            for mismatch in read.mismatches:
                if(mismatch.isTCMismatch(read.direction == ReadDirection.Reverse) and mismatch.referencePosition >= 0 and mismatch.referencePosition < chrLength):
                    if read.direction == ReadDirection.Reverse:
                        agCount[mismatch.referencePosition] += 1
                    else :
                        tcCount[mismatch.referencePosition] += 1

            for i in range(read.startRefPos, read.endRefPos):
                if(i >= 0 and i < chrLength):
                    if read.direction == ReadDirection.Reverse:
                        coverageMinus[i] += 1
                    else :
                        coveragePlus[i] += 1

        prevCoveragePlus = 0
        prevCoveragePlusPos = 0
        prevCoverageMinus = 0
        prevCoverageMinusPos = 0
        prevTCConversionRate = 0
        prevTCConversionRatePos = 0
        prevAGConversionRate = 0
        prevAGConversionRatePos = 0
        prevTCConversions = 0
        prevTCConversionPos = 0
        prevAGConversions = 0
        prevAGConversionPos = 0
        prevTCoverage = 0
        prevTCoveragePos = 0
        prevACoverage = 0
        prevACoveragePos = 0

        for pos in range(0, chrLength):
            if prevCoveragePlus != coveragePlus[pos]:
                print(chromosome + "\t" + str(prevCoveragePlusPos + 1) + "\t" + str(pos + 1) + "\t" + str(prevCoveragePlus), file = fileBedGraphCoveragePlus)
                prevCoveragePlus = coveragePlus[pos]
                prevCoveragePlusPos = pos
            if prevCoverageMinus != coverageMinus[pos]:
                print(chromosome + "\t" + str(prevCoverageMinusPos + 1) + "\t" + str(pos + 1) + "\t" + str(prevCoverageMinus), file = fileBedGraphCoverageMinus)
                prevCoverageMinus = coverageMinus[pos]
                prevCoverageMinusPos = pos

            tCoverage = 0

            if coveragePlus[pos] > 0:
                base = ref.fetch(reference=chromosome, start = pos + 1, end = pos + 2)
                if base.upper() == "T":
                    tCoverage = coveragePlus[pos]

            aCoverage = 0

            if coverageMinus[pos] > 0:
                base = ref.fetch(reference=chromosome, start = pos + 1, end = pos + 2)
                if base.upper() == "A":
                    aCoverage = coverageMinus[pos]

            if prevTCoverage != tCoverage:
                print(chromosome + "\t" + str(prevTCoveragePos + 1) + "\t" + str(pos + 1) + "\t" + str(prevTCoverage), file = fileBedGraphT)
                prevTCoverage = tCoverage
                prevTCoveragePos = pos

            if prevACoverage != aCoverage:
                print(chromosome + "\t" + str(prevACoveragePos + 1) + "\t" + str(pos + 1) + "\t" + str(prevACoverage), file = fileBedGraphA)
                prevACoverage = aCoverage
                prevACoveragePos = pos

            if prevTCConversions != tcCount[pos]:
                print(chromosome + "\t" + str(prevTCConversionPos + 1) + "\t" + str(pos + 1) + "\t" + str(prevTCConversions), file = fileBedGraphTCConversions)
                prevTCConversions = tcCount[pos]
                prevTCConversionPos = pos

            if prevAGConversions != agCount[pos]:
                print(chromosome + "\t" + str(prevAGConversionPos + 1) + "\t" + str(pos + 1) + "\t" + str(prevAGConversions), file = fileBedGraphAGConversions)
                prevAGConversions = agCount[pos]
                prevAGConversionPos = pos

            TCconversionRate = 0
            if coveragePlus[pos] > 0 and coveragePlus[pos] >= coverageCutoff:
                TCconversionRate = float(tcCount[pos]) / float(coveragePlus[pos])

            AGconversionRate = 0
            if coverageMinus[pos] > 0 and coverageMinus[pos] >= coverageCutoff:
                AGconversionRate = float(agCount[pos]) / float(coverageMinus[pos])

            if prevTCConversionRate != TCconversionRate:
                print(chromosome + "\t" + str(prevTCConversionRatePos + 1) + "\t" + str(pos + 1) + "\t" + str(prevTCConversionRate), file = fileBedGraphRatesPlus)
                prevTCConversionRate = TCconversionRate
                prevTCConversionRatePos = pos

            if prevAGConversionRate != AGconversionRate:
                print(chromosome + "\t" + str(prevAGConversionRatePos + 1) + "\t" + str(pos + 1) + "\t" + str(prevAGConversionRate), file = fileBedGraphRatesMinus)
                prevAGConversionRate = AGconversionRate
                prevAGConversionRatePos = pos

    fileBedGraphRatesPlus.close()
    fileBedGraphRatesMinus.close()
    fileBedGraphCoveragePlus.close()
    fileBedGraphCoverageMinus.close()
    fileBedGraphTCConversions.close()
    fileBedGraphAGConversions.close()
    fileBedGraphT.close()
    fileBedGraphA.close()

def genomewideReadSeparation(referenceFile, snpsFile, bam, minBaseQual, outputBAMPrefix, conversionThreshold, log):

    ref = pysam.FastaFile(referenceFile)

    snps = SNPtools.SNPDictionary(snpsFile)
    snps.read()

    # Go through one chr after the other
    testFile = SlamSeqBamFile(bam, referenceFile, snps)

    samFile = pysam.AlignmentFile(bam, "rb")

    chromosomes = testFile.getChromosomes()

    backgroundReadFileName = outputBAMPrefix + "_backgroundReads.bam"
    tcReadFileName = outputBAMPrefix + "_TCReads.bam"

    backgroundReadFile = pysam.AlignmentFile(backgroundReadFileName, "wb", template=samFile)
    tcReadFile = pysam.AlignmentFile(tcReadFileName, "wb", template=samFile)

    tcReadDict = dict()

    for chromosome in chromosomes:

        readIterator = testFile.readsInChromosome(chromosome, minBaseQual, conversionThreshold)

        for read in readIterator:
            if (read.isTcRead) :
                tcReadDict[read.name] = 0

    for read in samFile.fetch():
        if read.query_name in tcReadDict:
            tcReadFile.write(read)
        else:
            backgroundReadFile.write(read)

    backgroundReadFile.close()
    tcReadFile.close()

    pysamIndex(backgroundReadFileName)
    pysamIndex(tcReadFileName)
