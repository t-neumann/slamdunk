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

# Date located in: -
from __future__ import print_function
import random
import math
import os
import glob
import sys
import numpy
import pysam
import tempfile

from slamdunk.utils import SNPtools  # @UnresolvedImport
from slamdunk.utils.BedReader import BedIterator, bedToIntervallTree  # @UnresolvedImport
from slamdunk.utils.misc import shell, run, getBinary, getRNASeqReadSimulator, md5, getPlotter, callR  # @UnresolvedImport
from slamdunk.slamseq.SlamSeqFile import SlamSeqBamFile, SlamSeqInterval  # @UnresolvedImport
from slamdunk.version import __version__, __count_version__  # @UnresolvedImport

from Bio import SeqIO
from pybedtools import BedTool

projectPath = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
pathEvalHalfLifes = os.path.join(projectPath, "plot", "eval_halflife_per_gene_plots.R")
pathEvalConversionrates = os.path.join(projectPath, "plot", "eval_conversion_rate_plots.R")
pathEvalHalfLife = os.path.join(projectPath, "plot", "eval_halflifes_error_plot.R")

def getRndBaseWithoutDup(base):
    rndBase = getRndBase()
    while(rndBase == base):
        rndBase = getRndBase()
    return rndBase

def getRndBase():
    return ['A', 'T', 'G', 'C'][random.randrange(0, 3, 1)]

def getCmpBase(base):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}
    return complement[base]

def printUTR(utr, outBed, minLength):
    if utr.getLength() > minLength:
        print(utr.chromosome, utr.start, utr.stop, utr.name.replace("_", "-"), utr.score, utr.strand, sep="\t", file=outBed)

def prepareBED(bed, slamSimBed, minLength):
    utrs = []
    for utr in BedIterator(bed):
        utrs.append(utr)

    utrs.sort(key=lambda x: (x.name, -x.getLength()))

    outBed = open(slamSimBed, "w")

    partList = []
    lastUtr = None
    for utr in utrs:
        if utr.hasStrand() and utr.hasNonEmptyName():
            currentUtr = utr.name
            if currentUtr == lastUtr:
                partList.append(utr)
            else:
                if(not lastUtr is None):
                    printUTR(partList[0], outBed, minLength)
                partList = [utr]
            lastUtr = currentUtr
        else:
            print("Warning: Invalid BED entry found: " + str(utr))

    if(not lastUtr is None):
        printUTR(partList[0], outBed, minLength)

    outBed.close()

def simulateUTR(utrSeq, utr, polyALenght, snpRate, vcfFile):
    utrSeq = list(utrSeq)
    snpCount = 0
    for i in range(0, len(utrSeq)):
        if(random.uniform(0, 1) < snpRate):
            rndBase = getRndBaseWithoutDup(utrSeq[i])
            #print("Introducing SNP " + utrSeq[i] + " -> " + rndBase + " in UTR " + utr.name + " at position " + utr.chromosome + ":" + str(utr.start + i))
            snpPosition = 0
            snpCount += 1
            if(utr.strand == "+"):
                snpPosition = utr.start + i + 1
                snpRef = utrSeq[i]
                snpAlt = rndBase
            else:
                snpPosition = utr.stop - i
                snpRef = getCmpBase(utrSeq[i])
                snpAlt = getCmpBase(rndBase)

            print(utr.chromosome, snpPosition, utr.name + "_" + str(snpCount), snpRef, snpAlt, ".", "PASS", ".", sep="\t", file=vcfFile)
            utrSeq[i] = rndBase
    return "".join(utrSeq) + (polyALenght * 'A')

def prepareUTRs(bed, bed12, bed12Fasta, referenceFasta, readLength, polyALength, explv, snpRate, vcfFile):

    # Read utrs from BED file
    utrs = parseUtrBedFile(bed)

    vcf = open(vcfFile, "w")
    print("##fileformat=VCFv4.1", file=vcf)
    print("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO", file=vcf)

    bedFile = BedTool(bed)

    bedFasta = bedFile.sequence(fi=referenceFasta, s=True, name=True)

    f = tempfile.NamedTemporaryFile(mode='w', delete=False)

    for line in bedFasta.print_sequence().splitlines():
        if(line[0] == ">"):
            print(line.split("::")[0], file=f)
        else:
            print(line.rstrip(), file=f)

    f.close()

    bed12FastaFile = open(bed12Fasta, "w")
    utrName = None
    with open(f.name, 'r') as f:
        for line in f:
            if(line[0] == ">"):
                print(line.rstrip(), file=bed12FastaFile)
                utrName = line.rstrip()[1:]
            else:
                print(simulateUTR(line.rstrip(), utrs[utrName], polyALength, snpRate, vcf).rstrip(), file=bed12FastaFile)
    bed12FastaFile.close()
    vcf.close()

    bed12File = open(bed12, "w")

    totalLength = 0

    minFragmentLength = 150
    maxFragmentLength = 450
    for utr in BedIterator(bed):

        fragmentLength = random.randrange(minFragmentLength, maxFragmentLength, 1) #+ readLength
        fragmentLength = min(fragmentLength, utr.getLength())

        start = max(0, utr.getLength() - fragmentLength)
        end = utr.getLength() #- readLength

        totalLength += (end - start)
#         min(utr.getLength() + readLength / 4, fragmentLength + readLength / 4)
        print(utr.name, start, end, utr.name, utr.score, "+", start, end, "255,0,0", "1", (end - start), 0, sep="\t", file=bed12File)

    bed12File.close()

    output = shell(getRNASeqReadSimulator("genexplvprofile.py") + " --geometric 1 " + bed12 + " 2> /dev/null > " + explv)
    if len(output.strip()) > 5:
        print(output)

    return totalLength

def simulateReads(bed12, bed12Fasta, explv, bedReads, faReads, readLength, readCount, seqError):
    #output = shell(getBinary("gensimreads.py") + " -l " + str(readLength) + " -e " + explv + " -n " + str(readCount) + " -b " + rNASeqReadSimulatorPath + "demo/input/sampleposbias.txt --stranded " + bed12 + " > " + bedReads)

    output = shell(getRNASeqReadSimulator("gensimreads.py") + " -l " + str(readLength) + " -e " + explv + " -n " + str(readCount) + " --stranded " + bed12 + " 2> /dev/null > " + bedReads)
    if len(output.strip()) > 5:
        print(output)

    output = shell(getRNASeqReadSimulator("getseqfrombed.py") + " -f -r " + str(seqError) + " -l " + str(readLength) + " " + bedReads + " " + bed12Fasta + " 2> /dev/null > " + faReads)
    if len(output.strip()) > 5:
        print(output)

def getRndHalfLife(minHalfLife, maxHalfLife):
    return random.randrange(minHalfLife, maxHalfLife, 1)

def simulateTurnOver(bed, turnoverBed, minHalfLife, maxHalfLife):
    turnoverFile = open(turnoverBed, "w")
    for utr in BedIterator(bed):
        print(utr.chromosome, utr.start, utr.stop, utr.name, getRndHalfLife(minHalfLife, maxHalfLife), utr.strand, sep='\t', file=turnoverFile)
    turnoverFile.close()

def printFastaEntry(sequence, name, index, conversions, readOutSAM, conversionRate):
    #a = pysam.AlignedSegment()
    print(name + "_" + str(index) + "_" + str(conversions),
          "4",
          "*",
          "0",
          "0",
          "*",
          "*",
          "0",
          "0",
          sequence,
          "F" * len(sequence),
          "TC:i:" + str(conversions),
          "ID:i:" + str(index),
          "CR:f" + str(conversionRate),
           file=readOutSAM, sep="\t")

def convertRead(read, name, index, conversionRate, readOutSAM):

    tCount = 0
    TcCount = 0
    # TODO: Uncomment to go back to pysam
    #seq = list(read.sequence)
    seq = list(str(read.seq))
    for i in range(0, len(seq)):
        if seq[i] == 'T':
            tCount += 1
            #if conversionRate > 0 and random.uniform(0, 1) <= conversionRate:
            if conversionRate > 0 and numpy.random.binomial(1, conversionRate) == 1:
                seq[i] = 'C'
                TcCount += 1

    printFastaEntry("".join(seq), name, index, TcCount, readOutSAM, conversionRate)

    return tCount, TcCount

def getLambdaFromHalfLife(halfLife):
    return math.log(2) / float(halfLife)

def addTcConversionsToReads(utr, reads, readToConvertPercent, conversionRate, readOutSAM):
    #print(utr.name + " reads found: " + str(len(reads)))


    readsToConvert = int(len(reads) * readToConvertPercent)
    #print("Converting " + str(readsToConvert) + " reads (lambda = " + str(varLambda) + ")")

    totalTCount = 0
    totalTcCount = 0
    readSample = sorted(random.sample(range(0, len(reads)), readsToConvert))

    readsWithTc = 0
    for i in range(0, len(reads)):
        read = reads[i]
        if(i in readSample):
            tCount, TcCount = convertRead(read, utr.name, i, conversionRate, readOutSAM)
            if TcCount > 0:
                readsWithTc += 1
        else:
            tCount, TcCount = convertRead(read, utr.name, i, -1, readOutSAM)
        totalTcCount += TcCount
        totalTCount += tCount

    return readsWithTc, totalTCount, totalTcCount

def getUtrName(readName):
    return readName.split("_")[0]

def printUtrSummary(utr, totalReadCount, readsToConvert, totalTCount, totalTcCount, utrSummary, readsCPM, readToConvertPercent):
    #conversionRate = 0
    #if totalTCount > 0:
    #    conversionRate = totalTcCount * 1.0 / totalTCount
    #print(utr.chromosome, utr.start, utr.stop, utr.name, utr.score, utr.strand, totalReadCount, readsToConvert, totalTCount, totalTcCount, conversionRate, sep="\t", file=utrSummary)
    print(utr.chromosome,
          utr.start,
          utr.stop,
          utr.name, #utr.score,
          utr.stop - utr.start,
          utr.strand,
          readToConvertPercent,
          readsCPM,
          -1,
          totalTCount,
          totalTcCount,
          totalReadCount,
          readsToConvert,
          "-1" ,
          -1.0,
          -1.0,
         sep="\t", file=utrSummary)

def parseUtrBedFile(bed):
    utrs = {}
    for utr in BedIterator(bed):
        utrs[utr.name] = utr
    return utrs

def computeConversionRate(halflife, pulseTimePoint, chaseTimePoint, labeledTranscripts):
    if labeledTranscripts == -1.0:
        # Compute number of labeld transcript from half-life and timepoint
        varLambda = getLambdaFromHalfLife(halflife)
        readToConvertPercent = (1 - math.exp(-varLambda * pulseTimePoint))
        if(chaseTimePoint > 0):
            chaseStart = 1 - readToConvertPercent
            readToConvertPercent = math.exp(-varLambda * chaseTimePoint) - chaseStart
    else:
        # Nuber of labeled transcript provided directly
        readToConvertPercent = labeledTranscripts
    return readToConvertPercent

def addTcConversions(bed, readInFile, readOutFile, pulseTimePoint, chaseTimePoint, utrSummaryFile, conversionRate, librarySize, sampleInfo, labeledTranscripts = -1.0):

    # Read utrs from BED file
    utrs = parseUtrBedFile(bed)

    readOutTemp = readOutFile + "_tmp.sam"
    #bamheader = { 'HD': {'VN': '1.0'} }
    #readOutBAM = pysam.AlignmentFile(readOutTemp, "wb", header=bamheader, add_sq_text=False)
    readOutSAM = open(readOutTemp, "w")
    print("@HD\tVN:1.0\tSO:unsorted", file=readOutSAM)
    utrSummary = open(utrSummaryFile, "w")

    bedMD5 = md5(bed)
    print("#slamdunk v" + __version__, __count_version__, "sample info:", sampleInfo.Name, sampleInfo.ID, sampleInfo.Type, sampleInfo.Time, sep="\t", file=utrSummary)
    print("#annotation:", os.path.basename(bed), bedMD5, sep="\t", file=utrSummary)
    print(SlamSeqInterval.Header, file=utrSummary)

    reads = []
    lastUtrName = None
    utrName = None

    fasta_sequences = SeqIO.parse(open(readInFile),'fasta')

    for entry in fasta_sequences:

    # TODO: Uncomment to go back to pysam
    #with pysam.FastxFile(readInFile) as fh:
        #for entry in fh:
            #utrName = getUtrName(entry.name)
        utrName = getUtrName(entry.id)
        if(utrName == lastUtrName):
            reads.append(entry)
        elif(lastUtrName == None):
            reads.append(entry)
        else:
            readsCPM = len(reads)  * 1000000.0 / librarySize;
            readToConvertPercent = computeConversionRate(utrs[lastUtrName].score, pulseTimePoint, chaseTimePoint, labeledTranscripts)
            readsWithTC, totalTCount, totalTcCount = addTcConversionsToReads(utrs[lastUtrName], reads, readToConvertPercent, conversionRate, readOutSAM)
            printUtrSummary(utrs[lastUtrName], len(reads), readsWithTC, totalTCount, totalTcCount, utrSummary, readsCPM, readToConvertPercent)
            reads = []
        lastUtrName = utrName

    # Last UTR
    readsCPM = len(reads) * 1000000.0 / librarySize;
    readToConvertPercent = computeConversionRate(utrs[lastUtrName].score, pulseTimePoint, chaseTimePoint, labeledTranscripts)
    readsWithTC, totalTCount, totalTcCount = addTcConversionsToReads(utrs[lastUtrName], reads, readToConvertPercent, conversionRate, readOutSAM)
    printUtrSummary(utrs[lastUtrName], len(reads), readsWithTC, totalTCount, totalTcCount, utrSummary, readsCPM, readToConvertPercent)


    readOutSAM.close()
    utrSummary.close()


    readOutTempBAM = readOutFile + "_tmp.bam"
    # Convert to BAM
    run("samtools view -Sb " + readOutTemp + " > " + readOutTempBAM)
    #samFile = pysam.AlignmentFile(readOutTemp, "r", check_header = False, check_sq = False)
    #bamFile = pysam.AlignmentFile(readOutTempBAM, "wb", template=samFile)

    #for read in samFile:
    #    bamFile.write(read)
    #bamFile.close()
    #samFile.close()

    # Sort reads by query name (doesn't matter for mapping, but makes evaluation easier
    #pysam.sort("-o", readOutFile, readOutTempBAM)  # @UndefinedVariable
    run("samtools sort -o " + readOutFile + " " + readOutTempBAM)
    os.unlink(readOutTemp)
    os.unlink(readOutTempBAM)

def getTotalUtrLength(bed12File):
    totalUtrLength = 0
    for utr in BedIterator(bed12File):
        totalUtrLength += utr.getLength()
    return totalUtrLength

def evaluate(simulated, slamdunk, outputFile, log, printOnly=False, verbose=True, force=False):
    cmd = getPlotter("splash_eval_count_files") + " -s " + simulated + " -d " + slamdunk + " -o " + outputFile

    callR(cmd, log, dry=printOnly, verbose=verbose)

def evaluateReads(bam, referenceFile, bed, outputFile, mainOutput):

    print("Run " + bam)

    # Go through one chr after the other
    testFile = SlamSeqBamFile(bam, referenceFile, None)

    chromosomes = testFile.getChromosomes()

    bedTree = bedToIntervallTree(bed)
    #evalHist = [0] *

    outFile = open(outputFile, "w")
    print("read.name", "read.chromosome", "read.startRefPos", "sim.utr", "read.utr", "sim.tcCount", "read.tcCount", sep = "\t", file=outFile)

    total = 0
    correct = 0
    correcPosWrongTC = 0
    wrongPos = 0

    minBaseQual = 0
    for chromosome in chromosomes:
        readIterator = testFile.readsInChromosome(chromosome, minBaseQual)

        for read in readIterator:
            total += 1
            simInfo = read.name.split("_")
            utrSim = simInfo[0]
            tcCountSim = int(simInfo[2])

            utrFound = None
            if read.chromosome in bedTree:
                overlaps = list(bedTree[read.chromosome][read.startRefPos:read.endRefPos])
                if len(overlaps) > 0:
                    utrFound = overlaps[0].data

            if utrFound == utrSim:
                if tcCountSim == read.tcCount:
                    correct += 1
                else:
                    correcPosWrongTC += 1
            else:
                wrongPos += 1

            print(read.name, read.chromosome, read.startRefPos, utrSim, utrFound, tcCountSim, read.tcCount, sep = "\t", file=outFile)

    #print(correct * 100.0 / total, correcPosWrongTC * 100.0 / total, wrongPos * 100.0 / total, total)
    print(correct, correcPosWrongTC, wrongPos, total)

def plotconversiondifferences(simDir, slamDir, conversionRate, outputPDF):

    simFiles = sorted(glob.glob(simDir + "*_utrsummary.csv"))
    slamdunkFiles = sorted(glob.glob(slamDir + "*_reads_slamdunk_mapped_filtered_tcount.csv"))

    if(len(simFiles) == len(slamdunkFiles)):
        run("Rscript " + pathEvalConversionrates + " -c " + str(conversionRate) + " -s " + ",".join(simFiles) + " -f " + ",".join(slamdunkFiles) + " -o " + outputPDF, sys.stderr, dry=False, verbose=False)
    else:
        raise RuntimeError("Couldn't match files with timepoints")


def plotHalfLifes(bed, simDir, slamDir, timePointsStr, conversionRate, outputPDF):

    simFiles = sorted(glob.glob(simDir + "*_utrsummary.csv"))
    slamdunkFiles = sorted(glob.glob(slamDir + "*_reads_slamdunk_mapped_filtered_tcount.csv"))
    timePoints = timePointsStr.split(",")

    if(len(simFiles) == len(slamdunkFiles) and len(slamdunkFiles) == len(timePoints)):
        run("Rscript " + pathEvalHalfLifes + " -b " + bed + " -c " + str(conversionRate) + " -s " + ",".join(simFiles) + " -f " + ",".join(slamdunkFiles) + " -t " + ",".join(timePoints) + " -o " + outputPDF, sys.stderr, dry=False, verbose=False)
    else:
        raise RuntimeError("Couldn't match files with timepoints")

def evalHalfLifes(trueHLFile, simHLFile, predHLFile, outputPDF, erroutputCSV):
    run("Rscript " + pathEvalHalfLife + " -t " + trueHLFile + " -p " + predHLFile + " -s " + simHLFile + " -o " + outputPDF + " -m " + erroutputCSV, sys.stderr, dry=False, verbose=False)

def getConversionRateFromBam(bam, ref, chromosome, start, end, strand):

    testFile = SlamSeqBamFile(bam, ref, SNPtools.SNPDictionary(None))

    sumConversionRate = 0
    readCount = 0
    #for chromosome in testFile.getChromosomes():
    #    readIterator = testFile.readsInChromosome(chromosome)
    #readIterator = testFile.readInRegion("chr7", 3217778, 3221036, "+", 55)
    readIterator = testFile.readInRegion(chromosome, start, end, strand, 100)

    for read in readIterator:
        conversionRate = 0
        if(read.tCount > 0):
            conversionRate = read.tcCount * 1.0 / read.tCount

        #if(read.tcCount > 0):
        sumConversionRate += conversionRate
        readCount += 1

        #if(readCount % 1000 == 0 and readCount > 0):
        #    print(str(readCount) + ": " + str(sumConversionRate) + " / " + str(readCount) + " = " + str(sumConversionRate / readCount))
        #if(readCount >= 10000):
        #    break

    print("Read count: " + str(readCount))
    print("Avg. conversion rate: " + str(sumConversionRate / readCount))
