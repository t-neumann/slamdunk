#!/usr/bin/env python

from __future__ import print_function

import sys
import pysam
import os

from slamdunk.utils.misc import replaceExtension, getSampleInfo, SlamSeqInfo, md5, callR, getPlotter  # @UnresolvedImport
from slamdunk.utils.BedReader import BedIterator  # @UnresolvedImport

from slamdunk.utils import SNPtools  # @UnresolvedImport
from slamdunk.slamseq.SlamSeqFile import SlamSeqBamFile, ReadDirection, SlamSeqInterval  # @UnresolvedImport

from slamdunk.version import __version__, __bam_version__, __count_version__  # @UnresolvedImport

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

def computeTconversions(ref, bed, snpsFile, bam, maxReadLength, minQual, outputCSV, outputBedgraphPlus, outputBedgraphMinus, strictTCs, log, mle = False):
    
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
            refSeq = referenceFile.fetch(region=region).upper()
            
            if (utr.strand == "-") :
                #refSeq = complement(refSeq[::-1])
                Tcontent = refSeq.count("A")
            else :
                Tcontent = refSeq.count("T")
                
            
            slamSeqUtr._Tcontent = Tcontent

        readIterator = testFile.readInRegion(utr.chromosome, utr.start, utr.stop, utr.strand, maxReadLength, minQual)
      
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
            if (not read.isTcRead and strictTCs) :
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
            
            for i in xrange(read.startRefPos, read.endRefPos):
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
            
            for position in xrange(0, len(coverageUtr)):

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

