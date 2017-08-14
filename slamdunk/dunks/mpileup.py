#!/usr/bin/env python

from __future__ import print_function
import sys, os

from argparse import ArgumentParser, RawDescriptionHelpFormatter
    
import re, os, sys
import pysam, random
from collections import defaultdict
from __builtin__ import False

from fisher import pvalue

def getFisher(expReads1, expReads2, obsReads1, obsReads2):
        
    if expReads1 < 0 : expReads1 = 0
    if expReads2 < 0 : expReads2 = 0
    if obsReads1 < 0 : obsReads1 = 0
    if obsReads2 < 0 : obsReads2 = 0
    
    p = pvalue(obsReads1, obsReads2, expReads1, expReads2)
    
    if (p.right_tail > p.left_tail) :
        return p.left_tail
    else :
        return p.right_tail
    

def getSignificance(obsReads1, obsReads2) :
    
    pValue = 1.0
    baselineError = 0.001
    
    coverage = obsReads1 + obsReads2
    
    expReads2 = int(coverage * baselineError)
    expReads1 = coverage - expReads2
    
    pValue = getFisher(obsReads1, obsReads2, expReads1, expReads2)
    
    return pValue

def getShortIndel(gt):
    
    if ("INS" in gt or "DEL" in gt) :
        try:
            gtContents = gt.split("-")
            indelType = gtContents[0]
            indelSize = gtContents[1]
            indelBases = gtContents[2]
            
            if ("INS" in indelType) :
                return "+" + indelBases
            else :
                return "-" + indelBases
        except Exception as ex:
            print("Warning: error generating consensu from " + gt, file = sys.stderr)
            
    return "N"

def genotypeToCode(gt) :
    
    if (gt == "AA"): return "A"
    if (gt == "CC"): return "C"
    if (gt == "GG"): return "G"
    if (gt == "TT"): return "T"
    if (gt == "AC" or gt == "CA"): return "M"
    if (gt == "AG" or gt == "GA"): return "R"
    if (gt == "AT" or gt == "TA"): return "W"
    if (gt == "CG" or gt == "GC"): return "S"
    if (gt == "CT" or gt == "TC"): return "Y"
    if (gt == "GT" or gt == "TG"): return "K"
    if (gt[0] == "N") : return "N"
    
    return gt
     
def callPosition(refBase, readCounts, callType, minReads2, minVarFreq, minAvgQual, pValueThreshold, minFreqForHom) :
    
    callResult = ""
    
    reads1 = 0
    reads2 = 0
    readsWithIndels = 0
    strands1 = 0
    strands2 = 0
    avgQual1 = 0
    avgQual2 = 0
    avgMap1 = 0
    avgMap2 = 0
    reads1indel = 0
    reads1plus = 0
    reads1minus = 0
    reads2plus = 0
    reads2minus = 0
    pValue = 1
    varFreq = 0.00
    varAllele = ""
    
    try :
        
        try: 
            if (refBase in readCounts) :
                refBaseContent = readCounts[refBase].split("\t")
                reads1 = int(refBaseContent[0])
                strands1 = int(refBaseContent[1])
                avgQual1 = int(refBaseContent[2])
                avgMap1 = int(refBaseContent[3])
                reads1plus = int(refBaseContent[4])
                reads1minus = int(refBaseContent[5])
                
                if (len(refBaseContent) > 6) :
                    reads1indel = reads1minus = int(refBaseContent[6])
        except Exception as ex:
            print("refBase readcount error: " + readCounts[refBase], file = sys.stderr)
            
        alleleKeys = readCounts.keys()
        alleleKeys.sort()
        
        totalReadCounts = 0
        
        for allele in alleleKeys:
            alleleContents = readCounts[allele].split("\t")
            try :
                thisReads = int(alleleContents[0])
                totalReadCounts += thisReads
            except Exception as ex:
                pass
            
        for allele in alleleKeys:
            alleleContents = readCounts[allele].split("\t")
            
            if (allele != refBase) :
                thisReads1      = reads1
                thisReads2      = 0
                thisStrands2    = 0
                thisAvgQual2    = 0
                thisAvgMap2     = 0
                thisReads2plus  = 0
                thisReads2minus = 0
                
                try:
                    thisReads2 = int(alleleContents[0])
                    thisStrands2 = int(alleleContents[1])
                    thisAvgQual2 = int(alleleContents[2])
                    thisAvgMap2 = int(alleleContents[3])
                    thisReads2plus = int(alleleContents[4])
                    thisReads2minus = int(alleleContents[5])
                    
                    if ("INS" in allele or "DEL" in allele) :
                        readsWithIndels += thisReads2
                except Exception as ex:
                    pass
                
        if (callType != "CNS" or thisReads2 > reads2) :
            thisVarFreq = float(thisReads2) / float(totalReadCounts)
            thisPvalue = 1
            
            if ("INS" in allele or "DEL" in allele) :
                thisTotalReadCounts = totalReadCounts - reads1indel
                if (thisTotalReadCounts < thisReads2) :
                    thisTotalReadCounts = thisReads2
                    
                thisVarFreq = float(thisReads2) / float(totalReadCounts)
                
            if (pValueThreshold == 0.99) :
                thisPvalue = 0.98
            else :
                thisPvalue = getSignificance(reads1, thisReads2)
            
            if (thisReads2 > reads2 and thisAvgMap2 >= minAvgQual) :
                if ("INS" in allele or "DEL" in allele) :
                    varAllele = getShortIndel(allele)
                else :
                    varAllele = allele
                    
                reads2 = thisReads2
                strands2 = thisStrands2
                avgQual2 = thisAvgQual2
                avgMap2 = thisAvgMap2
                reads2plus = thisReads2plus
                reads2minus = thisReads2minus
                varFreq = thisVarFreq * 100.0
                pValue = thisPvalue
                
            if (thisReads2 >= minReads2 and thisAvgQual2 >= minAvgQual and thisVarFreq >= minVarFreq) :
                thisReads1 = reads1
                thisVarFreq = thisVarFreq * 100.0
                
                thisVarType = "SNP"
                
                if ("INS" in allele or "DEL" in allele) :
                    thisVarType = "INDEL"
                    thisReads1 = reads1
                    if (thisReads1 < 0) :
                        thisReads1 = 0
                    allele = getShortIndel(allele)
                    
                if (thisPvalue <= pValueThreshold) :
                    if (callType == "SNP" or callType == "INDEL") :
                        
                        reads2 = thisReads2
                        strands2 = thisStrands2
                        avgQual2 = thisAvgQual2
                        avgMap2 = thisAvgMap2
                        reads2plus = thisReads2plus
                        reads2minus = thisReads2minus
                        pValue = thisPvalue
                        
                        genotype = ""
                        
                        if (thisVarFreq >= (float(minFreqForHom) * 100)) :
                            genotype = allele + allele
                            if thisVarType == "INDEL" :
                                genotype = allele + "/" + allele
                        else:
                            genotype = refBase + allele
                            if thisVarType == "INDEL":
                                genotype = "*/" + allele
                                
                        if thisVarType == callType:
                            if (len(callResult) > 0) :
                                callResult += "\n"
                            
                            if (thisReads1 < 0) :
                                thisReads1 = 0
                                
                            if (reads2 < 0) :
                                reads2 = 0
                                
                            callResult += genotypeToCode(genotype) + "\t" + str(thisReads1) + "\t" + str(reads2) + "\t" + str(thisVarFreq) + "%\t" + str(strands1) + "\t" 
                            callResult += str(strands2) + "\t" + str(avgQual1) + "\t" + str(avgQual2) + "\t" + str(pValue) + "\t" + str(avgMap1) + "\t"
                            callResult += str(avgMap2) + "\t" + str(reads1plus) + "\t" + str(reads1minus) + "\t" + str(reads2plus) + "\t" 
                            callResult += str(reads2minus) + "\t" + str(varAllele)
                                
    except Exception as ex:
        pass
    
    if (len(callResult) == 0 and callType == "CNS") :
        if (reads1 > 0 and reads1 > minReads2) :
            callResult = str(refBase) + "\t" + str(reads1) + "\t" + str(reads2) + "\t" + str(varFreq) + "%\t" + str(strands1) + "\t" + str(strands2) + "\t" 
            callResult += str(avgQual1) + "\t" + str(avgQual2) + "\t" + str(pValue) + "\t" + str(avgMap1) + "\t" + str(avgMap2)
            callResult += "\t" + str(reads1plus) + "\t" + str(reads1minus) + "\t" + str(reads2plus) + "\t" + str(reads2minus) + "\t" + str(varAllele)
        else:
            callResult = "N" + "\t" + str(reads1) + "\t" + str(reads2) + "\t" + str(varFreq) + "%\t" + str(strands1) + "\t" + str(strands2)
            callResult += "\t" + str(avgQual1) + "\t" + str(avgQual2) + "\t" + str(pValue) + "\t" + str(avgMap1) + "\t" + str(avgMap2)             
            callResult += "\t" + str(reads1plus) + "\t" + str(reads1minus) + "\t" + str(reads2plus) + "\t" + str(reads2minus) + "\t" + str(varAllele)

    return callResult

def getReadCounts(refBase, readBases, readQuals, minAvgQual, mapQuals):
    readCounts = dict()
    readCountsPlus = dict()
    readCountsMinus = dict()
    qualitySum = dict()
    mapQualitySum = dict()
    strandsSeen = dict()
    
    reads1 = 0
    reads1indel = 0
    readBase = ""
    prevBase = ""
    nextBase = ""
    baseQ = 0
    prevBaseQ = 0
    mapQ = 1
    strand = ""
    
    readStart = False
    j = 0
    
    for i, c in enumerate(readBases):
        readBase = c
        if(i == 0 and len(readBase) == 0):
            
            i += 1
            readBase = readBases[i]
                        
        prevBase = ""
        if(i > 1 and i < len(readBases) - 1) :
            prevBase = readBases[i - 1]
        if (j > 1 and j < len(readQuals) - 1) :
            prevBaseQ = ord(readQuals[j - 1]) - 33
            
        nexBase = ""
        if (i < len(readBases) - 1) :
            nextBase = readBases[i + 1]
            
        if (j < len(readQuals)):
            baseQ = ord(readQuals[j]) - 33
            
        if (j < len(mapQuals)) :
            mapQ = ord(mapQuals[j]) - 33
            
        if (readBase == "." or readBase == "," and nextBase != "-" or nextBase == "+") :
            strand = "+"
            if (readBase == ",") :
                strand = "-"
            
            if (baseQ >= minAvgQual) :
                reads1 += 1
                
                if ("ref" in strandsSeen) :
                    alreadySeen = strandsSeen["ref"]
                    if (not(len(alreadySeen) >= 2 or alreadySeen == strand)):
                        strandsSeen["ref"] = strandsSeen["ref"] + strand
                else :
                    strandsSeen["ref"] = strand
                    
                if (strand == "+") :
                    if ("ref" in readCountsPlus) :
                        readCountsPlus["ref"] += 1
                    else :
                         readCountsPlus["ref"] = 1
                else :
                    if ("ref" in readCountsMinus) :
                        readCountsMinus["ref"] += 1
                    else :
                         readCountsMinus["ref"] = 1
                
                if ("ref" in qualitySum) :
                    qualitySum["ref"] += baseQ
                    mapQualitySum["ref"] += mapQ
                else :
                    qualitySum["ref"] = baseQ
                    mapQualitySum["ref"] = mapQ
                    
        elif (readBase.upper() == "A" or readBase.upper() == "G" or readBase.upper() == "C" or readBase.upper() == "T") :
            
            strand = "+"
            
            if (readBase == "a" or readBase == "g" or readBase == "c" or readBase == "t") :
                strand = "-"
                
            readBase = readBase.upper()
                
            if (baseQ >= minAvgQual) :
                if (readBase in readCounts) :
                    readCounts[readBase] += 1
                else :
                    readCounts[readBase] = 1
                    
                if (strand == "+") :
                    if (readBase in readCountsPlus) :
                        readCountsPlus[readBase] += 1
                    else :
                        readCountsPlus[readBase] = 1
                else :
                    if (readBase in readCountsMinus) :
                        readCountsMinus[readBase] += 1
                    else :
                        readCountsMinus[readBase] = 1
                        
                if (readBase in strandsSeen) :
                    alreadySeen = strandsSeen[readBase]
                    if (not(len(alreadySeen) >= 2 or alreadySeen == strand)) :
                        strandsSeen[readBase] += strand
                else :
                    strandsSeen[readBase] = strand
                    
                if (readBase in qualitySum) :
                    qualitySum[readBase] += baseQ
                    mapQualitySum[readBase] += mapQ
                else :
                    qualitySum[readBase] = baseQ
                    mapQualitySum[readBase] = mapQ
                    
            j += 1
            readStart = False
            
        elif (readBase == "+" or readBase == "-") :
            
            if (not readBases[i + 1].isdigit()) :
                i += 1
            else :
            
                indelType = ""
                
                if (readBase == "+") :
                    indelType = "INS"
                else :
                    indelType = "DEL"
                    
                if (prevBase == "." or prevBase == ",") :
                    if (prevBaseQ >= minAvgQual) :
                        reads1indel += 1
                        
                indelSize = 0
                maxParse = 1
                indelBases = ""
                
                try:
                
                    stringWithSize = readBases[i + 1] + readBases[i + 2] + readBases[i + 3]
                    stringWithSize = re.sub(r'[^0-9]', '', stringWithSize)
                    maxParse = int(indelSize) + len(indelSize)
                    
                    for basesParsed in range(0, maxParse) :
                        thisBase = readBases[i + 1 + basesParsed]
                        try :
                            catchNum = int(thisBase)
                        except Exception as ex:
                            if (thisBase == "." or thisBase == ",") :
                                basesParsed = maxParse
                            elif (thisBase.upper() == "A" or thisBase.upper() == "C" or thisBase.upper() == "G" or thisBase.upper() == "T" or thisBase.upper() == "N") :
                                indelBases = thisBase
                    i = i + maxParse
                except Exception as ex:
                    indelSize = int(readBases[i + 1])
                    for basesParsed in range(0, indelSize) :
                        indelBases += readBases[i + 2 + basesParsed]
                    i = i + 1 + indelSize
                    
                if (indelBases == indelBases.upper()) :
                    strand = "+"
                else:
                    strand = "-"
                    
                indelBases = indelBases.upper()
                
                indelKey = indelType + "-" + str(indelSize) + "-" + indelBases
                
                if (indelKey in readCounts) :
                    readCounts[indelKey] += 1
                else :
                    readCounts[indelKey] = 1
                    
                if (strand == "+") :
                    if (indelKey in readCountsPlus) :
                        readCountsPlus[indelKey] += 1
                    else :
                        readCountsPlus[indelKey] = 1
                else :
                    if (indelKey in readCountsMinus) :
                        readCountsMinus[indelKey] += 1
                    else :
                        readCountsMinus[indelKey] = 1
                
                if (indelKey in strandsSeen) :
                    alreadySeen = strandsSeen[indelKey]
                    if (not(len(alreadySeen) >= 2 or alreadySeen == strand)) :
                        strandsSeen[indelKey] += strand
                else :
                    strandsSeen[indelKey] = strand
                    
                if (j < len(readQuals)):
                    baseQ = ord(readQuals[j]) - 33
                    j += 1
                if (j < len(mapQuals)):
                    mapQ = ord(mapQuals[j]) - 33
                    
                if (indelKey in qualitySum) :
                    qualitySum[indelKey] += baseQ
                    mapQualitySum[indelKey] += mapQ
                else :
                    qualitySum[indelKey] = baseQ
                    mapQualitySum[indelKey] = mapQ
                    
                readStart = False
        elif (readBase.upper() == "N"):
            j += 1
        elif (readBase == "^") :
            i += 1
            readStart = True
        elif (readBase == "$") :
            readStart = False
        else :
            if (readBase == "." or readBase == ",") :
                pass
            else :
                j += 1
                    
    results = dict()
    
    strands1 = 0
    if ("ref" in strandsSeen) :
        strands1 = len(strandsSeen["ref"])
        
    avgQual1 = 0
    if (reads1 > 0) :
        avgQual1 = qualitySum["ref"] / reads1
        
    avgMapQual1 = 0
    if (reads1 > 0) :
        avgMapQual1 = mapQualitySum["ref"] / reads1
        
    reads1plus = 0
    reads1minus = 0
    if ("ref" in readCountsPlus):
       reads1plus = readCountsPlus["ref"]
    if ("ref" in readCountsMinus) :
        reads1minus = readCountsMinus["ref"]
        
    if (reads1 < 0) :
        reads1 = 0
        results[refBase] = str(reads1) + "\t" + str(strands1) + "\t" + str(avgMapQual1) + "\t" + str(avgMapQual1) + "\t" + str(reads1plus) + "\t" + str(reads1minus) + "\t" + str(reads1indel)
    
    variantKeys = readCounts.keys()
    variantKeys.sort()
    
    for key in variantKeys:
        reads2 = readCounts[key]
        
        reads2plus = 0
        reads2minus =0
        if (key in readCountsPlus) :
            reads2plus = readCountsPlus[key]
        if (key in readCountsMinus):
            reads2minus = readCountsMinus[key]

        strands2 = 0
        if (key in strandsSeen) :
            strands2 = len(strandsSeen[key])
            
        avgQual2 = qualitySum[key] / reads2
        
        avgmapQual2 = mapQualitySum[key] / reads2
        
        if (reads2 > 0) :
            results[key] = str(reads2) + "\t" + str(strands2) + "\t" + str(avgmapQual2) + "\t" + str(reads2plus) + "\t" + str(reads2minus) 
        
    return results

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

minCoverage = 10
minReads2 = 2
minAvgQual = 15
minVarFreq = 0.8
minFreqForHom = 0.75
pValueThreshold = 0.99
strandPvalueThreshold = 0.01
snpsOnly = True
strandFilter = False

verbose = True

numBases = 0
numVariantPositions = 0
numSNPpositions = 0
numIndelPositions = 0
numFailStrandFilter = 0
numVariantsReported = 0
numSNPsReported = 0
numIndelsReported = 0


numParsingExceptions = 0

vcfHeader = "##fileformat=VCFv4.1"
vcfHeader += "\n" + "##source=Slamdunk"
vcfHeader += "\n" + "##INFO=<ID=ADP,Number=1,Type=Integer,Description=\"Average per-sample depth of bases with Phred score >= " + str(minAvgQual) + "\">"
vcfHeader += "\n" + "##INFO=<ID=WT,Number=1,Type=Integer,Description=\"Number of samples called reference (wild-type)\">"
vcfHeader += "\n" + "##INFO=<ID=HET,Number=1,Type=Integer,Description=\"Number of samples called heterozygous-variant\">"
vcfHeader += "\n" + "##INFO=<ID=HOM,Number=1,Type=Integer,Description=\"Number of samples called homozygous-variant\">"
vcfHeader += "\n" + "##INFO=<ID=NC,Number=1,Type=Integer,Description=\"Number of samples not called\">"
vcfHeader += "\n" + "##FILTER=<ID=str10,Description=\"Less than 10% or more than 90% of variant supporting reads on one strand\">"
vcfHeader += "\n" + "##FILTER=<ID=indelError,Description=\"Likely artifact due to indel reads at this position\">"
vcfHeader += "\n" + "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">"
vcfHeader += "\n" + "##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype Quality\">"
vcfHeader += "\n" + "##FORMAT=<ID=SDP,Number=1,Type=Integer,Description=\"Raw Read Depth as reported by SAMtools\">"
vcfHeader += "\n" + "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Quality Read Depth of bases with Phred score >= " + str(minAvgQual) + "\">"
vcfHeader += "\n" + "##FORMAT=<ID=RD,Number=1,Type=Integer,Description=\"Depth of reference-supporting bases (reads1)\">"
vcfHeader += "\n" + "##FORMAT=<ID=AD,Number=1,Type=Integer,Description=\"Depth of variant-supporting bases (reads2)\">"
vcfHeader += "\n" + "##FORMAT=<ID=FREQ,Number=1,Type=String,Description=\"Variant allele frequency\">"
vcfHeader += "\n" + "##FORMAT=<ID=PVAL,Number=1,Type=String,Description=\"P-value from Fisher's Exact Test\">"
vcfHeader += "\n" + "##FORMAT=<ID=RBQ,Number=1,Type=Integer,Description=\"Average quality of reference-supporting bases (qual1)\">"
vcfHeader += "\n" + "##FORMAT=<ID=ABQ,Number=1,Type=Integer,Description=\"Average quality of variant-supporting bases (qual2)\">"
vcfHeader += "\n" + "##FORMAT=<ID=RDF,Number=1,Type=Integer,Description=\"Depth of reference-supporting bases on forward strand (reads1plus)\">"
vcfHeader += "\n" + "##FORMAT=<ID=RDR,Number=1,Type=Integer,Description=\"Depth of reference-supporting bases on reverse strand (reads1minus)\">"
vcfHeader += "\n" + "##FORMAT=<ID=ADF,Number=1,Type=Integer,Description=\"Depth of variant-supporting bases on forward strand (reads2plus)\">"
vcfHeader += "\n" + "##FORMAT=<ID=ADR,Number=1,Type=Integer,Description=\"Depth of variant-supporting bases on reverse strand (reads2minus)\">"
vcfHeader += "\n" + "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT"
#varscanCmd = "java -jar " + getBinary("VarScan.v2.4.1.jar") + " mpileup2snp  --strand-filter 0 --output-vcf --min-var-freq " + str(minVarFreq) + " --min-coverage " + str(minCov) + " --variants 1"

print(vcfHeader)

with open(args.pileup) as f:
    for line in f:
        
        fields = line.rstrip().split("\t")
        
        if (len(fields) > 5 and len(fields[0]) > 0 and len(fields[1]) > 0 and len(fields[2]) > 0 and len(fields[3]) > 0):
            chr = fields[0]
            pos = fields[1]
            ref = fields[2]
            depth = 0
            results = 0
            vcfResults = 0
            varAllel = 0
            snpFlag = False
            samplesRef = 0
            samplesHet = 0
            samplesHom = 0
            samplesUncalled = 0
            
            allReadDepth = 0
            allReads1plus = 0
            allReads1minus = 0
            allReads2plus = 0
            allReads2minus = 0
            strandPvalue = 1
            alReadBases = ""
            allReadQualities = ""
            
            depth = int(fields[3])
            readBases = fields[4]
            readQs = fields[5]
            mapQs = ""
            
            qualityDepth = 0
            for c in readQs :
                if (ord(c) - 33 >= minAvgQual) :
                    qualityDepth += 1
                    
            thisVcf = "./.:.:" + str(qualityDepth)
            
            if (depth >= minCoverage and qualityDepth >= minCoverage) :
                
                readCounts = getReadCounts(ref, readBases, readQs, minAvgQual, mapQs)
                positionCall = callPosition(ref, readCounts, "CNS", minReads2, minVarFreq, minAvgQual, pValueThreshold, minFreqForHom)
                #print(positionCall)
                                    
                    
        
        if (verbose and numBases % 100000 == 0 and numBases != 0) :
            print("Parsed " + str(numBases) + " positions.")
            
        numBases += 1

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
