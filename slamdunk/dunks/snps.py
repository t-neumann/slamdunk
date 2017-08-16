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
import subprocess
import csv, re
from slamdunk.utils.misc import checkStep, getBinary  # @UnresolvedImport

from fisher import pvalue
import numpy as np
from slamdunk.version import __version__

def getQualityDepth(readQuals, minAvgQual):
    
    qualityDepth = 0;
    
    j = 0
    
    while (j < len(readQuals)) :
        
        baseQuality = ord(readQuals[j]) - 33
        if (baseQuality >= minAvgQual) :
            qualityDepth += 1
        j += 1

    return qualityDepth

def isHomozygous(genotype):
    
    if "/" in genotype:
        alleles = genotype.split("/")
        if (alleles[0] == alleles[1]):
                 return True
    else:
        if (genotype == "A" or genotype == "C" or genotype == "G" or genotype == "T") :
            return True
    return False

def getFisher(expReads1, expReads2, obsReads1, obsReads2):
        
    if expReads1 < 0 : expReads1 = 0
    if expReads2 < 0 : expReads2 = 0
    if obsReads1 < 0 : obsReads1 = 0
    if obsReads2 < 0 : obsReads2 = 0
    
    p = pvalue(obsReads1, obsReads2, expReads1, expReads2)
    
    if (p.right_tail >= 0.999) :
        return p.left_tail
    else :
        return p.right_tail
    

def getSignificance(obsReads1, obsReads2) :
    
    pValue = 1.0
    baselineError = 0.001
    
    coverage = obsReads1 + obsReads2
    
    expReads2 = int(float(coverage) * float(baselineError))
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
            print("Warning: error generating consensus from " + gt, file=log)
            
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
        
        if (refBase in readCounts) :
            
            try:
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
                print("refBase readcount error: " + readCounts[refBase], file=log)

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
                    
                    if (thisReads2 > reads2 and thisAvgQual2 >= minAvgQual) :
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
                            
                            if (callType == "CNS" and thisReads2 >= reads2):
                                reads2      = thisReads2
                                strands2    = thisStrands2
                                avgQual2    = thisAvgQual2
                                avgMap2     = thisAvgMap2
                                reads2plus  = thisReads2plus
                                reads2minus = thisReads2minus
                                pValue      = thisPvalue

                                genotype = ""
                                
                                if (thisVarFreq >= (minFreqForHom * 100)):
                                    genotype = allele + allele
                                    if thisVarType == "INDEL":
                                        genotype = allele + "/" + allele
                                else:
                                    genotype = refBase + allele
                                    if thisVarType == "INDEL":
                                        genotype = "*/" + allele

                                callResult = genotypeToCode(genotype) + "\t" + str(thisReads1) + "\t" + str(reads2) + "\t" + ('%05.3f' % thisVarFreq) + "%\t" + str(strands1) + "\t"
                                callResult += str(strands2) + "\t" + str(avgQual1) + "\t" + str(avgQual2) + "\t" + str(pValue) + "\t" + str(avgMap1) + "\t" + str(avgMap2)
                                callResult += "\t" + str(reads1plus) + "\t" + str(reads1minus) + "\t" + str(reads2plus) + "\t" + str(reads2minus) + "\t" + str(varAllele)
    except Exception as ex:
        pass
    
    if (len(callResult) == 0 and callType == "CNS") :
        if (reads1 > 0 and reads1 > minReads2) :
            callResult = str(refBase) + "\t" + str(reads1) + "\t" + str(reads2) + "\t" + ('%05.3f' % varFreq) + "%\t" + str(strands1) + "\t" + str(strands2) + "\t" 
            callResult += str(avgQual1) + "\t" + str(avgQual2) + "\t" + str(pValue) + "\t" + str(avgMap1) + "\t" + str(avgMap2)
            callResult += "\t" + str(reads1plus) + "\t" + str(reads1minus) + "\t" + str(reads2plus) + "\t" + str(reads2minus) + "\t" + str(varAllele)
        else:
            callResult = "N" + "\t" + str(reads1) + "\t" + str(reads2) + "\t" + ('%05.3f' % varFreq) + "%\t" + str(strands1) + "\t" + str(strands2)
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
    
    i = 0
    while i < len(readBases) :
        readBase = readBases[i]
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
            
        if ((readBase == "." or readBase == ",") and nextBase != "-" and nextBase != "+") :
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
            j += 1
            readStart = False
                    
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
                    indelSize = int(stringWithSize);
                    maxParse = int(stringWithSize) + len(stringWithSize)
                            
                    for basesParsed in range(0, maxParse) :
                        thisBase = readBases[i + 1 + basesParsed]
                        try :
                            catchNum = int(thisBase)
                        except Exception as ex:
                            if (thisBase == "." or thisBase == ",") :
                                basesParsed = maxParse
                            elif (thisBase.upper() == "A" or thisBase.upper() == "C" or thisBase.upper() == "G" or thisBase.upper() == "T" or thisBase.upper() == "N") :
                                indelBases += thisBase

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
                
        i += 1
                    
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
    results[refBase] = str(reads1) + "\t" + str(strands1) + "\t" + str(avgQual1) + "\t" + str(avgMapQual1) + "\t" + str(reads1plus) + "\t" + str(reads1minus) + "\t" + str(reads1indel)
        
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
            results[key] = str(reads2) + "\t" + str(strands2) + "\t" + str(avgQual2) + "\t" + str(avgmapQual2) + "\t" + str(reads2plus) + "\t" + str(reads2minus) 
        
    return results

def varCallPileup(mpileup, fileSNP, minVarFreq, minCoverage, minAvgQual, log, verbose):
    
    minReads2 = 2
    minFreqForHom = 0.75
    pValueThreshold = 0.01
    strandPvalueThreshold = 0.01
    snpsOnly = True
    indelsOnly = False
    variantsOnly = False
    strandFilter = False
        
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
    vcfHeader += "\n" + "##source=Slamdunk_v" + __version__
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
    
    print(vcfHeader, file=fileSNP)
    
    with mpileup.stdout:
        for line in iter(mpileup.stdout.readline, b''):
        
            fields = line.rstrip().split("\t")
            
            if (len(fields) > 5 and len(fields[0]) > 0 and len(fields[1]) > 0 and len(fields[2]) > 0 and len(fields[3]) > 0):
                chr = fields[0]
                pos = fields[1]
                refBase = fields[2].upper()
                callDepths = ""
                callResults = ""
                vcfResults = ""
                varAlleles = dict()
                variantFlag = False
                snpFlag = False
                indelFlag = False
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
                allReadBases = ""
                allReadQualities = ""
                
                readDepth = int(fields[3])
                readBases = fields[4]
                readQualities = fields[5]
                mapQualities = ""
                
                allReadDepth += readDepth;
                allReadBases = allReadBases + readBases;
                allReadQualities = allReadQualities + readQualities;
                
                qualityDepth = 0
                for c in readQualities :
                    if (ord(c) - 33 >= minAvgQual) :
                        qualityDepth += 1
                        
                thisCall = "N" + ":" + str(qualityDepth) + ":-:-:-:-"                    
                thisVCF = "./.:.:" + str(qualityDepth)
                
                if (readDepth >= minCoverage and qualityDepth >= minCoverage) :
                    
                    readCounts = getReadCounts(refBase, readBases, readQualities, minAvgQual, mapQualities)
                    positionCall = callPosition(refBase, readCounts, "CNS", minReads2, minVarFreq, minAvgQual, pValueThreshold, minFreqForHom)
                                    
                    if (len(positionCall) > 0) :
                        callLines = positionCall.split("\n")
                        
                        for lineCounter in range(0, len(callLines)):
                            callContents = callLines[lineCounter].split("\t")
                            consBase = callContents[0]
                            reads1 = int(callContents[1])
                            reads2 = int(callContents[2])
                            varFreq = callContents[3]
                            strands1 = int(callContents[4])
                            strands2 = int(callContents[5])
                            qual1 = int(callContents[6])
                            qual2= int(callContents[7])
                            pValue = float(callContents[8])
                            reads1plus = int(callContents[11])
                            reads1minus = int(callContents[12])
                            reads2plus = int(callContents[13])
                            reads2minus = int(callContents[14])
                            varAllele = ""
                            
                            logP = 0.0
                            if (pValue > 0) :
                                logP = 0 - (10 * np.log10(pValue))
                                if (logP > 255) : logP = 255
                            
                            if (consBase != refBase and consBase != "N" and len(callContents) > 15) :
                                varAllele = callContents[15]
                                
                                varAlleleNumber = 0
                                
                                if (varAllele in varAlleles) :
                                    varAlleleNumber = varAlleles[varAllele]
                                else :
                                    varAlleleNumber = len(varAlleles) + 1
                                    varAlleles[varAllele] = varAlleleNumber
                                
                                if (isHomozygous(consBase)) :
                                    samplesHom += 1
                                    thisVCF = str(varAlleleNumber) + "/" + str(varAlleleNumber)
                                else :
                                    samplesHet += 1
                                    thisVCF = "0" + "/" + str(varAlleleNumber)
                                    
                                thisVCF += ":" + str(int(logP)) + ":" + str(readDepth) + ":" + str(qualityDepth)
                                thisVCF += ":" + str(reads1) + ":" + str(reads2) + ":" + str(varFreq) + ":" + str(pValue)
                                thisVCF += ":" + str(qual1) + ":" + str(qual2)
                                thisVCF += ":" + str(reads1plus) + ":" + str(reads1minus) + ":" + str(reads2plus) + ":" + str(reads2minus)
                                
                            elif (consBase == refBase) :
                                expReads1 = int((reads1 + reads2) / 2)
                                expReads2 = (reads1 + reads2) - expReads1
                                
                                newPvalue = getFisher(reads1, reads2, expReads1, expReads2)
                                newLogP = 0
                                if (newPvalue > 0) :
                                    newLogP = 0 - (10 * np.log10(newPvalue))
                                thisVCF = "0" + "/" + "0"
                                thisVCF += ":" + str(int(newLogP)) + ":" + str(readDepth) + ":" + str(qualityDepth)
                                thisVCF += ":" + str(reads1) + ":" + str(reads2) + ":" + str(varFreq) + ":" + str(pValue)
                                thisVCF += ":" + str(qual1) + ":" + str(qual2)
                                thisVCF += ":" + str(reads1plus) + ":" + str(reads1minus) + ":" + str(reads2plus) + ":" + str(reads2minus)
                                
                            thisCall = consBase + ":" + str(qualityDepth) + ":" + str(reads1) + ":" + str(reads2)
                            thisCall += ":" + str(varFreq) + ":" + str(pValue)
                            
                            if (consBase != refBase and consBase != "N") :
                                variantFlag = True
                                if (len(consBase) > 1) :
                                    indelFlag = True
                                else :
                                    snpFlag = True
                                    
                                allReads1plus += reads1plus
                                allReads1minus += reads1minus
                                allReads2plus += reads2plus
                                allReads2minus += reads2minus
                                
                            else :
                                samplesRef += 1
                    
                    else :
                        samplesUncalled += 1
                else :
                    samplesUncalled += 1
                                
                if len(callDepths) > 0:
                    callDepths += " "
                callDepths += str(readDepth)
                
                if (len(callResults) > 0) :
                    callResults += " "
                callResults += thisCall
                
                if (len(vcfResults) > 0) :
                    vcfResults += "\t"
                vcfResults += thisVCF
                
                qualityDepth = 0
                qualityDepth = getQualityDepth(allReadQualities, minAvgQual)
                
                allMapQualities = ""
                allConsensusCall = "N:" + str(qualityDepth) + ":-:-:-:-"
                
                varBases = ""
                sortedKeys = varAlleles.keys()
                alleleKeys = varAlleles.keys()
                
                for allele in sortedKeys:
                    arrayIndex = varAlleles[allele] - 1
                    alleleKeys[arrayIndex] = allele
                    
                for allele in alleleKeys :
                    if (len(varBases) > 0) :
                        varBases += ","
                    varBases += allele
                    
                if(len(varBases) == 0) :
                    varBases = "."
                    
                if (variantFlag) :
                    numVariantPositions += 1
                if (snpFlag) :
                    numSNPpositions += 1
                if (indelFlag) :
                    numIndelPositions += 1
                    
                strandFilterStatus = "Pass:" + str(allReads1plus) + ":" + str(allReads1minus) + ":" + str(allReads2plus) + ":" + str(allReads2minus) + ":" + str(strandPvalue)
                failedStrandFilter = False
                    
                outLine = chr + "\t" + pos + "\t"
                
                avgQualityDepth = 0
                
                if (samplesRef + samplesHet + samplesHom + samplesUncalled > 0) :
                    avgQualityDepth = qualityDepth / (samplesRef + samplesHet + samplesHom + samplesUncalled)
                    
                refColumn = ""
                varColumn = ""
                
                if ("," in varBases and "-" in varBases or "+" in varBases) :
                    maxDelSize = 0
                    maxDelBases = ""
                    varBaseContents = varBases.split(",")
                    for varAllele in varBaseContents:
                        if varAllele[0] == "-":
                            varAllele = re.sub(r'-', '', varAllele)
                            if (len(varAllele) > maxDelSize) :
                                maxDelBases = varAllele
                                maxDelSize = len(varAllele)
                
                
                    refColumn = refBase + maxDelBases
                    
                    varColumn = ""
                    
                    for varAllele in varBaseContents:
                        if len(varColumn) > 0 :
                            varColumn += ","
                            
                        if varAllele[0] == "-":
                            varAllele = re.sub(r'-', '', varAllele)
                        
                            if (len(varAllele) < maxDelSize) :
                                varEntry = re.sub(varAllele, '', varEntry)
                                varColumn = varColumn + refBase + varEntry
                            else :
                                varColumn += refBase
                                
                        elif varAllele[0] == "+":
                            varAllele = re.sub(r'\+', '', varAllele)
                            varEntry = refBase + varAllele + maxDelBases
                            varColumn += varEntry
                        else:
                            varEntry = varAllele + maxDelBases
                            varColumn += varEntry
                            
                elif (varBases[0] == "+") :
                    refColumn = refBase
                    varColumn = refBase + re.sub(r'\+', '', varBases)
                    
                elif (varBases[0] == "-") :
                    refColumn = refBase + re.sub('-', '', varBases)
                    varColumn = refBase
                    
                else:
                    refColumn = refBase
                    varColumn = varBases
                    
                varColumn = re.sub(r'\+', '', varColumn)
                varColumn = re.sub(r'-', '', varColumn)
                
                outLine += "." + "\t" + refColumn.upper() + "\t" + varColumn + "\t.\t";
                
                if ("Pass" in strandFilterStatus):
                    outLine += "PASS\t"
                else :
                    outLine += "str10\t"
                outLine += "ADP=" + str(avgQualityDepth) + ";WT=" + str(samplesRef) + ";HET=" + str(samplesHet) + ";HOM=" + str(samplesHom) + ";NC=" + str(samplesUncalled)
                outLine += "\t" + "GT:GQ:SDP:DP:RD:AD:FREQ:PVAL:RBQ:ABQ:RDF:RDR:ADF:ADR" + "\t" + vcfResults
                    
                reportFlag = False
                
                if (variantFlag and strandFilter and failedStrandFilter) :
                    if (not variantsOnly and not snpsOnly and not indelsOnly) :
                        reportFlag = True
                elif ((variantsOnly or snpsOnly or indelsOnly) and not variantFlag) :
                    pass
                elif (not variantsOnly and not snpsOnly and not indelsOnly) :
                    reportFlag = True
                elif (variantFlag and variantsOnly) :
                    reportFlag = True
                elif (snpFlag and snpsOnly) :
                    reportFlag = True
                elif (indelFlag and indelsOnly) :
                    reportFlag = True
                else :
                    pass
                            
                if reportFlag:
                    print(outLine, file=fileSNP)
                    
                    if variantFlag:
                        numVariantsReported += 1
                    if snpFlag:
                        numSNPsReported += 1
                    if indelFlag:
                        numIndelsReported += 1
                        
            else:
                if (len(fields) >= 4 and fields[3] == "0") :
                    pass
                else :
                    print("Warning: Line ignored: Invalid pileup format in line " + str(numBases) + "\n" + line + "\n", file=log)    
       
            
            if (verbose and numBases % 100000 == 0 and numBases != 0) :
                print("Parsed " + str(numBases) + " positions.", file=log)
                
            numBases += 1
    mpileup.wait()
      
    if (verbose):      
        print(str(numBases) + " bases in pileup file", file=log)
        print(str(numVariantPositions) + " variant positions (" + str(numSNPpositions) + " SNP, " + str(numIndelPositions) + " indel)", file=log)
        print(str(numFailStrandFilter) + " were failed by the strand-filter", file=log)
        print(str(numVariantsReported) + " variant positions reported (" + str(numSNPsReported) + " SNP, " + str(numIndelsReported) + " indel)", file=log)

def SNPs(inputBAM, outputSNP, referenceFile, minVarFreq, minCov, minQual, log, printOnly=False, verbose=True, force=False):
    if(checkStep([inputBAM, referenceFile], [outputSNP], force)):
        fileSNP = open(outputSNP, 'w')
        
        mpileupCmd = getBinary("samtools") + " mpileup -B -A -f " + referenceFile + " " + inputBAM
        if(verbose):
            print(mpileupCmd, file=log)
        if(not printOnly):
            mpileup = subprocess.Popen(mpileupCmd, shell=True, stdout=subprocess.PIPE, stderr=log)
           
        if (not printOnly) :
            varCallPileup(mpileup, fileSNP, minVarFreq, minCov, minQual, log, verbose) 
#         varscanCmd = "java -jar " + getBinary("VarScan.v2.4.1.jar") + " mpileup2snp  --strand-filter 0 --output-vcf --min-var-freq " + str(minVarFreq) + " --min-coverage " + str(minCov) + " --variants 1"
#         if(verbose):
#             print(varscanCmd, file=log)
#         if(not printOnly):
#             varscan = subprocess.Popen(varscanCmd, shell=True, stdin=mpileup.stdout, stdout=fileSNP, stderr=log)
#             varscan.wait()
        
        fileSNP.close()
    else:
        print("Skipping SNP calling", file=log)    
        
def countSNPsInFile(inputFile):
    snpCount = 0
    tcSnpCount = 0
    with open(inputFile, "r") as snpFile:
            snpReader = csv.reader(snpFile, delimiter='\t')
            for row in snpReader:
                if((row[2].upper() == "T" and row[3].upper() == "C") or (row[2].upper() == "A" and row[3].upper() == "G")):
                    tcSnpCount = tcSnpCount + 1
                snpCount = snpCount + 1
    return snpCount, tcSnpCount