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

from scipy import stats
import numpy as np
from slamdunk.version import __version__

hypergeom = stats.distributions.hypergeom

baseAlphabet = ["A","G","C","T"]
fullAlphabet = ["A","G","C","T","N"]

def fisherFromHypergeom(expReadsRef, expReadsVar, obsReadsRef, obsReadsVar):
        
    if expReadsRef < 0 : expReadsRef = 0
    if expReadsVar < 0 : expReadsVar = 0
    if obsReadsRef < 0 : obsReadsRef = 0
    if obsReadsVar < 0 : obsReadsVar = 0
    
    n1 = obsReadsRef + obsReadsVar
    n2 = expReadsRef + expReadsVar
    n = obsReadsRef + expReadsRef
    
    rightTail = hypergeom.cdf(obsReadsVar, n1 + n2, n1, obsReadsVar + expReadsVar)
    leftTail = hypergeom.cdf(obsReadsRef, n1 + n2, n1, n)
    
    if (rightTail >= 0.999) :
        return leftTail
    else :
        return rightTail
    

def fisherTest(obsReadsRef, obsReadsVar) :
    
    pValue = 1.0
    
    seqError = 0.001
    
    coverage = obsReadsRef + obsReadsVar
    
    expReadsVar = int(float(coverage) * float(seqError))
    expReadsRef = coverage - expReadsVar
        
    pValue = fisherFromHypergeom(obsReadsRef, obsReadsVar, expReadsRef, expReadsVar)
    
    return pValue

def genotypeCoding(genotype) :
    
    if (genotype == "AA"): return "A"
    if (genotype == "CC"): return "C"
    if (genotype == "GG"): return "G"
    if (genotype == "TT"): return "T"
    if (genotype == "AC" or genotype == "CA"): return "M"
    if (genotype == "AG" or genotype == "GA"): return "R"
    if (genotype == "AT" or genotype == "TA"): return "W"
    if (genotype == "CG" or genotype == "GC"): return "S"
    if (genotype == "CT" or genotype == "TC"): return "Y"
    if (genotype == "GT" or genotype == "TG"): return "K"
    if (genotype[0] == "N") : return "N"
    
    return genotype
     
def callVariant(ref, variants, minVarFreq, minQualThreshold) :
    
    pValueThreshold = 0.01
    varReadThreshold = 2
    
    variantCall = ""
    
    refReads = 0
    varReads = 0
    refQual = 0
    varQual = 0
    pValue = 1
    varFreq = 0.00
    
        
    if (ref in variants) :
        
        try:
            refStats = variants[ref].split("\t")
            refReads = int(refStats[0])
            refQual = int(refStats[1])     
        except Exception as ex:
            print("Parsing error in positionCall: " + variants[ref], file=log)
            
    totalReadCounts = 0
    
    for allele in variants:
        alleleStats = variants[allele].split("\t")
        try :
            totalReadCounts += int(alleleStats[0])
        except Exception as ex:
            print("Parsing error in positionCall: " + variants[allele], file=log)
                    
    for allele in variants:
        alleleStats = variants[allele].split("\t")
        
        
        if (allele != ref) :
            thisRefReads      = refReads
            thisVarReads      = 0
            thisVarQual    = 0
            
            try:
                thisVarReads = int(alleleStats[0])
                thisVarQual = int(alleleStats[1])
                
            except Exception as ex:
                 print("Parsing error in positionCall: " + variants[allele], file=log)
                 
            
            if (thisVarReads > varReads) :
                thisVarFreq = float(thisVarReads) / float(totalReadCounts)
                thisPvalue = fisherTest(refReads, thisVarReads)
                
                if (thisVarReads > varReads and thisVarQual >= minQualThreshold) :
                        
                    varReads = thisVarReads
                    varQual = thisVarQual
                    varFreq = thisVarFreq
                    pValue = thisPvalue
                    
                if (thisVarReads >= varReadThreshold and thisVarQual >= minQualThreshold and thisVarFreq >= minVarFreq) :
                    thisRefReads = refReads
                    thisVarFreq = thisVarFreq
                                                                                            
                    if (thisPvalue <= pValueThreshold) :
                        
                        if (thisVarReads >= varReads):
                            
                            varReads      = thisVarReads
                            varQual    = thisVarQual
                            pValue      = thisPvalue

                            genotype = ref + allele

                            variantCall = genotypeCoding(genotype) + "\t" + str(thisRefReads) + "\t" + str(varReads) + "\t" 
                            variantCall += ('%.3f' % thisVarFreq) + "\t" + str(refQual) + "\t" + str(varQual) + "\t" 
                            variantCall += str(pValue) + "\t" + str(allele)
    
    return variantCall

def parsePileup(ref, readBases, readQuals, minQual):
    
    variants = dict()
    qualitySum = dict()
    
    refReads = 0
    readBase = ""
    nextBase = ""
    baseQ = 0
    
    qualPointer = 0
    basePointer = 0
    
    while basePointer < len(readBases) :
        readBase = readBases[basePointer]
            
        if (qualPointer < len(readQuals)):
            baseQ = ord(readQuals[qualPointer]) - 33
            
        if (readBase == "." or readBase == ",") :
            
            if (baseQ >= minQual) :
                refReads += 1
                
                if ("ref" in qualitySum) :
                    qualitySum["ref"] += baseQ
                else :
                    qualitySum["ref"] = baseQ
            qualPointer += 1
                    
        elif (readBase.upper() in baseAlphabet) :
                            
            readBase = readBase.upper()
                
            if (baseQ >= minQual) :
                if (readBase in variants) :
                    variants[readBase] += 1
                else :
                    variants[readBase] = 1
                    
                if (readBase in qualitySum) :
                    qualitySum[readBase] += baseQ
                else :
                    qualitySum[readBase] = baseQ
                    
            qualPointer += 1
            
        elif (readBase == "+" or readBase == "-") :
            
            basePointer += 1
            if (readBases[basePointer].isdigit()) :
                
                parse = ""
                while (readBases[basePointer].isdigit()) :
                    parse += readBases[basePointer]
                    basePointer += 1
                skip = int(parse) - 1
                
                while (skip > 0) :
                    basePointer += 1
                    skip -= 1

        elif (readBase == "^") :
            basePointer += 1
        elif (readBase == "$") :
            pass
        else:
            qualPointer += 1
                
        basePointer += 1
                    
    results = dict()
        
    refQual = 0
    if (refReads > 0) :
        refQual = qualitySum["ref"] / refReads
        
    if (refReads < 0) :
        refReads = 0
    results[ref] = str(refReads) + "\t" + str(refQual)
        
    variantKeys = variants.keys()
    variantKeys.sort()
    
    for key in variantKeys:
        varReads = variants[key]
            
        varQual = qualitySum[key] / varReads
                
        if (varReads > 0) :
            results[key] = str(varReads) + "\t" + str(varQual)
                       
    return results

def printVCFHeader(fileSNP, minQual):
    
    vcfHeader = "##fileformat=VCFv4.1"
    vcfHeader += "\n" + "##source=Slamdunk v" + __version__
    vcfHeader += "\n" + "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Raw Read Depth as reported by SAMtools\">"
    vcfHeader += "\n" + "##INFO=<ID=QDP,Number=1,Type=Integer,Description=\"Depth of bases with Phred score >= " + str(minQual) + "\">"
    vcfHeader += "\n" + "##INFO=<ID=PVAL,Number=1,Type=String,Description=\"P-value from Fisher's Exact Test\">"
    vcfHeader += "\n" + "##INFO=<ID=FREQ,Number=1,Type=String,Description=\"Variant allele frequency\">"
    vcfHeader += "\n" + "##INFO=<ID=RD,Number=1,Type=Integer,Description=\"Depth of reference-supporting bases\">"
    vcfHeader += "\n" + "##INFO=<ID=AD,Number=1,Type=Integer,Description=\"Depth of variant-supporting bases\">"
    vcfHeader += "\n" + "##INFO=<ID=RBQ,Number=1,Type=Integer,Description=\"Average quality of reference-supporting bases\">"
    vcfHeader += "\n" + "##INFO=<ID=ABQ,Number=1,Type=Integer,Description=\"Average quality of variant-supporting bases\">"
    vcfHeader += "\n" + "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT"
    
    print(vcfHeader, file=fileSNP)
    

def varCallPileup(mpileup, fileSNP, minVarFreq, minCoverage, minQual, log, verbose):
        
    baseCount = 0
    variantCount = 0
    
    printVCFHeader(fileSNP, minQual)
        
    with mpileup.stdout:
        for line in iter(mpileup.stdout.readline, b''):
        
            fields = line.rstrip().split("\t")
            
            # Check if proper pileup format
            
            if (len(fields) > 5):
                chr = fields[0]
                pos = fields[1]
                ref = fields[2].upper()
                variantFlag = False
                
                readDP = int(fields[3])
                readBases = fields[4]
                readQuals = fields[5]
                
                varBase = "."
                refReads = 0
                varReads = 0
                varFreq = 0
                refQual = 0
                varQual= 0
                pValue = 0
                
                qualityDP = 0
                
                for c in readQuals :
                    if (ord(c) - 33 >= minQual) :
                        qualityDP += 1
                
                if (qualityDP >= minCoverage) :
                    
                    variants = parsePileup(ref, readBases, readQuals, minQual)
                    variantCall = callVariant(ref, variants, minVarFreq, minQual)
                             
                    if (len(variantCall) > 0) :
                        callLines = variantCall.split("\n")
                        
                        for lineCounter in range(0, len(callLines)):
                            callContents = callLines[lineCounter].split("\t")
                            consBase = callContents[0]
                            refReads = int(callContents[1])
                            varReads = int(callContents[2])
                            varFreq = callContents[3]
                            refQual = int(callContents[4])
                            varQual= int(callContents[5])
                            pValue = float(callContents[6])
                            
                            if (consBase != ref and consBase != "N" and len(callContents) > 7) :
                                varBase = callContents[7]             
                                variantFlag = True
                    
                outLine = chr + "\t" + pos + "\t"            
                                   
                outLine += "." + "\t" + ref.upper() + "\t" + varBase + "\t.\tPASS\t";
                outLine += "DP=" + str(readDP) + ";QDP=" + str(qualityDP) + ";PVAL=" + str(pValue) + ";FREQ=" + str(varFreq)
                outLine += ";RD=" + str(refReads) + ";AD=" + str(varReads) + ";RBQ=" + str(refQual) + ";ABQ=" + str(varQual)
                                    
                if variantFlag:
                    print(outLine, file=fileSNP)
                    variantCount += 1
                        
            else:
                if (len(fields) >= 4 and fields[3] == "0") :
                    pass
                else :
                    print("Warning: Line ignored: Invalid pileup format in line " + str(baseCount) + "\n" + line + "\n", file=log)    
       
            
            if (verbose and baseCount % 100000 == 0 and baseCount != 0) :
                print("Parsed " + str(baseCount) + " positions.", file=log)
                
            baseCount += 1
    mpileup.wait()
      
    if (verbose):      
        print(str(baseCount) + " bases in pileup file", file=log)
        print(str(variantCount) + " variant positions.", file=log)

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