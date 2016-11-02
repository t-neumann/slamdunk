'''
Created on Jan 30, 2016

@author: philipp_
'''

import os

from pybedtools import BedTool


class SNPDictionary(object):
    
    def __init__(self, vcfFile):
        self._vcfFile = vcfFile
        self._tcSNPs = {}
        self._agSNPs = {}

    def _addSNP(self, snp):
        if(snp[3].upper() == "T" and snp[4].upper() == "C"):
            key = snp[0] + snp[1]
            self._tcSNPs[key] = True
        
        if(snp[3].upper() == "A" and snp[4].upper() == "G"):
            key = snp[0] + str(int(snp[1]) - 1)
            self._agSNPs[key] = True
        
    def read(self):        
        if (self._vcfFile != None):
            if(os.path.exists(self._vcfFile)):
                vcfReader = BedTool(self._vcfFile)
                
                if(vcfReader.file_type != "vcf"):
                    print("Wrong file type. Empty or not a vcf file.")
         
                for snp in vcfReader:
                    self._addSNP(snp)
            else:
                print("Warning: SNP file " + self._vcfFile + " not found.")
            
    def isAGSnp(self, chromosome, position):
        key = chromosome + str(int(position) + 1)
        return key in self._agSNPs
    
    
    def isTCSnp(self, chromosome, position):
        key = chromosome + str(int(position) + 1)
        return key in self._tcSNPs

    def getAGSNPsInUTR(self, chromosome, start, stop, snpType):
        count = 0
        for i in range(start, stop):
                if(self.isAGSnp(chromosome, i)):
                    count += 1
        return count

    def getTCSNPsInUTR(self, chromosome, start, stop, snpType):
        count = 0
        for i in range(start, stop):
                if(self.isTCSnp(chromosome, i)):
                    count += 1
        return count
