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
            key = snp[0] + snp[1]
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
