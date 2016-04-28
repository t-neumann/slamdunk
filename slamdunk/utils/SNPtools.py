'''
Created on Jan 30, 2016

@author: philipp_
'''

from pybedtools import BedTool

class SNPDictionary(object):
    '''
    classdocs
    '''

    _vcfFile = None
    
    _tcSNPs = {}
    _agSNPs = {}

    def __init__(self, vcfFile):
        '''
        Constructor
        '''
        self._vcfFile = vcfFile

    def _addSNP(self, snp):
        if(snp[3].upper() == "T" and snp[4].upper() == "C"):
        #if(snp[3].upper() == "T" and snp[4].upper() == "G"):
            
            #key = snp[0] + (snp[1] - 1)
            key = snp[0] + snp[1]
            self._tcSNPs[key] = True
        
        if(snp[3].upper() == "A" and snp[4].upper() == "G"):
            #key = snp[0] + snp[1]
            key = snp[0] + (snp[1] - 1)
            self._agSNPs[key] = True
    
    #     if not snps is None:
#         with open(snps, "r") as f:
#             for line in f:
#                 if(len(line) > 0 and line[0] != "#"):
#                     # Parse VarScan VCF format
#                     cols = line.rstrip().split("\t")
#                     if(cols[3] == "A" and cols[4] == "G"):
#                         key = cols[0] + cols[1]
#                         snpDict[key] = 2
#                     if(cols[3] == 'T' and cols[4] == 'C'):
#                         key = cols[0] + cols[1]
#                         snpDict[key] = 1

        
    def read(self):
        vcfReader = BedTool(self._vcfFile)
    
        if(vcfReader.file_type != "vcf"):
            raise RuntimeError("Wrong file type. Not a vcf file.")
 
        for snp in vcfReader:
            self._addSNP(snp)
    
            
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


#x = SNPDictionary("/project/ngs/philipp/slamdunk-analysis/veronika/ngm-filtered/snps/34362_An312_wt-2n_mRNA-slamseq-autoquant_12h-R3.fq_slamdunk_mapped_filtered_snp.vcf")                   
#x.read()        