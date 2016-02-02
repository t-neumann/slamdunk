'''
Created on Jan 29, 2016

@author: philipp_
'''
import pysam
import re

class ReadDirection:
    Forward = 1
    Reverse = 2

class SlamSeqRead:
    
    readName = None
    tcCount = None
    conversionRates = None
    direction = None 
    
#     def getRate(rates, refBase, readBase):
#         return rates[baseNumber * encodeBase(refBase) + encodeBase(readBase)]

    
class SlamSeqIterator:
    
    _readIterator = None
#     _refSeq = None
    _snps = None
    _maxReadLength = 0
    _minQual = 0
    _chromosome = None
    _startPosition = 0
    
    _baseNumber = 5
    _toBase = [ 'A', 'C', 'G', 'T', 'N' ]

    
    def encodeBase(self, base):
        if(base.upper() == 'A'):
            return 0
        if(base.upper() == 'C'):
            return 1
        if(base.upper() == 'G'):
            return 2
        if(base.upper() == 'T'):
            return 3
        
        return 4

    
    def incRate(self, rates, refBase, readBase):
        rates[self._baseNumber * self.encodeBase(refBase) + self.encodeBase(readBase)] += 1
        return rates

    
    def computeRatesForRead(self, read):
        rates = [0] * 25
    
        for pair in read.get_aligned_pairs(matches_only=True, with_seq=False):
            readPos = pair[0]
#             refPos = pair[1]
            refPos = pair[1] - int(self._startPosition) + self._maxReadLength + 1
            refBase = self._refSeq[refPos]
#             refBase = pair[2]
            readBase = read.query_sequence[readPos]
            readQlty = read.query_qualities[readPos]
#             print(readPos, refPos, readBase, refBase)
            if(readQlty >= self._minQual):            
                rates = self.incRate(rates, refBase, readBase)
                        
        return rates

    def computeRatesForReadNGM(self, read):
        ratesNgm = None
        if(read.has_tag("RA")):
            ratesNgm = map(int, read.get_tag("RA").split(","))
        return ratesNgm
        
    
    def getTCNgm(self, read):
        return int(read.get_tag("TC"))
    
    def getTC(self, read):
        tcCount = 0;
        for pair in read.get_aligned_pairs(matches_only=True, with_seq=False):
            readPos = pair[0]
            refPos = pair[1] - int(self._startPosition) + self._maxReadLength + 1
#             refBase = pair[2]
#             print(len(self._refSeq), refPos)
            refBase = self._refSeq[refPos]
            readBase = read.query_sequence[readPos]
            readQlty = read.query_qualities[readPos]
            if readBase != refBase:
    #             print(pair, readPos, refPos, readBase, refBase, readQlty, file=sys.stderr)
                if(readQlty >= self._minQual):
                    if(read.is_reverse):
                        if(refBase == 'A' and readBase == 'G'):
                            if(self._snps == None or not self._snps.isAGSnp(self._chromosome, int(pair[1]))):
                                tcCount += 1
    #                         print(pair, readPos, refPos, readBase, refBase, readQlty, file=sys.stderr)
                    else:
                        if(refBase == 'T' and readBase == 'C'):
                            if(self._snps == None or not self._snps.isTCSnp(self._chromosome, int(pair[1]))):
                                tcCount += 1
    #                         print(pair, readPos, refPos, readBase, refBase, readQlty, file=sys.stderr)
                        
        return tcCount

    
    def __init__(self, readIterator, refSeq, chromosome, startPosition, maxReadLength, snps):
        self._readIterator = readIterator
        self._refSeq = refSeq
        self._chromosome = chromosome
        self._startPosition = startPosition
        self._maxReadLength = maxReadLength
        self._snps = snps
        
    def __iter__(self):
        return self
 
    #Check if two rates arrays are equal
    def compareLists(self, a, b):
        if(len(a) != len(b)):
            return False
        for x,y in zip(a, b):
            if(a != b):
                return False
    
        return True

 
    def next(self):
        read = self._readIterator.next()
        try:
            slamSeqRead = SlamSeqRead()
            slamSeqRead.readName = read.query_name
            if(read.is_reverse):
                slamSeqRead.direction = ReadDirection.Reverse
            else:
                slamSeqRead.direction = ReadDirection.Forward
            ngmTC = self.getTC(read)
            slamSeqRead.tcCount = self.getTC(read)
            
            slamSeqRead.conversionRates = self.computeRatesForRead(read)
            ngmRates = self.computeRatesForReadNGM(read)    
            
            
            if(not self.compareLists(slamSeqRead.conversionRates, ngmRates) or slamSeqRead.tcCount != ngmTC):
                print("Difference found:")
                print(read)
                print(ngmRates)
                print(slamSeqRead.conversionRates)
                print("TC (ngm): " + str(ngmTC))
                print("TC (pys): " + str(slamSeqRead.tcCount))
                #sys.stdin.read(1)
                raise RuntimeError("Difference found between NGM and Py.")
            
            return slamSeqRead
        
        except IndexError as e:
            #Error is: IndexError: string index out of range
            #TODO: use with_seq=False for get_aligned_pairs instead of reading ref sequence manually
            print("Error computing rates for read " + read.query_name)
            print(e)
            print(read)
            return self.next()
        
        
class SlamSeqFile:
    '''
    classdocs
    '''
    _bamFile = pysam.AlignmentFile
    _referenceFile = pysam.FastaFile
    _snps = None

    def __init__(self, bamFile, referenceFile, snps):
        self._bamFile = pysam.AlignmentFile(bamFile, "rb")
        self._referenceFile = pysam.FastaFile(referenceFile)   
        self._snps = snps
        
    def readInRegion(self, chromosome, start, stop, maxReadLength):
        refRegion = chromosome + ":" + str(int(start) - maxReadLength) + "-" + str(int(stop) + maxReadLength)
        
        region = chromosome + ":" + start + "-" + stop
        refSeq = self._referenceFile.fetch(region=refRegion)
        return SlamSeqIterator(self._bamFile.fetch(region=region), refSeq, chromosome, start, maxReadLength, self._snps)
    
    def readsInChromosome(self, chromosome):
        refSeq = self._referenceFile.fetch(region=chromosome)
        return SlamSeqIterator(self._bamFile.fetch(region=chromosome), refSeq, chromosome, 1, 0, self._snps)
    
    
    def atoi(self, text):
        return int(text) if text.isdigit() else text
    
    def natural_keys(self, text):
        '''
        alist.sort(key=natural_keys) sorts in human order
        http://nedbatchelder.com/blog/200712/human_sorting.html
        (See Toothy's implementation in the comments)
        '''
        return [ self.atoi(c) for c in re.split('(\d+)', text) ]

    def getChromosomes(self):
        refs = list(self._referenceFile.references)
        refs.sort(key=self.natural_keys)
        return refs

#     def __iter__(self):
#         return self
# 
#     def next(self):
#         if self.i < self.n:
#             i = self.i
#             self.i += 1
#             return i
#         else:
#             raise StopIteration()
        
# snps = SNPtools.SNPDictionary("/project/ngs/philipp/slamdunk-analysis/debug/snps/ngm/26338_mESC-wt_0.5h-4SU_trimmed_fixed_downsample_slamdunk_mapped_filtered_snp.vcf")
# testFile = SlamSeqFile("/project/ngs/philipp/slamdunk-analysis/debug/filtered/ngm/26338_mESC-wt_0.5h-4SU_trimmed_fixed_downsample_slamdunk_mapped_filtered.bam", "/project/ngs/philipp/slamseq/ref/GRCm38.fa", 55, snps)
# 
# chromosomes = testFile.getChromosomes()
# 
# for chromosome in chromosomes:
#     print(chromosome)
#     readIterator = testFile.readsInChromosome(chromosome)
# 
#     for read in readIterator:
#         print(read.readName, read.direction, read.tcCount, read.conversionRates)

