from __future__ import print_function
import pysam
import re
import sys

class ReadDirection:
    Forward = 1
    Reverse = 2

#               Read
#             A     C     G     T     N
#      A      0     1     2     3     4
# R    C      5     6     7     8     9
# e    G     10    11    12    13    14
# f    T     15    16    17    18    19
#      N     20    21    22    23    24
class SlamSeqConversionRates:

    _baseNumber = 5
    _toBase = [ 'A', 'C', 'G', 'T', 'N' ]
    
    _data = []
    
    # Make the object act like a list
    def __len__(self):
        return len(self._data)
    
    def __getitem__(self, index):
        return self._data[index]

    def __repr__(self):
        return ",".join(str(x) for x in self._data)
    
    def __iter__(self):
        return self._data.__iter__()

    def __init__(self):
        self._data = [0] * (self._baseNumber * self._baseNumber)
         
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
    
    def incRate(self, refBase, readBase):
        self._data[self._baseNumber * self.encodeBase(refBase) + self.encodeBase(readBase)] += 1

    def getRate(self, refBase, readBase):
        return self._data[self._baseNumber * self._encodeBase(refBase) + self._encodeBase(readBase)]
    

class SlamSeqAlignmentPosition:
    # Position on the read. 
    # ReadPos 0 on a reverse read is the last position of the read reported 
    # in the BAM file
    readPosition = None
    # Position in the reference sequence were the mismatch occurred
    referencePosition = None
    # Base in read
    readBase = None
    # Bas in reference
    referenceBase = None
    # Base quality
    readBaseQlty = None
    # If the mismatch overlaps with a SNP position (if vcf file is specified)
    isSnpPosition = None
    
    def __init__(self, readPosition, referencePosition, readBase, referenceBase, readBaseQlty, isSnpPos):
        self.readBase = readBase
        self.readBaseQlty = readBaseQlty
        self.readPosition = readPosition
        self.referenceBase = referenceBase
        self.referencePosition = referencePosition
        self.isSnpPosition = isSnpPos
        
    def isMismatch(self):
        return self.readBase != self.referenceBase
    
    def __repr__(self):
        return str(self.referencePosition) + "," + self.referenceBase + "," + str(self.readPosition) + "," + self.readBase + "," + str(self.readBaseQlty) + "," +  str(self.isSnpPosition)
        #return self.referenceBase + "," + str(self.readPosition) + "," + self.readBase + "," + str(self.readBaseQlty) + "," +  str(self.isSnpPosition)
    
    def isTCMismatch(self, isReverse):
        if(isReverse):
            return self.referenceBase == "A" and self.readBase == "G" and not self.isSnpPosition
        else:
            return self.referenceBase == "T" and self.readBase == "C" and not self.isSnpPosition
        
    def isT(self, isReverse):
        if(isReverse):
            return self.referenceBase == "A"
        else:
            return self.referenceBase == "T"
        

class SlamSeqRead:
    
    # Name of the parsed read
    name = None
    # Nuber of Ts that were converted 
    # to a C on forward reads and A to G on revse reads
    tcCount = None
    # Number Ts in the reference
    tCount = None
    # Percentage of converted Ts/As
    tcRate = None
    # Number of all possible conversions for one reads
    conversionRates = None
    # Direction of the reads
    direction = None 
    # Read sequence
    sequence = None
    # List of mismatches in alignment
    mismatches = []
    
    def __repr__(self):
        return "\t".join([self.name, str(self.direction), self.sequence, str(self.tcCount), str(self.tCount), str(self.tcRate), self.conversionRates.__repr__(), self.mismatches.__repr__()])

class SlamSeqWriter:
    
    _seperator = '\t'
    _file = None
    
    def __init__(self, fileName):
        self._file = open(fileName, "w")
        self._printHeader()
    
    def _printHeader(self):
        print("Name", file=self._file, end=self._seperator)
        print("Direction",  file=self._file, end=self._seperator)
        print("Sequence",  file=self._file, end=self._seperator)
        print("Mismatches",  file=self._file, end=self._seperator)
        print("tCount",  file=self._file, end=self._seperator)
        print("tcCount",  file=self._file, end=self._seperator)
        print("tcRate",  file=self._file, end=self._seperator)
        print("ConversionRates",  file=self._file)
        
    
    def write(self, slamSeqRead):
        print(slamSeqRead.name,  file=self._file, end=self._seperator)
        print(slamSeqRead.direction,  file=self._file, end=self._seperator)
        print(slamSeqRead.sequence,  file=self._file, end=self._seperator)
        print(slamSeqRead.tCount,  file=self._file, end=self._seperator)
        print(slamSeqRead.tcCount,  file=self._file, end=self._seperator)
        print(slamSeqRead.tcRate,  file=self._file, end=self._seperator)
        print(slamSeqRead.conversionRates,  file=self._file, end=self._seperator)
        
        for mismatch in slamSeqRead.mismatches:
            print(mismatch,  file=self._file, end=";")
        
        print(file=self._file)
        
    def close(self):
        self._file.close()
        
    
    

class SlamSeqBamIterator:
    
    _readIterator = None
    _refSeq = None
    _snps = None
    _maxReadLength = 0
    _minQual = 0
    _chromosome = None
    _startPosition = 0

    # TODO: merge with fillMismatches
    # pysam code (get_aligend_pairs) should be present only once!        
    def computeRatesForRead(self, read):
        rates = SlamSeqConversionRates()
    
        for pair in read.get_aligned_pairs(matches_only=True):
            alnPosition = self.toAlignmentPos(pair, self._refSeq, read)
            if(alnPosition.readBaseQlty >= self._minQual):            
                rates.incRate(alnPosition.referenceBase, alnPosition.readBase)
                        
        return rates

    def computeRatesForReadNGM(self, read):
        ratesNgm = None
        if(read.has_tag("RA")):
            ratesNgm = map(int, read.get_tag("RA").split(","))
        return ratesNgm
        
    def getTCNgm(self, read):
#         tCount = 0
#         if(read.is_reverse):
#             tCount = read.query_sequence.lower().count("A")
#         else:
#             tCount = read.query_sequence.lower().count("T")
#       TODO: let NGM output T count in reference and return it here
        return int(read.get_tag("TC")), int(read.get_tag("TC")) 
    
    def toAlignmentPos(self, pysamPosition, referenceSequence, read):
        readPos = pysamPosition[0]
        refPos = pysamPosition[1] - int(self._startPosition) + self._maxReadLength + 1
        
        if(refPos < len(referenceSequence)):
            refBase = self._refSeq[refPos]
        else:
            refBase = 'N'
            
        # After retreiving reference base, shift back by read length to retain position within UTR
        refPos = refPos - self._maxReadLength
        
        readBase = read.query_sequence[readPos]
        readQlty = read.query_qualities[readPos]
        
        if(read.is_reverse):
            isSnpPos = self._snps != None and self._snps.isAGSnp(self._chromosome, int(pysamPosition[1]))
            return SlamSeqAlignmentPosition(read.query_length - readPos - 1, refPos, readBase, refBase, readQlty, isSnpPos)
        else:
            isSnpPos = self._snps != None and self._snps.isTCSnp(self._chromosome, int(pysamPosition[1]))
            return SlamSeqAlignmentPosition(readPos, refPos, readBase, refBase, readQlty, isSnpPos)
    
    
    def fillMismatches(self, read):
        mismatchList = []
        tCount = 0
        for pair in read.get_aligned_pairs(matches_only=True):
            alnPosition = self.toAlignmentPos(pair, self._refSeq, read)
            if(alnPosition.isT(read.is_reverse)):
                tCount += 1
            if(alnPosition.isMismatch() and alnPosition.readBaseQlty >= self._minQual):
                mismatchList.append(alnPosition)
        return mismatchList, tCount
          
    
    def getTC(self, mismatches, isReverse):
        tcCount = 0
        
        for mismatch in mismatches:
            if(mismatch.isTCMismatch(isReverse)):
                tcCount += 1
        
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
            if(x != y):
                return False
    
        return True

 
    def next(self):
        read = self._readIterator.next()

        # Create SlamSeqRead
        slamSeqRead = SlamSeqRead()
        slamSeqRead.name = read.query_name
        slamSeqRead.sequence = read.query_sequence
        if(read.is_reverse):
            slamSeqRead.direction = ReadDirection.Reverse
        else:
            slamSeqRead.direction = ReadDirection.Forward
        slamSeqRead.mismatches, slamSeqRead.tCount = self.fillMismatches(read)
        slamSeqRead.tcCount = self.getTC(slamSeqRead.mismatches, read.is_reverse) 
        slamSeqRead.conversionRates = self.computeRatesForRead(read)
        slamSeqRead.tcRate = 0.0
        if(slamSeqRead.tCount > 0):
            slamSeqRead.tcRate = slamSeqRead.tcCount * 100.0 / slamSeqRead.tCount     
        
        
        
        # Get TC count and rates from NGM bam file
#        ngmTC, ngmTCount = self.getTCNgm(read)
#       ngmRates = self.computeRatesForReadNGM(read)
#        ngmRate = 0.0
#        if(ngmTCount > 0):
#            ngmRate = ngmTC * 100.0 / ngmTCount
        
        # Check if results from pysam and NGM are the same
        # TODO: remove at some point
#        if(not self.compareLists(slamSeqRead.conversionRates, ngmRates) or slamSeqRead.tcCount != ngmTC):# or ngmRate != slamSeqRead.tcRate):
#            print("Difference found:")
#            print(read)
#            print(ngmRates)
#            print(slamSeqRead.conversionRates)
#            print("TC (ngm): " + str(ngmTC))
#            print("TC (pys): " + str(slamSeqRead.tcCount))
#            print("TC rate (ngm): " + str(ngmRate))
#            print("TC rate (pys): " + str(slamSeqRead.tcRate))
            #sys.stdin.read(1)
#            raise RuntimeError("Difference found between NGM and Py.")
        
        return slamSeqRead
                
class SlamSeqBamFile:
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
        
        region = chromosome + ":" + str(start) + "-" + str(stop)
        
        if(self.isInReferenceFile(chromosome)):
            refSeq = self._referenceFile.fetch(region=refRegion)
            return SlamSeqBamIterator(self._bamFile.fetch(region=region), refSeq, chromosome, start, maxReadLength, self._snps)
        else:
            return iter([])
    
    def readsInChromosome(self, chromosome):
        refSeq = self._referenceFile.fetch(region=chromosome)
        return SlamSeqBamIterator(self._bamFile.fetch(region=chromosome), refSeq, chromosome, 1, 0, self._snps)
    
    
    def atoi(self, text):
        return int(text) if text.isdigit() else text
    
    def natural_keys(self, text):
        '''
        alist.sort(key=natural_keys) sorts in human order
        http://nedbatchelder.com/blog/200712/human_sorting.html
        (See Toothy's implementation in the comments)
        '''
        return [ self.atoi(c) for c in re.split('(\d+)', text) ]

    def isInReferenceFile(self, chromosome):
        return chromosome in list(self._referenceFile.references)
    
    def getChromosomes(self):
        refs = list(self._referenceFile.references)
        refs.sort(key=self.natural_keys)
        return refs

        
# snps = SNPtools.SNPDictionary("/project/ngs/philipp/slamdunk-analysis/debug/snps/ngm/26338_mESC-wt_0.5h-4SU_trimmed_fixed_downsample_slamdunk_mapped_filtered_snp.vcf")
# testFile = SlamSeqBamFile("/project/ngs/philipp/slamdunk-analysis/debug/filtered/ngm/26338_mESC-wt_0.5h-4SU_trimmed_fixed_downsample_slamdunk_mapped_filtered.bam", "/project/ngs/philipp/slamseq/ref/GRCm38.fa", 55, snps)
# 
# chromosomes = testFile.getChromosomes()
# 
# for chromosome in chromosomes:
#     print(chromosome)
#     readIterator = testFile.readsInChromosome(chromosome)
# 
#     for read in readIterator:
#         print(read.name, read.direction, read.tcCount, read.conversionRates)

