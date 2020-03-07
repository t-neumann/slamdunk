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
import pysam
import re, sys

class ReadDirection:
    Forward = 1
    Reverse = 2

class SlamSeqConversionRates:

    _baseNumber = 5
    _toBase = [ 'A', 'C', 'G', 'T', 'N' ]

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

    def decRate(self, refBase, readBase):
        self._data[self._baseNumber * self.encodeBase(refBase) + self.encodeBase(readBase)] -= 1

    def getRate(self, refBase, readBase):
        return self._data[self._baseNumber * self.encodeBase(refBase) + self.encodeBase(readBase)]

    def setRate(self, refBase, readBase, count):
        self._data[self._baseNumber * self.encodeBase(refBase) + self.encodeBase(readBase)] = count


    def getData(self):
        return self._data

    def setData(self, data):
        self._data = data

class SlamSeqInterval:

    Header = "\t".join(["Chromosome", "Start", "End", "Name", "Length", "Strand", "ConversionRate", "ReadsCPM", "Tcontent", "CoverageOnTs", "ConversionsOnTs", "ReadCount", "TcReadCount", "multimapCount", "ConversionRateLower", "ConversionRateUpper"])

    def __init__(self, chromosome, start, stop, strand, name, Tcontent, readsCPM, coverageOnTs, conversionsOnTs, conversionRate, readCount, tcReadCount, multimapCount, conversionRateLower = -1.0, conversionRateUpper = -1.0):
        self._chromosome = chromosome
        self._start = start
        self._stop = stop
        self._strand = strand
        self._name = name
        self._Tcontent = Tcontent
        self._readsCPM = readsCPM
        self._coverageOnTs = coverageOnTs
        self._conversionsOnTs = conversionsOnTs
        self._conversionRate = conversionRate
        self._conversionRateLower = conversionRateLower
        self._conversionRateUpper = conversionRateUpper
        self._readCount = readCount
        self._tcReadCount = tcReadCount
        self._multimapCount = multimapCount

    def __repr__(self):
        return (self._chromosome + "\t" + str(self._start) + "\t" + str(self._stop) + "\t" + self._name + "\t" + str(self._stop - self._start) + "\t" + self._strand + "\t" + str(self._conversionRate) + "\t" + str(self._readsCPM) + "\t" + str(self._Tcontent) + "\t" + str(self._coverageOnTs) + "\t" + str(self._conversionsOnTs) + "\t" + str(self._readCount) + "\t" + str(self._tcReadCount) + "\t" + str(self._multimapCount) + "\t" + str(self._conversionRateLower) + "\t" + str(self._conversionRateUpper))

class SlamSeqAlignmentPosition:

    def __init__(self, readPosition, referencePosition, readBase, referenceBase, readBaseQlty, isSnpPos, referenceBase5PrimeContext, referenceBase3PrimeContext):
        # Base in read
        self.readBase = readBase
        # Base quality
        self.readBaseQlty = readBaseQlty
        # Position on the read.
        # ReadPos 0 on a reverse read is the last position of the read reported
        # in the BAM file
        self.readPosition = readPosition
        # Base in reference
        self.referenceBase = referenceBase
        # Position in the reference sequence were the mismatch occurred
        self.referencePosition = referencePosition
        # If the mismatch overlaps with a SNP position (if vcf file is specified)
        self.isSnpPosition = isSnpPos
        # Context of reference base
        self.referenceBase5PrimeContext = referenceBase5PrimeContext
        self.referenceBase3PrimeContext = referenceBase3PrimeContext

    def isMismatch(self):
        return self.readBase != self.referenceBase

    def __repr__(self):
        return str(self.referencePosition) + "," + self.referenceBase + "," + self.referenceBase5PrimeContext + "," + self.referenceBase3PrimeContext + "," + str(self.readPosition) + "," + self.readBase + "," + str(self.readBaseQlty) + "," +  str(self.isSnpPosition)

    def __eq__(self, other):
        return (isinstance(other, SlamSeqAlignmentPosition) and self.readBase == other.readBase and self.referenceBase == other.referenceBase
            and self.readBaseQlty == other.readBaseQlty and self.readPosition == other.readPosition and self.referencePosition == other.referencePosition
            and self.isSnpPosition == other.isSnpPosition)

    def __ne__(self, other):
        return not self.__eq__(other)

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

    def __init__(self):
        # Name of the parsed read
        self.name = None
        # Number of Ts that were converted
        # to a C on forward reads and A to G on reverse reads
        self.tcCount = None
        # Number Ts in the reference
        self.tCount = None
        # Percentage of converted Ts/As
        self.tcRate = None
        # Number of all possible conversions for one reads
        self.conversionRates = None
        # Direction of the reads
        self.direction = None
        # Read sequence
        self.sequence = None
        # List of mismatches in alignment
        self.mismatches = []
        # Start position in the reference
        self.startRefPos = None
        # End position in the reference
        self.endRefPos = None
        # Chr
        self.chromosome = None
        # Multiple TC-conversion flag
        self.isTcRead = None
        # Read is Multimapper
        self.isMultimapper = None

    def getTcount(self):
        if(self.direction == ReadDirection.Reverse):
            return self.sequence.count("a") + self.sequence.count("A")
        else:
            return self.sequence.count("t") + self.sequence.count("T")

    def __repr__(self):
        return "\t".join([self.name, str(self.direction), self.sequence, str(self.tcCount), str(self.tCount), str(self.tcRate), self.conversionRates.__repr__(), str(self.startRefPos), str(self.endRefPos), self.mismatches.__repr__(), str(self.isTcRead), str(self.isMultimapper)])

class SlamSeqWriter:

    _seperator = '\t'

    def __init__(self, fileName):
        self._file = open(fileName, "w")
        self._printHeader()

    def _printHeader(self):
        print("Name", file=self._file, end=self._seperator)
        print("Direction",  file=self._file, end=self._seperator)
        print("Sequence",  file=self._file, end=self._seperator)
        print("Mismatches",  file=self._file, end=self._seperator)
        print("tcCount",  file=self._file, end=self._seperator)
        print("ConversionRates",  file=self._file)


    def write(self, slamSeqRead):
        print(slamSeqRead.name,  file=self._file, end=self._seperator)
        print(slamSeqRead.direction,  file=self._file, end=self._seperator)
        print(slamSeqRead.sequence,  file=self._file, end=self._seperator)
        print(slamSeqRead.tcCount,  file=self._file, end=self._seperator)
        print(slamSeqRead.conversionRates,  file=self._file, end=self._seperator)

        for mismatch in slamSeqRead.mismatches:
            print(mismatch,  file=self._file, end=";")

        print(file=self._file)

    def close(self):
        self._file.close()




class SlamSeqBamIterator:

    def getRefSeq(self):
        return self._refSeq[self._maxReadLength:-self._maxReadLength]

    # @deprecated: Pysam calculation of conversion rates
    def computeRatesForRead(self, read, mismatches):
        rates = SlamSeqConversionRates()

        for base in SlamSeqConversionRates._toBase:
            baseCount = read.query_alignment_sequence.count(base)
#             print(base, baseCount)
            rates.setRate(base, base, baseCount)

        for mismatch in mismatches:

            if not mismatch.isSnpPosition and mismatch.readBaseQlty >= self._minQual:
                rates.incRate(mismatch.referenceBase, mismatch.readBase)
                rates.decRate(mismatch.readBase, mismatch.readBase)
            else:
                rates.incRate(mismatch.referenceBase, mismatch.referenceBase)
                rates.decRate(mismatch.readBase, mismatch.readBase)


        return rates

    #               Read
    #             A     C     G     T     N
    #      A      0     1     2     3     4
    # R    C      5     6     7     8     9
    # e    G     10    11    12    13    14
    # f    T     15    16    17    18    19
    #      N     20    21    22    23    24

    def MPTagToConversion(self, MPTag):
        if (MPTag == "0") :
            return "A","A"
        if (MPTag == "1") :
            return "A","C"
        if (MPTag == "2") :
            return "A","G"
        if (MPTag == "3") :
            return "A","T"
        if (MPTag == "4") :
            return "A","N"
        if (MPTag == "5") :
            return "C","A"
        if (MPTag == "6") :
            return "C","C"
        if (MPTag == "7") :
            return "C","G"
        if (MPTag == "8") :
            return "C","T"
        if (MPTag == "9") :
            return "C","N"
        if (MPTag == "10") :
            return "G","A"
        if (MPTag == "11") :
            return "G","C"
        if (MPTag == "12") :
            return "G","G"
        if (MPTag == "13") :
            return "G","T"
        if (MPTag == "14") :
            return "G","N"
        if (MPTag == "15") :
            return "T","A"
        if (MPTag == "16") :
            return "T","C"
        if (MPTag == "17") :
            return "T","G"
        if (MPTag == "18") :
            return "T","T"
        if (MPTag == "19") :
            return "T","N"
        if (MPTag == "20") :
            return "N","A"
        if (MPTag == "21") :
            return "N","C"
        if (MPTag == "22") :
            return "N","G"
        if (MPTag == "23") :
            return "N","T"
        if (MPTag == "24") :
            return "N","N"

        return None

    def fillMismatchesNGM(self, read):

        tcCount = 0

        mismatchList = []
        if (read.has_tag("MP")) :

            ngmMismatches = read.get_tag("MP").split(",")
            for mismatch in ngmMismatches:
                conversion, readPos, refPos = mismatch.split(":")
                refBase, readBase = self.MPTagToConversion(conversion)
                readPos = int(readPos) - 1

                if readPos >= len(read.query_qualities):
                    readQlty = 0
                else :
                    readQlty = read.query_qualities[readPos]
                if readQlty >= self._minQual:
                    refPos = int(refPos) - 1

                    if(read.is_reverse):
                        isSnpPos = self._snps != None and self._snps.isAGSnp(self._chromosome, read.reference_start + int(refPos))
                        readPos = read.query_length - readPos - 1
                    else :
                        isSnpPos = self._snps != None and self._snps.isTCSnp(self._chromosome, read.reference_start + int(refPos))

                    refPos = read.reference_start - self._startPosition + refPos

                    alnPos = SlamSeqAlignmentPosition(readPos, refPos, readBase, refBase, readQlty, isSnpPos, "N","N")

                    if (alnPos.isTCMismatch(read.is_reverse)) :
                        tcCount += 1

                    mismatchList.append(alnPos)

        return mismatchList, tcCount

    def __init__(self, readIterator, refSeq, chromosome, startPosition, strand, maxReadLength, snps, minQual, conversionThreshold = 1):
        self._readIterator = readIterator
        self._refSeq = refSeq
        self._chromosome = chromosome
        self._startPosition = startPosition
        self._strand = strand
        self._maxReadLength = maxReadLength
        self._snps = snps
        self._minQual = minQual
        self._conversionThreshold = conversionThreshold

    def __iter__(self):
        return self

    def __next__(self):

        read = self._readIterator.__next__()

        # Strand-specific assay - skip all reads from antisense-strand
        while((self._strand == "+" and read.is_reverse) or (self._strand == "-" and not read.is_reverse)) :
            read = self._readIterator.__next__()

        # Create SlamSeqRead
        slamSeqRead = SlamSeqRead()
        slamSeqRead.name = read.query_name
        slamSeqRead.sequence = read.query_sequence
        if(read.is_reverse):
            slamSeqRead.direction = ReadDirection.Reverse
        else:
            slamSeqRead.direction = ReadDirection.Forward

        if (read.mapping_quality == 0) :
            slamSeqRead.isMultimapper = True
        else :
            slamSeqRead.isMultimapper = False

        slamSeqRead.mismatches, slamSeqRead.tcCount = self.fillMismatchesNGM(read)
        slamSeqRead.conversionRates = self.computeRatesForRead(read, slamSeqRead.mismatches)
        slamSeqRead.startRefPos = read.reference_start - int(self._startPosition)
        slamSeqRead.endRefPos = read.reference_end - int(self._startPosition)

        # "reference_name" not available in pysam < 0.9.0
        if hasattr(read, 'reference_name'):
            slamSeqRead.chromosome = read.reference_name
        else:
            slamSeqRead.chromosome = None

        if(slamSeqRead.tcCount >= self._conversionThreshold) :
            slamSeqRead.isTcRead = True
        else :
            slamSeqRead.isTcRead = False

        return slamSeqRead

class SlamSeqBamFile:

    def __init__(self, bamFile, referenceFile, snps):
        self._bamFile = pysam.AlignmentFile(bamFile, "rb")

        # Get version from BAM file
        self.bamVersion = "None"
        if('PG' in self._bamFile.header):
            for pg in self._bamFile.header['PG']:
                if pg['ID'] == "slamdunk":
                    self.bamVersion = pg['VN']

        self._referenceFile = pysam.FastaFile(referenceFile)
        self._snps = snps

    def readInRegion(self, chromosome, start, stop, strand, maxReadLength, minQual = 0, conversionThreshold = 1):

        if(self.isInReferenceFile(chromosome)):
            chromosomeLength = self._referenceFile.get_reference_length(chromosome)

            fillupLeft = 0
            leftBorder = int(start) - maxReadLength
            if (leftBorder < 0) :
                fillupLeft = abs(leftBorder)
                leftBorder = 0

            fillupRight = 0
            rightBorder = int(stop) + maxReadLength
            if(rightBorder > chromosomeLength):
                fillupRight =  rightBorder - chromosomeLength
                rightBorder = chromosomeLength

            refSeq = self._referenceFile.fetch(reference=chromosome, start=leftBorder, end=rightBorder).upper()
            # If start or top is less than maxReadLength bp away from chromosome start or end, fill up with Ns
            refSeq = 'N' * fillupLeft + refSeq + 'N' * fillupRight

            if (chromosome in self._bamFile.references) :
                return SlamSeqBamIterator(self._bamFile.fetch(reference=chromosome, start=max(0, start), end=min(chromosomeLength, stop)), refSeq, chromosome, start, strand, maxReadLength, self._snps, minQual, conversionThreshold)
            else:
                return iter([])
        else:
            return iter([])

    def readsInChromosome(self, chromosome, minQual = 0, conversionThreshold = 1):

        if (chromosome in self._bamFile.references) :
            refSeq = self._referenceFile.fetch(reference=chromosome).upper()
            return SlamSeqBamIterator(self._bamFile.fetch(reference=chromosome), refSeq, chromosome, 1, ".", 0, self._snps, minQual, conversionThreshold)
        else :
            return iter([])

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

    def getChromosomeLength(self, chromosome):
        return self._referenceFile.get_reference_length(chromosome)
