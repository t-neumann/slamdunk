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

from slamdunk.utils.misc import checkStep  # @UnresolvedImport
from slamdunk.slamseq.SlamSeqFile import SlamSeqBamFile, SlamSeqWriter  # @UnresolvedImport
from slamdunk.utils import SNPtools  # @UnresolvedImport

def dumpReadInfo(referenceFile, bam, minQual, outputCSV, snpsFile, log, printOnly=False, verbose=True, force=False):
    
    if(not checkStep([bam, referenceFile], [outputCSV], force)):
        print("Skipped computing T->C per reads position for file " + bam, file=log)
    else:
                
        snps = SNPtools.SNPDictionary(snpsFile)
        snps.read()
    
        outputFile = SlamSeqWriter(outputCSV)
        
        #Go through one chr after the other
        testFile = SlamSeqBamFile(bam, referenceFile, snps)
        
        chromosomes = testFile.getChromosomes()
        
        for chromosome in chromosomes:
            readIterator = testFile.readsInChromosome(chromosome)
            for read in readIterator:
                outputFile.write(read)

        
        outputFile.close()