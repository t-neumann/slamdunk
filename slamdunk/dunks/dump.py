#!/usr/bin/env python

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