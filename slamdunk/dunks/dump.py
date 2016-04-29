#!/usr/bin/env python

from __future__ import print_function
import os

from utils.misc import checkStep
from slamseq.SlamSeqFile import SlamSeqBamFile, SlamSeqWriter
from utils import SNPtools

projectPath = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
pathComputeOverallRates = os.path.join(projectPath, "plot", "compute_overall_rates.R")
pathConversionPerReadPos = os.path.join(projectPath, "plot", "conversion_per_read_position.R")

            
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