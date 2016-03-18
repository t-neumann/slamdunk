#!/usr/bin/env python

# Date located in: -
from __future__ import print_function
import pysam

from utils.misc import checkStep, run, runIndexBam, runFlagstat


def Filter_old(inputBAM, outputBAM, log, MQ=2, printOnly=False, verbose=True, force=True):
    if(printOnly or checkStep([inputBAM], [outputBAM], force)):
        run(" ".join([ "samtools", "view -q", str(MQ), "-b", inputBAM, ">", outputBAM]), log, verbose=verbose, dry=printOnly)
    else:
        print("Skipped filtering for " + inputBAM, file=log)

    runIndexBam(outputBAM, log, verbose=verbose, dry=printOnly)
    runFlagstat(outputBAM, log, verbose=verbose, dry=printOnly)


def pysamIndex(outputBam):
    pysam.index(outputBam)

#def pysamFlagstat(outputBam):
#    pysam.flagstat(outputBam)

def Filter(inputBAM, outputBAM, log, MQ=2, minIdentity=0.8, NM=-1, printOnly=False, verbose=True, force=True):
    if(printOnly or checkStep([inputBAM], [outputBAM], force)):
    
        
        infile = pysam.AlignmentFile(inputBAM, "rb")
        outfile = pysam.AlignmentFile(outputBAM, "wb", template=infile)
        for read in infile:
            if(read.is_unmapped):
                continue
            if(read.mapping_quality < MQ):
                continue
            if(float(read.get_tag("XI")) < minIdentity):
                continue
            if(NM > -1 and int(read.get_tag("NM")) > NM):
                continue
            outfile.write(read)
        
        infile.close()
        outfile.close()
        
        pysamIndex(outputBAM)
        #pysamFlagstat(outputBAM)
        runFlagstat(outputBAM, log, verbose=verbose, dry=printOnly)
        
    
    else:
        print("Skipped filtering for " + inputBAM, file=log)

    