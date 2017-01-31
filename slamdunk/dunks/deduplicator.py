#!/usr/bin/env python

from __future__ import print_function
import pysam

from slamdunk.utils.misc import checkStep, pysamIndex  # @UnresolvedImport

def Dedup(inputBAM, outputBAM, tcMutations, log, printOnly=False, verbose = True, force=False):
    
    if(printOnly or checkStep([inputBAM], [outputBAM], force)):
        
        samfile = pysam.AlignmentFile(inputBAM, "rb")
        outfile = pysam.AlignmentFile(outputBAM, "wb", template=samfile)
        
        processedReads = 0
        retainedReads = 0

        prevChr = ""
        prevStart = ""
        
        duplicateBuffer = {}
        
        for read in samfile:
            
            flag = read.cigarstring
            chr = read.reference_id
            start = read.reference_start
            seq = read.query_sequence
            if (read.has_tag("TC")) :
                tcflag = read.get_tag("TC")
            else :
                tcflag = 0
            
            if (tcflag >= tcMutations) :
                
                if (chr != prevChr or start != prevStart) :
                                
                    if (prevChr != "") :
                        for curSeq in duplicateBuffer :
                            for curFlag in duplicateBuffer[curSeq]:
                                outfile.write(duplicateBuffer[curSeq][curFlag])
                                retainedReads += 1
                        duplicateBuffer.clear()
                
                duplicateBuffer[seq] = {}
                duplicateBuffer[seq][flag] = read
                 
                prevChr = chr
                prevStart = start
            
            processedReads += 1
            
        for seq in duplicateBuffer:
            for flag in duplicateBuffer[seq] :
                outfile.write(duplicateBuffer[seq][flag])
                retainedReads += 1
        duplicateBuffer.clear()
                
        print("Retained " + str(retainedReads) + " of " + str(processedReads) + " reads (", file=log, end = "")
        print("{0:.2f}".format(float(retainedReads) / float(processedReads)),file=log,end="")
        print(" compression rate)", file=log)
        
        pysamIndex(outputBAM)
        
    else:
        print("Skipped deduplication for " + inputBAM, file=log)