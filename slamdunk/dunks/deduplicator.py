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
                                for readEntry in duplicateBuffer[curSeq][curFlag]:
                                    if not readEntry.is_duplicate:
                                       retainedReads += 1 
                                    outfile.write(readEntry)
                        duplicateBuffer.clear()
                
                if not seq in duplicateBuffer:
                    duplicateBuffer[seq] = {}
                if not flag in duplicateBuffer[seq]:
                    duplicateBuffer[seq][flag] = list()
                if len(duplicateBuffer[seq][flag]) > 0 :
                    read.is_duplicate = True
                duplicateBuffer[seq][flag].append(read)
                 
                prevChr = chr
                prevStart = start
            
                processedReads += 1
            
        for seq in duplicateBuffer:
            for flag in duplicateBuffer[seq] :
                for readEntry in duplicateBuffer[seq][flag]:
                    if not readEntry.is_duplicate:
                        retainedReads += 1 
                    outfile.write(readEntry)
        duplicateBuffer.clear()
        
        outfile.close()
                
        print("Retained " + str(retainedReads) + " of " + str(processedReads) + " reads (", file=log, end = "")
        print("{0:.2f}".format(float(retainedReads) / float(processedReads)),file=log,end="")
        print(" compression rate)", file=log)
        
        pysamIndex(outputBAM)
        
    else:
        print("Skipped deduplication for " + inputBAM, file=log)