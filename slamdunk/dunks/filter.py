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

# Date located in: -
from __future__ import print_function
import pysam, random, os

from slamdunk.version import __version__, __bam_version__  # @UnresolvedImport

from slamdunk.utils.BedReader import bedToIntervallTree  # @UnresolvedImport
from slamdunk.utils.misc import checkStep, run, removeFile, getBinary, pysamIndex, SlamSeqInfo, md5  # @UnresolvedImport

# def Filter_old(inputBAM, outputBAM, log, MQ=2, printOnly=False, verbose=True, force=True):
#     if(printOnly or checkStep([inputBAM], [outputBAM], force)):
#         run(" ".join([ getBinary("samtools"), "view -q", str(MQ), "-b", inputBAM, ">", outputBAM]), log, verbose=verbose, dry=printOnly)
#     else:
#         print("Skipped filtering for " + inputBAM, file=log)
#
#     runIndexBam(outputBAM, log, verbose=verbose, dry=printOnly)
#     runFlagstat(outputBAM, log, verbose=verbose, dry=printOnly)

def bamSort(outputBAM, log, newHeader, verbose):

    tmp = outputBAM + "_tmp"
    if(newHeader != None):
        pyOutputBAM = pysam.AlignmentFile(outputBAM, "rb")
        pyTmp = pysam.AlignmentFile(tmp, "wb", header=newHeader)
        for read in pyOutputBAM:
            pyTmp.write(read)
        pyOutputBAM.close()
        pyTmp.close()
    else:
        os.rename(outputBAM, tmp)

    #run(" ".join(["samtools", "sort", "-@", str(threads) , tmp, replaceExtension(outFile, "")]), log, verbose=verbose, dry=dry)
    run(" ".join(["samtools sort", "-o", outputBAM, tmp]), log, verbose=verbose, dry=False)
    #pysam.sort(tmp, outputBAM)  # @UndefinedVariable
    removeFile(tmp)

def dumpBufferToBam (buffer, multimapList, outbam, infile):
    # Randomly write hit from read
    #read = random.choice(buffer.values()).pop()
    read = list(buffer.values()).pop().pop()

#     printer = read.query_name + "\t" + infile.getrname(read.reference_id) + "\t" + str(read.reference_start) + "\t" + str(read.reference_end) + "\tPRINT\tTrue"
    read.set_tag("RD", multimapList.rstrip(" "), "Z")
    read.is_secondary = False
    read.is_supplementary = False
    outbam.write(read)

#     return printer
#     for key in buffer.keys():
#         for read in buffer[key]:
#             outbam.write(read)

def multimapUTRRetainment (infile, outfile, bed, minIdentity, NM, log):

    mappedReads = 0
    unmappedReads = 0
    filteredReads = 0

    mqFiltered = 0
    idFiltered = 0
    nmFiltered = 0

    utrIntervallTreeDict = bedToIntervallTree(bed)

#     debugLog = os.path.join("multimapdebug.log")
#
#     fo = open(debugLog, "w")

    # Buffers for multimappers
    multimapBuffer = {}
    prevRead = ""
    # If read maps to another than previously recorded UTR -> do not dump reads to file
    dumpBuffer = True
    # This string tracks all multiple alignments
    multimapList = ""
#     logList = []

    for read in infile:
        if(not read.is_secondary and not read.is_supplementary):
            if(read.is_unmapped):
                unmappedReads += 1
            else:
                mappedReads += 1

        # First pass general filters
        if(read.is_unmapped):
            continue
        if(float(read.get_tag("XI")) < minIdentity):
            idFiltered += 1
            continue
        if(NM > -1 and int(read.get_tag("NM")) > NM):
            nmFiltered += 1
            continue
        if (read.mapping_quality == 0) :
            # Previous read was also multimapper
            if (read.query_name != prevRead and prevRead != "") :

                #if (dumpBuffer and (len(multimapBuffer) > 1 or len(multimapBuffer["nonUTR"]) > 0)) :
                if (dumpBuffer and len(multimapBuffer) > 0) :
                    dumpBufferToBam(multimapBuffer, multimapList, outfile, infile)
                    filteredReads += 1

#                     ret = dumpBufferToBam(multimapBuffer, outfile, infile)
#                     print(ret,file = fo)
                    #multimapBuffer = {}
                    #multimapBuffer["nonUTR"] = []

#                 for entry in logList:
#                     print(prevRead + "\t" + entry + "\t" + str(dumpBuffer), file = fo)
#                 logList = []

                dumpBuffer = True
                multimapList = ""
                multimapBuffer = {}

            # Query Intervall tree for given chromosome for UTs
            chr = infile.getrname(read.reference_id)
            start = read.reference_start
            end = read.reference_end

            if (chr in utrIntervallTreeDict) :
                query = utrIntervallTreeDict[chr][start:end]
            else :
                query = set()

            if len(query) > 0:
                # First UTR hit is recorded without checks
                if (len(multimapBuffer) == 0) :
                    for result in query :
                        if (not result.data in multimapBuffer) :
                            multimapBuffer[result.data] = []
                        multimapBuffer[result.data].append(read)
                # Second UTR hit looks at previous UTR hits -> no dump if hit on different UTR
                else :
                    for result in query :
                        if (not result.data in multimapBuffer) :
                            multimapBuffer[result.data] = []
                            multimapBuffer[result.data].append(read)
                            dumpBuffer = False
                        else :
                            multimapBuffer[result.data].append(read)

#             else :
#                 # If no overlap -> nonUTR
#                 multimapBuffer["nonUTR"].append(read)
#                 for result in query :
#                     logList.append(chr + "\t" + str(start) + "\t" + str(end) + "\t" + result.data)
#             else :
#                 logList.append(chr + "\t" + str(start) + "\t" + str(end) + "\t" + "OFF")

            multimapList = multimapList + chr + ":" + str(start) + "-" + str(end) + " "

            prevRead = read.query_name
        else :
            # Dump any multimappers before a unique mapper
            #if (len(multimapBuffer) > 1 or len(multimapBuffer["nonUTR"]) > 0) :
            if (len(multimapBuffer) > 0) :
                if (dumpBuffer) :
                    dumpBufferToBam(multimapBuffer, multimapList, outfile, infile)
                    filteredReads += 1
#                     ret = dumpBufferToBam(multimapBuffer, outfile, infile)
#                     print(ret,file = fo)
                multimapBuffer = {}
#                 for entry in logList:
#                     print(prevRead + "\t" + entry + "\t" + str(dumpBuffer), file = fo)
#                 logList = []
                #multimapBuffer["nonUTR"] = []
                dumpBuffer = True
                multimapList = ""

            # Record all unique mappers
            prevRead = read.query_name
            outfile.write(read)
            filteredReads += 1

    # Dump last portion if it was multimapper
    #if (dumpBuffer and (len(multimapBuffer) > 1 or len(multimapBuffer["nonUTR"]) > 0)) :
    if (dumpBuffer and len(multimapBuffer) > 0) :
        dumpBufferToBam(multimapBuffer, multimapList, outfile, infile)
        filteredReads += 1

    multimapper = mappedReads - filteredReads - idFiltered - nmFiltered

    print("Criterion\tFiltered reads",file=log)
    print("MQ < 0\t0",file=log)
    print("ID < " + str(minIdentity) + "\t" + str(idFiltered),file=log)
    print("NM > " + str(NM) + "\t" + str(nmFiltered),file=log)
    print("MM\t" + str(multimapper),file=log)

#     fo.close()
    return mappedReads, unmappedReads, filteredReads, mqFiltered, idFiltered, nmFiltered, multimapper


def Filter(inputBAM, outputBAM, log, bed, MQ=2, minIdentity=0.8, NM=-1, printOnly=False, verbose=True, force=False):
    if(printOnly or checkStep([inputBAM], [outputBAM], force)):

        mappedReads = 0
        unmappedReads = 0
        filteredReads = 0

        mqFiltered = 0
        idFiltered = 0
        nmFiltered = 0
        multimapper = 0

        infile = pysam.AlignmentFile(inputBAM, "rb")
        outfile = pysam.AlignmentFile(outputBAM, "wb", template=infile)

        # Default filtering without bed
        if (bed == None) :

            print("#No bed-file supplied. Running default filtering on " + inputBAM + ".",file=log)

            for read in infile:

                if(not read.is_secondary and not read.is_supplementary):
                    if(read.is_unmapped):
                        unmappedReads += 1
                    else:
                        mappedReads += 1

                if(read.is_unmapped):
                    continue
                if(read.mapping_quality < MQ):
                    mqFiltered += 1
                    continue
                if(float(read.get_tag("XI")) < minIdentity):
                    idFiltered += 1
                    continue
                if(NM > -1 and int(read.get_tag("NM")) > NM):
                    nmFiltered += 1
                    continue

                if(not read.is_secondary and not read.is_supplementary):
                    filteredReads += 1

                outfile.write(read)

            print("Criterion\tFiltered reads",file=log)
            print("MQ < " + str(MQ) + "\t" + str(mqFiltered),file=log)
            print("ID < " + str(minIdentity) + "\t" + str(idFiltered),file=log)
            print("NM > " + str(NM) + "\t" + str(nmFiltered),file=log)
            print("MM\t0",file=log)
        else :
            # Multimap retention strategy filtering when bed is supplied

            random.seed(1)

            print("#Bed-file supplied. Running multimap retention filtering strategy on " + inputBAM + ".",file=log)

            mappedReads, unmappedReads, filteredReads, mqFiltered, idFiltered, nmFiltered, multimapper = multimapUTRRetainment (infile, outfile, bed, minIdentity, NM, log)
            #mappedReads, unmappedReads, filteredReads = multimapUTRRetainment (infile, outfile, bed, minIdentity, NM, log)

        # Add number of sequenced and number of mapped reads to the read group description
        # Used for creating summary file
        inFileBamHeader = outfile.header
        if('RG' in inFileBamHeader and len(inFileBamHeader['RG']) > 0):
            slamseqInfo = SlamSeqInfo()
            slamseqInfo.SequencedReads = mappedReads + unmappedReads
            slamseqInfo.MappedReads = mappedReads
            slamseqInfo.FilteredReads = filteredReads
            slamseqInfo.MQFilteredReads = mqFiltered
            slamseqInfo.IdFilteredReads = idFiltered
            slamseqInfo.NmFilteredReads = nmFiltered
            slamseqInfo.MultimapperReads = multimapper

            if (bed != None) :
                slamseqInfo.AnnotationName = os.path.basename(bed)
                slamseqInfo.AnnotationMD5 = md5(bed)
            else :
                slamseqInfo.AnnotationName = ""
                slamseqInfo.AnnotationMD5 = ""

            if not isinstance(inFileBamHeader, dict):
                inFileBamHeader = inFileBamHeader.to_dict()
            inFileBamHeader['RG'][0]['DS'] = str(slamseqInfo)
            #inFileBamHeader['RG'][0]['DS'] = "{'sequenced':" + str(mappedReads + unmappedReads) + "," + "'mapped':" + str(mappedReads) + "," + "'filtered':" + str(filteredReads) + "}"

        slamDunkPG = { 'ID': 'slamdunk', 'PN': 'slamdunk filter v' + __version__, 'VN': __bam_version__ }
        if('PG' in inFileBamHeader):
            inFileBamHeader['PG'].append(slamDunkPG)
        else:
            inFileBamHeader['PG'] = [ slamDunkPG ]

        infile.close()
        outfile.close()

        # Sort afterwards
        bamSort(outputBAM, log, inFileBamHeader, verbose)

        pysamIndex(outputBAM)
        #pysamFlagstat(outputBAM)
        #runFlagstat(outputBAM, log, verbose=verbose, dry=printOnly)

    else:
        print("Skipped filtering for " + inputBAM, file=log)
