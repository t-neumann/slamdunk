# #!/usr/bin/env python
# 
# from __future__ import print_function
# import sys, os, re
# 
# #Import utils from parent folder
# #Extract absolute path of phtools project
# incPath = os.path.abspath(__file__)
# incName = "phtools"
# incI = incPath.find(incName)
# projectBase = incPath[0:(incI + len(incName) + 1)]
# #Add utils folder to sys path
# sys.path.append(projectBase + "utils/")
# from utils import *
# from slamseq import *
# 
# 
# #sys.path.insert(0,"/home/CIBIV/philipp_/.local/lib/python2.7/site-packages/")
# import pysam
# 
# from argparse import ArgumentParser, RawDescriptionHelpFormatter
#     
# from os.path import basename
# 
# usage = "Computes (T -> C) conversion rates from BAM file"
# parser = ArgumentParser(description=usage, formatter_class=RawDescriptionHelpFormatter)
# parser.add_argument("-i", "--input", type=str, required=True, dest="bam", help="BAM file" )
# parser.add_argument("-r", "--reference", type=str, required=True, dest="ref", help="Reference fasta file" )
# parser.add_argument("-q", "--min-base-qual", type=int, default=0, required=False, dest="minQual", help="Min base quality for T -> C conversions" )
# parser.add_argument("-s", "--snps", type=str, required=False, dest="snp", help="SNP file" )
# 
# parser.add_argument("-t", "--tc-bam", type=str, required=True, dest="tcbam", help="BAM file with TC reads" )
# parser.add_argument("-o", "--nontc-bam", type=str, required=True, dest="nontcbam", help="BAM file with non TC reads" )
# 
# args = parser.parse_args()
# 
# baseNumber = 5
# toBase = [ 'A', 'C', 'G', 'T', 'N' ]
# 
# #verbose
# v = True
# #dry
# d = False
#             
# ref = args.ref
# bam = args.bam
# minQual = args.minQual
# snps = args.snp
# 
# tcbam = args.tcbam
# nontcbam = args.nontcbam
# 
# snpDict = readSNPs(snps)
# 
# reffile = pysam.FastaFile(ref)
# samfile = pysam.AlignmentFile(bam, "rb")
# 
# outfileTC = pysam.AlignmentFile(tcbam, "wb", template=samfile)
# outfileNonTC = pysam.AlignmentFile(nontcbam, "wb", template=samfile)
# 
# totalTC = 0
# totalNonTC = 0
# 
# #Init
# refCount = 1
# totalRefCount = len(reffile.references)
# 
# refs = list(reffile.references)
# refs.sort(key=natural_keys)
# #Go through one chr after the other
# for ref in refs:
#     readCount = 0
#     
#     #Read chr into memory
#     refSeq = reffile.fetch(ref)
#     #Get total number of reads mapping to chr
#     totalCount = samfile.count(ref)
#     #Get all reads on chr
#     for read in samfile.fetch(ref):
#         
#         chr = samfile.getrname(read.reference_id)
#         #try:
# 
# 
#         tc, ag = getTCchr(read, refSeq, chr, minQual, snpDict)
# 
#         if(read.is_reverse):
#             tc = ag
#         
#         if(tc > 0):
#             outfileTC.write(read)
#             totalTC += 1
#         else:
#             outfileNonTC.write(read)
#             totalNonTC += 1
#             
#         #except:
#         #    print("")
#         #    print("Error computing rates for read " + read.query_name)
#             #print("Msg: " + sys.exc_info()[0])
#         #    print(read)
#         #Progress info
#         readCount += 1
#         if(readCount % 1000 == 0):
#             sys.stderr.write("\rProcessing %s (%d/%d): %d%%                        " % (ref, refCount, totalRefCount, int(readCount * 100.0 / totalCount)))
#             sys.stderr.flush()
#             
#     refCount += 1
# 
# #Cleanup
# sys.stderr.write("\n")
# samfile.close()
# reffile.close()
# 
# outfileTC.close()
# outfileNonTC.close()
# 
# runIndexBam(tcbam, True)
# runIndexBam(nontcbam, True)
# 
# print("Reads written: " + str(totalTC) + " / " + str(totalNonTC), file=sys.stderr)
