#!/usr/bin/env python

#########################################################################
# Main routine for the SLAMdunk analyzer
#########################################################################
# Imports
#########################################################################

from __future__ import print_function
import sys, os, subprocess

from argparse import ArgumentParser, RawDescriptionHelpFormatter
    
from os.path import basename

########################################################################
# Global variables
########################################################################

logFile = "slamdunk.log"

# Open log
log = open(logFile,'w')

########################################################################
# Routine definitions
########################################################################

def runMap() :
    print("slamdunk map",end="")
    for bam in args.bam :
        print(" " + bam,end="")
    print()
        
def runFilter() :
    print("slamdunk filter",end="")
    for bam in args.bam :
        print(" " + bam,end="")
    print()

def runSnp() :
    
    fasta = args.fasta
    minCov = args.cov
    minVarFreq = args.var

    for bam in args.bam :
        
        mpileupCmd = "samtools mpileup -f" + fasta + " " + bam
        mpileup = subprocess.Popen(mpileupCmd, shell=True, stdout = subprocess.PIPE, stderr = log)
        
        varscanCmd = "java -jar ~/bin/VarScan.v2.4.1.jar mpileup2snp  --strand-filter 0 --min-var-freq " + str(minVarFreq) + " --min-coverage " + str(minCov) + " --variants 1"
        varscan = subprocess.Popen(varscanCmd, shell=True, stdin = mpileup.stdout, stdout = sys.stdout,stderr = log)
        varscan.wait()
        
def runDedup() :
    print("slamdunk dedup",end="")
    for bam in args.bam :
        print(" " + bam,end="")
    print()
        
def runCount() :
    print("slamdunk count",end="")
    for bam in args.bam :
        print(" " + bam,end="")
    print()
        
def runStats() :
    print("slamdunk stats",end="")
    for bam in args.bam :
        print(" " + bam,end="")
    print()
        
def runAll() :
    print("slamdunk all",end="")
    for bam in args.bam :
        print(" " + bam,end="")
    print()

########################################################################
# Argument parsing
########################################################################

# Info
usage = "SLAMdunk software for analyzing SLAM-seq data"
version = "1.0"

# Main Parsers
parser = ArgumentParser(description=usage, formatter_class=RawDescriptionHelpFormatter,version=version)

# Initialize Subparsers
subparsers = parser.add_subparsers(help="",dest="command")

# map command

mapparser = subparsers.add_parser('map', help='Map SLAM-seq read data')
mapparser.add_argument('bam', action='store', help='Bam file(s)' ,nargs="+")
#mapparser.add_argument("-i", "--input", type=str, required=True, dest="bam", help="BAM file" )

# filter command

filterparser = subparsers.add_parser('filter', help='Filter SLAM-seq aligned data')
filterparser.add_argument('bam', action='store', help='Bam file(s)' ,nargs="+")

# snp command

snpparser = subparsers.add_parser('snp', help='Call SNPs on SLAM-seq aligned data')
snpparser.add_argument('bam', action='store', help='Bam file(s)' ,nargs="+")
snpparser.add_argument("-f", "--fasta", required=True, dest="fasta", type=str, help="Reference fasta file")
snpparser.add_argument("-c", "--min-coverage", required=False, dest="cov", type=int, help="Minimimum coverage to call variant",default = 10)
snpparser.add_argument("-a", "--var-fraction", required=False, dest="var", type=float, help="Minimimum variant fraction variant",default = 0.8)

# dedup command

dedupparser = subparsers.add_parser('dedup', help='Deduplicate SLAM-seq aligned data')
dedupparser.add_argument('bam', action='store', help='Bam file(s)' ,nargs="+")

# count command

countparser = subparsers.add_parser('count', help='Count T/C conversions in SLAM-seq aligned data')
countparser.add_argument('bam', action='store', help='Bam file(s)' ,nargs="+")

# stats command

countparser = subparsers.add_parser('stats', help='Calculate stats on SLAM-seq datasets')
countparser.add_argument('bam', action='store', help='Bam file(s)' ,nargs="+")

# all command

countparser = subparsers.add_parser('all', help='Run entire SLAMdunk analysis')
countparser.add_argument('bam', action='store', help='Bam file(s)' ,nargs="+")

args = parser.parse_args()

########################################################################
# Routine selection
########################################################################

command = args.command

if (command == "map") :
    runMap()
elif (command == "filter") :
    runFilter()
elif (command == "snp") :
    runSnp()
elif (command == "dedup") :
    runDedup()
elif (command == "count") :
    runCount()
elif (command == "stats") :
    runStats()
elif (command == "all") :
    runAll()
    
#########################################################################
# Cleanup
########################################################################

log.close()
    
sys.exit(0)