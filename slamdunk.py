#!/usr/bin/env python

#########################################################################
# Main routine for the SLAMdunk analyzer
#########################################################################

from __future__ import print_function
import sys, os

from argparse import ArgumentParser, RawDescriptionHelpFormatter
    
from os.path import basename

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
    print("slamdunk snp",end="")
    for bam in args.bam :
        print(" " + bam,end="")
    print()
        
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
    
sys.exit(0)