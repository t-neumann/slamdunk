#!/usr/bin/env python

# Date located in: -
from __future__ import print_function
#import sys, os
#import subprocess

from utils import checkStep, run, runIndexBam, runFlagstat

def Filter(inputBAM, outputBAM, MQ=2, printOnly=False, verbose=True, force=True):
    if(printOnly or checkStep([inputBAM], [outputBAM], force)):
        run(" ".join([ "samtools", "view -q", str(MQ), "-b", inputBAM, ">", outputBAM]), verbose=verbose, dry=printOnly)
    else:
        print("Skipped filtering for " + inputBAM)

    runIndexBam(outputBAM, verbose=verbose, dry=printOnly)
    runFlagstat(outputBAM, verbose=verbose, dry=printOnly)
    