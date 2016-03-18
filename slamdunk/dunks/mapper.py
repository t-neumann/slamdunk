#!/usr/bin/env python

# Date located in: -
from __future__ import print_function
import os

from utils.misc import replaceExtension, files_exist, checkStep, run, runIndexBam, removeFile, runFlagstat

projectPath = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
ngmPath = os.path.join(projectPath, "bin", "NextGenMap", "bin", "ngm-0.4.13", "ngm")


def runSam2bam(inFile, outFile, log, index=True, sort=True, delinFile=False, onlyUnique=False, onlyProperPaired=False, filterMQ=0, L=None, verbose=False, dry=False):
    if(delinFile and files_exist(outFile) and not files_exist(inFile)):
        print("Skipping sam2bam for " + outFile, file=log)
    else:
        if(onlyUnique and filterMQ == 0):
            filterMQ = 1;
            
        success = True    
        if(dry or checkStep([inFile], [outFile])):        
            cmd = ["samtools", "view", "-Sb", "-o", outFile, inFile]
            if filterMQ > 0:
                cmd+=["-q", str(filterMQ)]
            if onlyProperPaired:
                cmd+=["-f", "2"]
            if not L is None:
                cmd+=["-L", L]
            run(" ".join(cmd), log, verbose=verbose, dry=dry)
            
            if(sort):         
                tmp = outFile + "_tmp"
                if(not dry):
                    os.rename(outFile, tmp)                      
                run(" ".join(["samtools", "sort", tmp, replaceExtension(outFile, "")]), log, verbose=verbose, dry=dry)
                if(success):
                    removeFile(tmp)
            if(success and delinFile):
                if(not dry):
                    removeFile(inFile)
        else:
            print("Skipping sam2bam. " + outFile + " already exists.", file=log);
        
    if(index):
        runIndexBam(outFile, log, verbose=verbose, dry=dry)


def Map(inputBAM, inputReference, outputSAM, log, threads=1, parameter="--no-progress --slam-seq 2 -e" , outputSuffix="_ngm_slamdunk", trim5p=0, printOnly=False, verbose=True, force=True):    
    if(trim5p > 0):
        parameter = parameter + " -5 " + str(trim5p)
                    
    if(checkStep([inputReference, inputBAM], [outputSAM], force)):
        run(ngmPath + " -r " + inputReference + " -q " + inputBAM + " -t " + str(threads) + " " + parameter + " -o " + outputSAM, log, verbose=verbose, dry=printOnly)
    else:
        print("Skipped mapping for " + inputBAM, file=log)
        

def sort(inputSAM, outputBAM, log, keepSam=True, dry=False, verbose=True):    

    runSam2bam(inputSAM, outputBAM, log, True, True, not keepSam, dry=dry, verbose=verbose)
    runFlagstat(outputBAM, log, dry=dry, verbose=verbose)

