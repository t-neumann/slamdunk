#!/usr/bin/env python

# Date located in: -
from __future__ import print_function
import os

from utils.misc import replaceExtension, files_exist, checkStep, run, runIndexBam, removeFile, runFlagstat

projectPath = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
ngmPath = os.path.join(projectPath, "bin", "ngm")

def runSam2bam(inFile, outFile, log, index=True, sort=True, delinFile=False, onlyUnique=False, onlyProperPaired=False, filterMQ=0, L=None, threads=1, verbose=False, dry=False):
    if(delinFile and files_exist(outFile) and not files_exist(inFile)):
        print("Skipping sam2bam for " + outFile, file=log)
    else:
        if(onlyUnique and filterMQ == 0):
            filterMQ = 1;
            
        success = True    
        cmd = ["samtools", "view", "-@", str(threads), "-Sb", "-o", outFile, inFile]
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
            run(" ".join(["samtools", "sort", "-@", str(threads), tmp, replaceExtension(outFile, "")]), log, verbose=verbose, dry=dry)
            if(success):
                removeFile(tmp)
        if(success and delinFile):
            if(not dry):
                removeFile(inFile)
        
    if(index):
        runIndexBam(outFile, log, verbose=verbose, dry=dry)


def Map(inputBAM, inputReference, outputSAM, log, quantseqMapping, localMapping, threads=1, parameter="--no-progress --slam-seq 2" , outputSuffix="_ngm_slamdunk", trim5p=0, topn=1, printOnly=False, verbose=True, force=False):
    
    if(quantseqMapping is True) :
        parameter = "--no-progress"
            
    if(trim5p > 0):
        parameter = parameter + " -5 " + str(trim5p)
                    
    if(localMapping is True):
        parameter = parameter + " -l "
    else:
        parameter = parameter + " -e "
        
    if(topn > 1):
        parameter = parameter + " -n " + str(topn) + " --strata "
        
    if(checkStep([inputReference, inputBAM], [outputSAM], force)):
        run(ngmPath + " -r " + inputReference + " -q " + inputBAM + " -t " + str(threads) + " " + parameter + " -o " + outputSAM, log, verbose=verbose, dry=printOnly)
    else:
        print("Skipped mapping for " + inputBAM, file=log)
        

def sort(inputSAM, outputBAM, log, threads=1, keepSam=True, dry=False, verbose=True):    

    if(files_exist(inputSAM) and checkStep([inputSAM], [outputBAM + ".flagstat"])):
        #runSam2bam(inputSAM, outputBAM, log, True, True, not keepSam, threads=threads, dry=dry, verbose=verbose)
        runSam2bam(inputSAM, outputBAM, log, False, False, not keepSam, threads=threads, dry=dry, verbose=verbose)
        runFlagstat(outputBAM, log, dry=dry, verbose=verbose)
    else:
        print("Skipped sorting for " + inputSAM, file=log)

