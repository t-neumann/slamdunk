#!/usr/bin/env python

# Date located in: -
from __future__ import print_function
import sys, os
import subprocess

from os.path import basename
from utils import replaceExtension, files_exist, checkStep, run, runIndexBam, removeExtension, removeFile, runFlagstat

projectPath = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
ngmPath = os.path.join(projectPath, "bin", "NextGenMap", "bin", "ngm-0.4.13", "ngm")


# def runNGMSEBAM(ngmpath, ref, reads, output, verbose=False, dry=False, threads="1", parameter="", log=""):
#     time = ""
#     if(log != ""):
#         time = time + " -a -o " + log
#         log = " &> " + log
#     bamO = replaceExtension(output, ".bam")
#     tmpOutput = " | samtools view -hSb - > " + output + "_tmp"
#     runNGM = True
#     if(files_exist(bamO)):
#         runNGM = checkStep([ref, reads], [bamO]) 
#     else:
#         runNGM = checkStep([ref, reads], [output])
#         
#     if(runNGM):
#         run(time + " " + ngmpath + " -r " + ref + " -q " + reads + " -t " + threads + " " + parameter + tmpOutput + log, verbose=verbose, dry=dry)
#         run("samtools sort " + tmpOutput + removeExtension(output))
#         os.unlink(tmpOutput)
#         runIndexBam(output, verbose=verbose, dry=dry)
#     else:
#         print("Skipped mapping for " + reads)


def runSam2bam(inFile, outFile, index=True, sort=True, delinFile=False, onlyUnique=False, onlyProperPaired=False, filterMQ=0, L=None, verbose=False, dry=False):
    if(delinFile and files_exist(outFile) and not files_exist(inFile)):
        print("Skipping sam2bam for " + outFile)
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
            run(" ".join(cmd), verbose=verbose, dry=dry)
            
            if(sort):         
                tmp = outFile + "_tmp"
                if(not dry):
                    os.rename(outFile, tmp)                      
                run(" ".join(["samtools", "sort", tmp, replaceExtension(outFile, "")]), verbose=verbose, dry=dry)
                if(success):
                    removeFile(tmp)
            if(success and delinFile):
                if(not dry):
                    removeFile(inFile)
        else:
            print("Skipping sam2bam. " + outFile + " already exists.");
        
    if(index):
        runIndexBam(outFile, verbose=verbose, dry=dry)


def Map(inputBAM, inputReference, outputSAM, logFile="", threads=1, parameter="--no-progress --slam-seq 2 -e -i 0.95" , outputSuffix="_ngm_slamdunk", trim5p=0, printOnly=False, verbose=True, force=True):
    time = ""
    log = ""
    if(logFile != ""):
        time = time + " -a -o " + logFile
        log = " &> " + logFile
        
    if(trim5p > 0):
        parameter = parameter + " -5 " + str(trim5p)
                    
    if(checkStep([inputReference, inputBAM], [outputSAM], force)):
        run(time + " " + ngmPath + " -r " + inputReference + " -q " + inputBAM + " -t " + str(threads) + " " + parameter + " -o " + outputSAM + log, verbose=verbose, dry=printOnly)
    else:
        print("Skipped mapping for " + inputBAM)
        

def sort(inputSAM, outputBAM, keepSam=True, dry=False, verbose=True):    

    runSam2bam(inputSAM, outputBAM, True, True, not keepSam, dry=dry, verbose=verbose)
    runFlagstat(outputBAM, dry=dry, verbose=verbose)

