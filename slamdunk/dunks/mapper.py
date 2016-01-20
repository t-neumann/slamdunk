#!/usr/bin/env python

# Date located in: -
from __future__ import print_function
import sys, os
import subprocess

from os.path import basename

projectPath = os.path.dirname(os.path.abspath(__file__))
ngmPath = os.path.join(projectPath, "tools", "NextGenMap", "bin", "ngm-0.4.13") + "ngm"

#Replaces the file extension of inFile to with <newExtension> and adds a suffix
#Example replaceExtension("reads.fq", ".sam", suffix="_namg") => reads_ngm.sam
def replaceExtension(inFile, newExtension, suffix=""):    
    return os.path.splitext(inFile)[0] + suffix + newExtension

#Removes right-most extension from file name
def removeExtension(inFile):        
    return os.path.splitext(inFile)[0]

def files_exist(files):
    if (type(files) is list) :
        for f in files:
            if not os.path.exists(f):
                return False
    else:
        if not os.path.exists(files):
            return False
    return True

def checkStep(inFiles, outFiles, force=False):    
    if not files_exist(inFiles):
        print(inFiles, outFiles)
        raise RuntimeError("One or more input files don't exist.")
    inFileDate = os.path.getmtime(inFiles[0])
    for x in inFiles[1:]:
        inFileDate = max(inFileDate, os.path.getmtime(x))    
    
    if files_exist(outFiles):
        outFileDate = os.path.getmtime(outFiles[0])
        for x in outFiles[1:]:
            outFileDate = min(outFileDate, os.path.getmtime(x))        
        if outFileDate > inFileDate:
            if(force == True):
                return True
            else:
                #if(os.path.getmtime(os.path.realpath(__file__)) > outFileDate):
                #    logging.info("Script was modified. Rerunning computation.")
                #else:
                #logging.critical("outFiles exist and are newer than inFiles. Skipping execution")
                #sys.exit(0)
                return False
    
    return True

def run(cmd, verbose=False, dry=False):
    if(verbose or dry):
        print(cmd, file=sys.stderr)
    
    if(not dry):
        #ret = os.system(cmd)
        p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, shell=True)
        lines_iterator = iter(p.stdout.readline, b"")
        for line in lines_iterator:
            print(line, end="") # yield line
        p.wait();
        if(p.returncode != 0 != 0):
            raise RuntimeError("Error while executing command!")

def runIndexBam(inFileBam, verbose=False, dry=False):
    idxFile = inFileBam + ".bai"
    if(dry or checkStep([inFileBam], [idxFile])):
        run(" ".join(["samtools", "index", inFileBam]), verbose=verbose, dry=dry)



def runNGMSEBAM(ngmpath, ref, reads, output, verbose=False, dry=False, threads="1", parameter="", log=""):
    time = ""
    if(log != ""):
        time = time + " -a -o " + log
        log = " &> " + log
    bamO = replaceExtension(output, ".bam")
    tmpOutput = " | samtools view -hSb - > " + output + "_tmp"
    runNGM = True
    if(files_exist(bamO)):
        runNGM = checkStep([ref, reads], [bamO]) 
    else:
        runNGM = checkStep([ref, reads], [output])
        
    if(runNGM):
        run(time + " " + ngmpath + " -r " + ref + " -q " + reads + " -t " + threads + " " + parameter + tmpOutput + log, verbose=verbose, dry=dry)
        run("samtools sort " + tmpOutput + removeExtension(output))
        os.unlink(tmpOutput)
        runIndexBam(output, verbose=verbose, dry=dry)
    else:
        print("Skipped mapping for " + reads)



def map(inputBAM, inputReference, inputUTRs, outputBAM, logFile="", threads=1, parameter="--no-progress --slam-seq 2 -e -i 0.95" , outputSuffix="_ngm_slamdunk", trim5p=0, d=False, v=True):
    time = ""
    log = ""
    if(logFile != ""):
        time = time + " -a -o " + logFile
        log = " &> " + logFile
        
    if(trim5p > 0):
        parameter = parameter + " -5 " + str(trim5p)
                    
    if(checkStep([inputReference, inputBAM], [outputBAM])):
        run(time + " " + ngmPath + " -r " + inputReference + " -q " + inputBAM + " -t " + threads + " " + parameter + " -o " + outputBAM + log, verbose=v, dry=d)
    else:
        print("Skipped mapping for " + inputBAM)
        
    
    
    #runNGMSE(ngmslam, refFile, cfile, ngmslamO + ".sam", dry=d, verbose=v, threads="16", parameter="--no-progress " + " -5 " + str(trim5p) + " --slam-seq 2 -e -i 0.95 ", log=ngmslamO + ".log")
