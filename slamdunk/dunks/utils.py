#!/usr/bin/env python

# Date located in: -
from __future__ import print_function
import sys, os
import subprocess

#Replaces the file extension of inFile to with <newExtension> and adds a suffix
#Example replaceExtension("reads.fq", ".sam", suffix="_namg") => reads_ngm.sam
def replaceExtension(inFile, newExtension, suffix=""):    
    return os.path.splitext(inFile)[0] + suffix + newExtension

#Removes right-most extension from file name
def removeExtension(inFile):
    name = os.path.splitext(inFile)[0]
    ext = os.path.splitext(inFile)[1]
    if(ext == ".gz"):
        name = os.path.splitext(name)[0]
    return name

def files_exist(files):
    if (type(files) is list) :
        for f in files:
            if not os.path.exists(f):
                return False
    else:
        if not os.path.exists(files):
            return False
    return True

# remove a (list of) file(s) (if it/they exists)
def removeFile(files):
    if (type(files) is list) :
        for f in files:
            if os.path.exists(f):
                os.remove(f)
    else:
        if os.path.exists(files):
            os.remove(files)


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

def run(cmd, log=sys.stderr, verbose=False, dry=False):
    if(verbose or dry):
        print(cmd, file=log)
    
    if(not dry):
        #ret = os.system(cmd)
        p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, shell=True)
        lines_iterator = iter(p.stdout.readline, b"")
        for line in lines_iterator:
            print(line, end="", file=log) # yield line
        p.wait();
        if(p.returncode != 0):
            raise RuntimeError("Error while executing command!")

def runIndexBam(inFileBam, log=sys.stderr, verbose=False, dry=False):
    idxFile = inFileBam + ".bai"
    if(dry or checkStep([inFileBam], [idxFile])):
        run(" ".join(["samtools", "index", inFileBam]), log, verbose=verbose, dry=dry)

def runFlagstat(bam, log=sys.stderr, verbose=False, dry=False):
    flagstat = bam + ".flagstat"
    if(dry or checkStep([bam], [flagstat])):
        run(" ".join([ "samtools", "flagstat", bam, ">", flagstat]), log, verbose=verbose, dry=dry)
    else:
        print("Skipped flagstat for " + bam)


