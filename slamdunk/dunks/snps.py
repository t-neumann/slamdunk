#!/usr/bin/env python

# Date located in: -
from __future__ import print_function
import os
import subprocess
import csv
from slamdunk.utils.misc import checkStep, getBinary

projectPath = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
varScanPath = os.path.join(projectPath, "bin", "VarScan.v2.4.1.jar")

def SNPs(inputBAM, outputSNP, referenceFile, minVarFreq, minCov, log, printOnly=False, verbose=True, force=False):
    if(checkStep([inputBAM, referenceFile], [outputSNP], force)):
        fileSNP = open(outputSNP, 'w')
        
        mpileupCmd = getBinary("samtools") + " mpileup -f" + referenceFile + " " + inputBAM
        if(verbose):
            print(mpileupCmd, file=log)
        if(not printOnly):
            mpileup = subprocess.Popen(mpileupCmd, shell=True, stdout=subprocess.PIPE, stderr=log)
            
        varscanCmd = "java -jar " + getBinary("VarScan.v2.4.1.jar") + " mpileup2snp  --strand-filter 0 --output-vcf --min-var-freq " + str(minVarFreq) + " --min-coverage " + str(minCov) + " --variants 1"
        if(verbose):
            print(varscanCmd, file=log)
        if(not printOnly):
            varscan = subprocess.Popen(varscanCmd, shell=True, stdin=mpileup.stdout, stdout=fileSNP, stderr=log)
            varscan.wait()
        
        fileSNP.close()
    else:
        print("Skipping SNP calling", file=log)    
        
def countSNPsInFile(inputFile):
    snpCount = 0
    tcSnpCount = 0
    with open(inputFile, "r") as snpFile:
            snpReader = csv.reader(snpFile, delimiter='\t')
            for row in snpReader:
                if((row[2].upper() == "T" and row[3].upper() == "C") or (row[2].upper() == "A" and row[3].upper() == "G")):
                    tcSnpCount = tcSnpCount + 1
                snpCount = snpCount + 1
    return snpCount, tcSnpCount