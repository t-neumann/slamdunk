#!/usr/bin/env python

# Copyright (c) 2015 Tobias Neumann, Philipp Rescheneder.
#
# This file is part of Slamdunk.
#
# Slamdunk is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as
# published by the Free Software Foundation, either version 3 of the
# License, or (at your option) any later version.
#
# Slamdunk is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Affero General Public License for more details.
#
# You should have received a copy of the GNU Affero General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

# Date located in: -
from __future__ import print_function
import sys, os
import pysam
import subprocess
import collections
import csv
import ast
import hashlib

ReadStat = collections.namedtuple('ReadStat' , 'SequencedReads MappedReads DedupReads FilteredReads SNPs AnnotationName AnnotationMD5')
SampleInfo = collections.namedtuple('SampleInfo' , 'ID Name Type Time')

class SlamSeqInfo:

    ID_SequencedRead = "sequenced"
    ID_MappedReads = "mapped"
    ID_FilteredReads = "filtered"
    ID_DedupReads = "dedup"
    ID_MQFilteredReads = "mqfiltered"
    ID_IdFilteredReads = "idfiltered"
    ID_NmFilteredReads = "nmfiltered"
    ID_MultimapperReads = "multimapper"
    ID_SNPs = "snps"
    ID_AnnotationName = "annotation"
    ID_AnnotationMD5 = "annotationmd5"

    def getFromReadStat(self, name, stats):
        if(name in stats):
            return stats[name]
        else:
            return "NA"

    def __init__(self, bam = None):
        if bam is None:
            self.SequencedReads = 0
            self.MappedReads = 0
            self.DedupReads = 0
            self.FilteredReads = 0
            self.MQFilteredReads = 0
            self.IdFilteredReads = 0
            self.NmFilteredReads = 0
            self.MultimapperReads = 0
            self.SNPs = 0
            self.AnnotationName = "NA"
            self.AnnotationMD5 = "NA"
        else:
            DS = ast.literal_eval(getReadGroup(bam)['DS'])

            self.SequencedReads = self.getFromReadStat(self.ID_SequencedRead, DS)
            self.MappedReads = self.getFromReadStat(self.ID_MappedReads, DS)
            self.DedupReads = self.getFromReadStat(self.ID_DedupReads, DS)
            self.FilteredReads = self.getFromReadStat(self.ID_FilteredReads, DS)
            self.MQFilteredReads = self.getFromReadStat(self.ID_MQFilteredReads, DS)
            self.IdFilteredReads = self.getFromReadStat(self.ID_IdFilteredReads, DS)
            self.NmFilteredReads = self.getFromReadStat(self.ID_NmFilteredReads, DS)
            self.MultimapperReads = self.getFromReadStat(self.ID_MultimapperReads, DS)
            self.SNPs = self.getFromReadStat(self.ID_SNPs, DS)
            self.AnnotationName = self.getFromReadStat(self.ID_AnnotationName, DS)
            self.AnnotationMD5 = self.getFromReadStat(self.ID_AnnotationMD5, DS)

    def __repr__(self):
        return "{" + "'" + self.ID_SequencedRead + "':" + str(self.SequencedReads) + "," + "'" + self.ID_MappedReads + "':" + str(self.MappedReads) + "," + "'" + self.ID_FilteredReads + "':" + str(self.FilteredReads) + "," + "'" + self.ID_MQFilteredReads + "':" + str(self.MQFilteredReads) + "," + "'" + self.ID_IdFilteredReads + "':" + str(self.IdFilteredReads) + "," + "'" + self.ID_NmFilteredReads + "':" + str(self.NmFilteredReads) + "," + "'" + self.ID_MultimapperReads + "':" + str(self.MultimapperReads) + "," + "'" + self.ID_DedupReads + "':" + str(self.DedupReads) + "," + "'" + self.ID_SNPs + "':" + str(self.SNPs) + "," + "'" + self.ID_AnnotationName + "':'" + str(self.AnnotationName) + "'," + "'" + self.ID_AnnotationMD5 + "':'" + str(self.AnnotationMD5) +  "'}"

def md5(fname):
    hash_md5 = hashlib.md5()
    with open(fname, "rb") as f:
        for chunk in iter(lambda: f.read(4096), b""):
            hash_md5.update(chunk)
    return hash_md5.hexdigest()

def estimateMaxReadLength(bam):

    readfile = pysam.AlignmentFile(bam, "rb")

    minLength = sys.maxsize
    maxLength = 0

    for read in readfile.head(n = 1000) :
        minLength = min(minLength, read.query_length + read.get_tag("XA"))
        maxLength = max(maxLength, read.query_length + read.get_tag("XA"))

    range = maxLength - minLength

    if (range <= 10) :
        return(maxLength + 10)
    else:
        return(-1)

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

def getchar():
    print("Waiting for input", file=sys.stderr)
    sys.stdin.readline()

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
        raise RuntimeError("One or more input files don't exist: " + str(inFiles))
    inFileDate = os.path.getmtime(inFiles[0])
    for x in inFiles[1:]:
        inFileDate = max(inFileDate, os.path.getmtime(x))

    if len(outFiles) > 0 and files_exist(outFiles):
        outFileDate = os.path.getmtime(outFiles[0])
        for x in outFiles[1:]:
            outFileDate = min(outFileDate, os.path.getmtime(x))
        if outFileDate > inFileDate:
            if(force == True):
                return True
            else:
                return False

    return True

def getBinary(name):

    projectPath = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

    return os.path.join(projectPath, "contrib", name)

def getRNASeqReadSimulator(name):

    projectPath = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

    return os.path.join(projectPath, "contrib", "RNASeqReadSimulator", "src", name)

def getPlotter(name):

    projectPath = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

    return os.path.join(projectPath, "plot", name + ".R")

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
            raise RuntimeError("Error while executing command: \"" + cmd + "\"")

def callR(cmd, log=sys.stderr, verbose=False, dry=False):

    if(verbose or dry):
        print(cmd, file=log)

    if(not dry):

        p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, shell=True)
        lines_iterator = iter(p.stdout.readline, b"")
        for line in lines_iterator:
            print(line, end="", file=log) # yield line
        p.wait();
        if(p.returncode != 0):
            raise RuntimeError("Error while executing command: \"" + cmd + "\"")

def pysamIndex(outputBam):
    pysam.index(outputBam)  # @UndefinedVariable

def countReads(bam):
    bamFile = pysam.AlignmentFile(bam)
    mapped = 0
    unmapped = 0
    for read in bamFile.fetch(until_eof=True):
        if(not read.is_secondary and not read.is_supplementary):
            if(read.is_unmapped):
                unmapped += 1
            else:
                mapped += 1
    bamFile.close()
    return mapped, unmapped

def getReadGroup(bam):
    bamFile = pysam.AlignmentFile(bam)
    header = bamFile.header
    bamFile.close()
    if('RG' in header and len(header['RG']) > 0):
        return header['RG'][0]
    else:
        raise RuntimeError("Could not get mapped/unmapped/filtered read counts from BAM file. RG is missing. Please rerun slamdunk filter.")

def getSampleInfo(bam):
    sampleInfo = getReadGroup(bam)
    sampleInfos = sampleInfo['SM'].split(":")
    return SampleInfo(ID = sampleInfo['ID'], Name = sampleInfos[0], Type = sampleInfos[1], Time = sampleInfos[2])

def readSampleNames(sampleNames, bams):
    samples = None

    if(sampleNames != None and files_exist(sampleNames)):
        samples = {}
        with open(sampleNames, "r") as sampleFile:
            samplesReader = csv.reader(sampleFile, delimiter='\t')
            for row in samplesReader:
                samples[removeExtension(row[0])] = row[1]

    return samples

def getSampleName(fileName, samples):
    if samples == None:
        return removeExtension(fileName)
    else:
        for key in samples:
            if(key in fileName):
                return samples[key]

    return

def matchFile(sample, files):
    fileName = None
    for item in files:
        if(sample in item):
            if(fileName == None):
                fileName = item
            else:
                raise RuntimeError("Found more than one matching file in list.")

    return fileName

def complement(seq):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N' : 'N'}
    bases = list(seq)
    bases = [complement[base] for base in bases]
    return ''.join(bases)

def shell(cmd):
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE, shell=True)
    p.wait()
    if(p.returncode != 0):
        raise RuntimeError("Error while executing command: " + cmd)
    else:
        return p.communicate()[0]

def shellerr(cmd, raiseError=True):
    p = subprocess.Popen(cmd, stderr=subprocess.PIPE, shell=True)
    p.wait()
    if(p.returncode != 0 and raiseError == True):
        raise RuntimeError("Error while executing command: " + cmd)
    else:
        return p.communicate()[1]
