#!/bin/bash

#-------------------------------------------
# Parameters
# Required parameters:
# Transcript annotation (BED file)
BED=input/sample.bed

# output FASTA prefix
FASTAFILE=output/single.fa

# reference chromosome
REFERENCE=input/reference.fa

# Optional parameters
# Read length
READLEN=75

# Number of reads generated
NREAD=100000

# positional bias file
POSBIAS=input/sampleposbias.txt

# Read error position profile
READERR=input/samplereaderror.txt

# Intermediate files
# File for random expression level assignment
RANDEXPLV=output/explvprofile.txt


#-----------------------------------------------
# Add paths if users don't install the script
export PATH=../src/:$PATH
# Commands to randomly assign weights to each transcript

if [ ! -d "output" ]; then
  mkdir output
fi

CMD0="genexplvprofile.py $BED > $RANDEXPLV"

echo "Commands to randomly assign weights to each transcript:"
echo $CMD0

genexplvprofile.py $BED > $RANDEXPLV

# Commands to simulate reads (output to STDOUT in BED format)  
# If you want single-end reads, don't use the "-p" option.
CMD1="gensimreads.py -e $RANDEXPLV  -n $NREAD -b $POSBIAS -l $READLEN  $BED "

echo "Commands to generate simulated paired-reads in BED format:"
echo $CMD1


# Commands to convert BED file to fasta file
CMD2="getseqfrombed.py -b $READERR -f A -r 0.01 -l $READLEN -  $REFERENCE"

echo "Commands to generate FASTA file from the last command:"
echo $CMD2
echo "Output FASTA prefix: $FASTAFILE"



# Execute two commands simultaneously
# If you want single-end reads, do not use the splitfasta.py. Instead, execute:
$CMD1 | $CMD2 > $FASTAFILE




