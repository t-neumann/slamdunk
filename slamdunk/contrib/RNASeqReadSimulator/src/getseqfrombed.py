#!/usr/bin/env python
"""
This script is used to extract sequences from bed file.

USAGE 
  getseqfrombed.py {OPTIONS} <.bed file|-> <reference .fasta file> 

OPTIONS

  -b/--seqerror [error file]\tSpecify the positional error profile to be used. The file should include at least 100 lines, each containing a positive number. The number at line x is the weight that an error is occured at x% position of the read. If no positional error file specified, uniform weight is assumed.

  -r/--errorrate [error rate]\tSpecify the overall error rate, a real positive number. The number of errors of each read will follow a Poisson distribution with its mean value specified by --errorrate.  Default 0 (no errors).

  -l/--readlen [read length]\tSpecify the read length. Default 75. 

  -f/--fill [seq]\tFill at the end of each read by the sequence seq, if the read is shorter than the read length. Default A (to simulate poly-A tails in RNA-Seq reads).

NOTE

  1. The input .bed file is best to sort according to chromosome names. Use - to input from STDIN.
  2. Biopython and numpy package are required.
  
  3. When applying models, we assume that all sequences are in the same length. The length information is given by the -l parameter. If the sequence length is greater than read length, nucleotides outside the read length will not be simulated for error.

HISTORY

	02/01/2013:
	  Fix a bug with no read errors generated.
	  Fix a bug with error profiles in the minus strand.
	08/25/2011:
    	  Rename makebedseq.py to getseqfrombed.py.
    	  Print results to stdout.
"""

import sys;
import pydoc;
import os;
import random;
import bisect;
import math;
import numpy;
from Bio import SeqIO;
from Bio.SeqRecord import SeqRecord;

# import argparse;
# parser=argparse.ArgumentParser('Extract sequences from bed file');
# parser.add_argument('-b','--seqerror',help='Specify the positional error profile to be used. The file should include at least 100 lines, each containing a positive number. The number at line x is the weight that an error is occured at x% position of the read. If no positional error file specified, uniform weight is assumed.');
# parser.add_argument('-r','--errorrate',type=float,default=0.0,help='Specify the overall error rate, a number between 0 and 1. Default 0 (no errors).');
# parser.add_argument('-l','--readlen',type=int,default=75,help='Specify the read length. Default 75.');
# parser.add_argument('-f','--fill',default='A',help='Fill at the end of each read by the sequence seq, if the read is shorter than the read length. Default A (to simulate poly-A tails in RNA-Seq reads).');

if len(sys.argv)<2:
  print>>sys.stderr, (pydoc.render_doc(sys.modules[__name__]));
  sys.exit();

# analyzing parameters
posweight=[];
errrate=0.00;
readlength=75;
forcelength=False;
filledseq='A';

for i in range(len(sys.argv)):
  if i<len(sys.argv)-1:
    if sys.argv[i]=='-b' or sys.argv[i]=='--seqerror':
      bline=0;
      tbweight=0;
      print>>sys.stderr, ('Using pos bias file'+sys.argv[i+1]);
      for lines in open(sys.argv[i+1]):
        bline=bline+1;
        if bline>100:
          break;
        tbweight=float(lines.strip());
        posweight.append(tbweight);
      if len(posweight)!=100:
        print>>sys.stderr, ('Error: the bias file should include at least 100 lines.');
        sys.exit();
    if sys.argv[i]=='-r' or sys.argv[i]=='--errorrate':
      errrate=float(sys.argv[i+1]);
      if errrate<0: # or errrate>1:
        print>>sys.stderr, ('Error: the error rate should be between 0-1.');
        sys.exit();
      print>>sys.stderr,('Error rate: '+str(errrate));
    if sys.argv[i]=='-l' or sys.argv[i]=='--readlen':
      readlength=int(sys.argv[i+1]);
      print>>sys.stderr,('Read length:'+str(readlength));
    if sys.argv[i]=='-f' or sys.argv[i]=='--fill':
      forcelength=True;
      filledseq=sys.argv[i+1];
      print>>sys.stderr,('Force same read length with filled :'+(filledseq));
    


# construct weight probability for read length, if possible
rlenweight=[];
if len(posweight)!=0:
  kweight=0;
  for i in xrange(readlength):
    nfrac=i*100.0/readlength;
    lower=int(math.floor(nfrac));
    higher=int(math.ceil(nfrac));
    if higher==lower: higher=lower+1;
    #print('higher:'+str(higher)+',lower:'+str(lower));
    if higher<100:
      val=posweight[lower]*(nfrac-lower)+posweight[higher]*(higher-nfrac);
    else:
      val=posweight[99];
    kweight+=val;
    rlenweight.append(kweight);

bedfile=sys.argv[-2];
reffile=sys.argv[-1];
#ofastafile=sys.argv[-1];

# build reference
seqref=SeqIO.index(reffile,'fasta');
refkeys = list(seqref.keys())

# read bed file, and ready for writing
if bedfile!="-":
  fid=open(bedfile);
else:
  fid=sys.stdin;
#ofid=open(ofastafile,'w');
ofid=sys.stdout

nlines=0;

prevchr='';
previndex='';

for lines in fid:
  # update line counter
  nlines=nlines+1;
  if nlines %10000==1:
    print>>sys.stderr,('Processing '+str(nlines)+' lines...');
  # parse lines
  bedfield=lines.strip().split('\t');
  if len(bedfield)!=12:
    print>>sys.stderr,('Error: incorrect number of fields at line %d (should be 12, observed %d)' % (nlines, len(bedfield)) );
    continue;
  # clustering
  fieldrange=[int(bedfield[1]),int(bedfield[2])];
  # parse all exons
  exonlen=[int(x) for x in bedfield[10][:-1].split(',')];
  exonstart=[int(x)+fieldrange[0] for x in bedfield[11][:-1].split(',')];
  if not bedfield[0] in refkeys:
    print>>sys.stderr,('Warning: '+bedfield[0]+ ' not in the reference. Ignore...' );
    continue;
  if bedfield[0]!=prevchr:
    print>>sys.stderr, ('Switching to %s ...' % bedfield[0]);
    prevchr=bedfield[0];
    previndex=seqref[bedfield[0]];
  # extract sequences
  thisseq=SeqRecord('');
  for i in range(len(exonlen)):
    thisseq+=previndex[exonstart[i]:(exonstart[i]+exonlen[i])];
  if forcelength:
    if sum(exonlen)<readlength:
      thisseq+=filledseq*(readlength-sum(exonlen));
  thisseq.id=bedfield[3];
  thisseq.description='';
  # mutation
  nmut=numpy.random.poisson(errrate);
  if nmut>0:
    newseq=thisseq.seq;
    for n in xrange(nmut):
      if len(posweight)==0:
        # uniform distrib
        modifyposition=random.choice(xrange(len(newseq)));
      else:
        rchosen=random.random()*kweight;
        modifyposition=bisect.bisect_right(posweight,rchosen);
      # mutate the position
      if len(newseq)>modifyposition:
        topos=random.choice('ATGC');
        while topos==newseq[modifyposition]:
          topos=random.choice('ATGC');
        print >>sys.stderr,('MUTATION at position '+str(modifyposition)+','+newseq[modifyposition]+'->'+topos);
        # print >>sys.stderr,('SEQ:'+newseq);
        newseq=newseq[:modifyposition]+topos+newseq[(modifyposition+1):];
        # print >>sys.stderr,('SEQ:'+newseq);
    #print>>sys.stderr,('NMUTATION:'+str(nmut));
    #print>>sys.stderr, (str(thisseq.seq));
    #print>>sys.stderr,(newseq);
    thisseq.seq=newseq;
  # reverse-complement the sequence if it is on the negative strand
  if bedfield[5]=='-':
    #print >>sys.stderr,('SEQ:'+thisseq.seq);
    thisseq.seq=thisseq.seq.reverse_complement();
    #print >>sys.stderr,('RVCSEQ:'+thisseq.seq);
  # write to record
  try:
    SeqIO.write(thisseq,ofid,'fasta');
  except ValueError:
    print>>sys.stderr, ('Skip at line '+str(nlines)+', sequence object:');
    print>>sys.stderr, (thisseq);
  
        

# ofid.close();
if bedfile!="-":
  fid.close();
    



