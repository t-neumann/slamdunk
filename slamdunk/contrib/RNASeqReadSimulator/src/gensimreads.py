#!/usr/bin/env python 
"""
This script generates simulated RNA-Seq reads (in .bed format) from known gene annotations.

USAGE 

  gensimreads.py {OPTIONS} <BED-File|->

PARAMETER

  BED-File\tThe gene annotation file (in BED format). Use '-' for STDIN input

OPTIONS

  -e/--expression [expression level file] \tSpecify the weight of each transcript. Each line in the file should have at least (NFIELD+1)  fields, with field 0 the annotation id, and field NFIELD the weight of this annoation. If this file is not provided, uniform weight is applied. 

  -n/--nreads readcnt \tSpecify the number of reads to be generated. Default 100000.

  -b/--posbias [positional bias file] \tSpecify the positional bias file. The file should include at least 100 lines, each contains only one integer number, showing the preference of the positional bias at this position. If no positional bias file is specified, use uniform distribution bias.

  -l/--readlen [read length] \tSpecify the read length. Default 32.

  -o/--output [output .bed file] \tSpecify the output file. Default STDOUT 

  -f/--field [NFIELD] \tThe field of each line as weight input. Default 7 (beginning from field 0) to compatible to genexplvprofile.py.

  -p/--pairend [PELENMEAN,PELENSTD]\t Generate paired-end reads with specified insert length mean and standard derivation. The default is 200,20.

  --stranded \tThe reads are strand specific.

NOTE 

  	1. The bed file is required to sort according to the chromosome name and position. In Unix systems, use "sort -k 1,1 -k 2,2n in.BED > out.BED" to get a sorted version (out.BED) of the bed file (in.BED).  

  	2. No problem to handle reads spanning multiple exons. 

HISTORY

	04/30/2012
	  Support generating stranded RNA-Seq reads

	02/16/2012
	  Now runs on python 2.7

	02/08/2012 
	  Change default value of NFIELD from 4 to 7 to be compatible with default genexplvprofile values.

	01/29/2012 
	  Add paired-end support.

	01/09/2012 
	  Add -f option.

AUTHOR
	Wei Li (li.david.wei AT gmail.com)
"""

from __future__ import print_function
import sys;
import subprocess;
import pydoc;
import os;
import random;
import bisect;
import math;
from getSegs import *;

import pdb;

# read length
readlen=32;
# number of reads to sample
readcnt=100000;

nfield=7;

if len(sys.argv)<2:
  print(pydoc.render_doc(sys.modules[__name__]));
  sys.exit();

allids={};
allidl=[];
allexp=[];

posweight=[];

#onbedfile=sys.argv[-1]+'.reads.bed';
onbedfile="-";

genpereads=False;
pemean=200;
pestd=20;

stranded=False;

for i in range(len(sys.argv)):
  if i<len(sys.argv)-1:
    if sys.argv[i]=='-e' or sys.argv[i]=='--expression':
      # parse the annoatation file, and sum up the weights
      nline=0;
      totalweight=0;
      print('Reading annoatation file...',file=sys.stderr);
      for lines in open(sys.argv[i+1]):
        nline=nline+1;
        if lines[0]=='#':
          continue;
        fields=lines.strip().split();
        if len(fields)<nfield+1:
          print('Error: the annotation file should include at least '+str(nfield+1)+' fields.',file=sys.stderr);
          sys.exit();
        allids[fields[0]]=0;
        totalweight+=float(fields[nfield]);
        allexp.append(totalweight);
        allidl.append(fields[0]);
      print('Read %d lines of the annoatation' % nline,file=sys.stderr);
      #print('Total weight: %f' % sum(totalweight));
    if sys.argv[i]=='-b' or sys.argv[i]=='--posbias':
      bline=0;
      tbweight=0;
      for lines in open(sys.argv[i+1]):
        bline=bline+1;
        if bline>100:
          break;
        tbweight=float(lines.strip());
        posweight.append(tbweight);
      if len(posweight)!=100:
        print('Error: the bias file should include at least 100 lines.',file=sys.stderr);
        sys.exit();
    if sys.argv[i]=='-n' or sys.argv[i]=='--nreads':
      readcnt=int(sys.argv[i+1]);
      print('Read count:',readcnt,file=sys.stderr);
    if sys.argv[i]=='-l' or sys.argv[i]=='--readlen':
      readlen=int(sys.argv[i+1]);
      print('Read length:',readlen,file=sys.stderr);
    if sys.argv[i]=='-o' or sys.argv[i]=='--output':
      onbedfile=sys.argv[i+1];
      print('Output bed file:',onbedfile,file=sys.stderr);
    if sys.argv[i]=='-f' or sys.argv[i]=='--field':
      nfield=int(sys.argv[i+1]);
      print('Field:',nfield,file=sys.stderr);
    if sys.argv[i]=='-p' or sys.argv[i]=='--pairend':
      genpereads=True;
      pef=sys.argv[i+1].split(',');
      pemean=int(pef[0]);
      pestd=int(pef[1]);
      print('Generate paired-end reads with mean and std '+str(pemean)+','+str(pestd),file=sys.stderr);
  if sys.argv[i]=='-h' or sys.argv[i]=='--help':
    print(pydoc.render_doc(sys.modules[__name__]));
    sys.exit();
  if sys.argv[i]=='--stranded':
    stranded=True;

      

bedfile=sys.argv[-1];

# if no annotation file is specified, use uniform distri.
print('Assigning weights...',file=sys.stderr);
if len(allexp)==0:
  totalweight=0;
  for lines in open(bedfile):
    bedfield=lines.strip().split();
    allids[bedfield[3]]=0;
    totalweight+=1;
    allexp.append(totalweight);
    allidl.append(bedfield[3]);

# sampling process
print('Sampling...',file=sys.stderr);
for j in range(readcnt):
  k=random.random()*totalweight;
  sel=bisect.bisect_right(allexp,k);
  allids[allidl[sel]]=allids[allidl[sel]]+1;

# if no bias file specified, use uniform distrib

print('Total assigned reads:',sum(allids.values()),file=sys.stderr);

  
#debug info:
#for k in allidl:
#  print (k, allids[k]);

#sys.exit();

if onbedfile!="-":
  onfid=open(onbedfile,'w');
else:
  onfid=sys.stdout;


nlines=0;

totalgenreads=0;
# read bed file
for lines in open(bedfile):
  # update line counter
  nlines=nlines+1;
  if nlines %10000==1:
    print('Processing '+str(nlines)+' lines...',file=sys.stderr);
  # parse lines
  bedfield=lines.strip().split();
  if len(bedfield)!=12:
    print('Error: incorrect number of fields (should be 12)',file=sys.stderr);
    continue;
  if bedfield[5]=='+':
    direction=1;
  elif bedfield[5]=='-':
    direction=-1;
  else:
    print('Error: incorrect field in field[5] %s:' %bedfield[5],file=sys.stderr);
  if bedfield[3] not in allids:
    # the current id not found, continue
    continue;
  nreads=allids[bedfield[3]];
  if nreads<1:
    continue;
  # parse all segments
  fieldrange=(int(bedfield[1]),int(bedfield[2]));
  if bedfield[10][-1]==',':
    bedfield[10]=bedfield[10][:-1];
  if bedfield[11][-1]==',':
    bedfield[11]=bedfield[11][:-1];
  exonlen=[int(x) for x in bedfield[10].split(',')];
  exonstart=[int(x)+fieldrange[0] for x in bedfield[11].split(',')];
  # old code: for each possible position in the transcript, build its segments
  # for ne in range(len(exonlen)):
  #  for pos in range(exonstart[ne],exonstart[ne]+exonlen[ne]):
  # create a position
  totallen=sum(exonlen);
  # here, we randomly choose one position
  if genpereads==False:
    selrange=totallen-readlen+1;
  else:
    selrange=totallen-pemean+2*pestd;
  if selrange<1:
    if genpereads==False:
      print('Ignore annoatation',bedfield[3],'of length',totallen,'Reads:',allids[bedfield[3]],file=sys.stderr);
    else:
      print('Ignore annoatation',bedfield[3],'of length',totallen,'since its shorter than paired-end mean insert length. Reads:',allids[bedfield[3]],file=sys.stderr);
    continue;
  totalgenreads+=nreads;
  cumlen=[];cumlen.extend(exonlen);
  for i in range(1,len(cumlen)):
    cumlen[i]=cumlen[i]+cumlen[i-1];
  # for nun-uniform distribution, construct a new array for selection
  thistbweight=[];
  if len(posweight)!=0:
    kweight=0;
    for i in range(selrange):
      nfrac=i*100.0/selrange; # a value between 0-100
      nlower=int(math.floor(nfrac)); # 0-100
      nhigher=int(math.ceil(nfrac)); # 0-100
      if nhigher==nlower: nhigher=nlower+1;
      if nhigher<100:
        val=posweight[nlower]*(nfrac-nlower)+posweight[nhigher]*(nhigher-nfrac);
      else:
        val=posweight[99];
      kweight+=val;
      thistbweight.append(kweight);
  for t in range(nreads):
    if len(posweight)==0:
      tpos=random.choice(range(selrange));
    else:
      rd=random.random()*kweight;
      bsl=bisect.bisect_right(thistbweight,rd);
      # for reverse transcripts: flip the position
      if direction==-1:
        bsl=selrange-1-bsl;
      tpos=bsl;
    pos=tpos2pos(tpos,cumlen,exonstart);
    if genpereads==True:
      tpos2=tpos+int(random.normalvariate(pemean-readlen+1,pestd));
      pos2=tpos2pos(tpos2,cumlen,exonstart);
    # get the segments
    if True:
      (startrange,lenrange,status)=getSegs(pos,readlen,1,exonstart,exonlen);
      if status!=0:
        print('Status:',status,', pos:', pos,'out of',len(cumlen),file=sys.stderr);
        #pdb.set_trace();
        continue;
      # generate another pair
      if genpereads==True:
        (startrange2,lenrange2,status2)=getSegs(pos2,readlen,1,exonstart,exonlen);
        if status==1:
          print('Status:',status,', pos:', pos,'out of',len(cumlen),file=sys.stderr);
      if genpereads==False:
        lineid="%s_e_%d_%s_%d" % (bedfield[3],t,bedfield[0],pos);
      else:
        lineid="%s_e_%d_%s_%d/1" % (bedfield[3],t,bedfield[0],pos);
        lineid2="%s_e_%d_%s_%d/2" % (bedfield[3],t,bedfield[0],pos);
      # random direction
      if stranded==False or direction==0:
        thisdir=random.choice([1,-1]);
      else:
        thisdir=direction;
      writeBedline(onfid,lineid,bedfield[0],thisdir,startrange,lenrange);
      if genpereads==True:
        writeBedline(onfid,lineid2,bedfield[0],thisdir*(-1),startrange2,lenrange2);
    else:
      print(bedfield[0],file=sys.stdout);

#print('Pospool:');
#for k in sorted(pospool.keys()):
#  print(str(k)+":"+str(pospool[k]),end=",");
#print();          
        

print('Total '+str(nlines)+' lines...',file=sys.stderr);
print('Total '+str(totalgenreads)+' reads...',file=sys.stderr);
if onbedfile!="-":
  onfid.close();

