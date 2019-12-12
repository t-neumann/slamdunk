#!/usr/bin/env python
"""
This script is used to add splicing variations from STDIN .BED file.

Usage: addvariation2splicingbed.py {OPTIONS} 

OPTIONS

ATTENTION: 

HISTORY
	01/09/2012 

"""
from __future__ import print_function;
import sys;
import subprocess;
import pydoc;
import os;
import random;
import bisect;
import math;
from getSegs import *;

import pdb;


errorrate=0.2;
onbedfile="-";


for i in range(len(sys.argv)):
  if sys.argv[i]=='-h':
    print(pydoc.render_doc(sys.modules[__name__]));
    sys.exit();
  if i<len(sys.argv)-1:
    if sys.argv[i]=='-o':
      onbedfile=sys.argv[i+1];
      print('Output bed file:',onbedfile,file=sys.stderr);



nlines=0;


for lines in sys.stdin:
  # update line counter
  nlines=nlines+1;
  if nlines %10000==1:
    print('Processing '+str(nlines)+' lines...',file=sys.stderr);
  # parse lines
  bedfield=lines.strip().split();
  if len(bedfield)!=12:
    print('Error: incorrect number of fields (should be 12)',file=sys.stderr);
    continue;
  if int(bedfield[9])!=2:
    print(lines,end='');
    continue;
  if bedfield[5]=='+':
    direction=1;
  elif bedfield[5]=='-':
    direction=-1;
  else:
    print('Error: incorrect field in field[5] %s:' %bedfield[5],file=sys.stderr);
  # parse all segments
  fieldrange=(int(bedfield[1]),int(bedfield[2]));
  if bedfield[10][-1]==',':
    bedfield[10]=bedfield[10][:-1];
  if bedfield[11][-1]==',':
    bedfield[11]=bedfield[11][:-1];
  exonlen=[int(x) for x in bedfield[10].split(',')];
  exonstart=[int(x)+fieldrange[0] for x in bedfield[11].split(',')];
  # get the segments
  if random.random()<errorrate:
    nshift=random.choice([-3,-2,-1,1,2,3]);
    ndir=random.choice([-1,1]);
    if ndir==-1:
      exonstart[0]=exonstart[0]+nshift;
    else:
      exonstart[1]=exonstart[1]+nshift;
  if True:
    # random direction
    writeBedline(sys.stdout,bedfield[3],bedfield[0],direction,exonstart,exonlen);
       

print('Total '+str(nlines)+' lines...',file=sys.stderr);

