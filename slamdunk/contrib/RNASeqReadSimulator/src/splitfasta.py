#!/usr/bin/env python
"""
Split fasta files including paired-end reads into two separate files
"""
from __future__ import print_function;
import sys;
import re;

outfile1="";
outfile2="";

for i in range(len(sys.argv)):
  if sys.argv[i]=="-o":
    outfile1=sys.argv[i+1]+"_1.fa";
    outfile2=sys.argv[i+1]+"_2.fa";

if outfile1=="":
  sys.exit(-1);

ofid1=open(outfile1,"w");
ofid2=open(outfile2,"w");

isleft=True;
for lines in sys.stdin:
  if lines[0]=='>':
    if lines.strip()[-1]=='1':
      isleft=True;
    else:
      isleft=False;
  lines=re.sub("/[12]","",lines);
  if isleft:
    print(lines,file=ofid1,end='');
  else:
    print(lines,file=ofid2,end='');





ofid1.close();
ofid2.close();
