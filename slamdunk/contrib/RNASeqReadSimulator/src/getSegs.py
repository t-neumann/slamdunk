#!/usr/bin/env python3

import bisect;

def getSegs(pos,seglen,direction,exonstart,exonlen):
  """
  This function returns the corresponding segstart and seg length (as in the list)
  parameters:
  pos: the position to be queried
  seglen: the length of the segment
  direction: 1 (forward), -1 (backward)
  exonstart,exonlen: the exon positions
  Return value: (segstart,seglength,status)
  status is used to indicate the possible errors
  0: normal
  1: the position is not in the exon range
  2: the segment exceeds the boundary
  """
  segstart=[];
  seglength=[];
  # find the corresponding exon which includes pos
  posexonid=-1;
  status=0;
  for ne in range(len(exonstart)):
    if pos in range(exonstart[ne],exonstart[ne]+exonlen[ne]):
      posexonid=ne;
      break;
  if posexonid==-1:
    status=1;
    return (segstart,seglength,status);
  if direction==1:
    while seglen>0:
      lentoadd=min(seglen,exonlen[posexonid]+exonstart[posexonid]-pos);
      segstart+=[pos];
      seglength+=[lentoadd];
      posexonid=posexonid+1;
      seglen=seglen-lentoadd;
      if posexonid>=len(exonstart):
        if seglen>0:
          status=2;
        return (segstart,seglength,status);
      pos=exonstart[posexonid];
  if direction==-1:
    while seglen>0:
      lentoadd=min(seglen,pos-exonstart[posexonid]+1);
      segstart.insert(0,pos-lentoadd+1); # insert in the front
      seglength.insert(0,lentoadd);
      posexonid=posexonid-1;
      seglen=seglen-lentoadd;
      if posexonid<0:
        if seglen>0:
          status=2;
        return (segstart,seglength,status);
      pos=exonstart[posexonid]+exonlen[posexonid]-1;
  return (segstart,seglength,status);
    
def tpos2pos(tpos,cumlen,exonstart):
  """
  Convertion from coordinates in a transcript to coordinates in a reference.
  Need to provide exon start position and the cumulate exon length as input.
  """        
  selseg=bisect.bisect_right(cumlen,tpos);
  # if the position is exceeding the boundary, set as the last position of the boundary
  if selseg>=len(cumlen):
    selseg=len(cumlen)-1;
    tpos=cumlen[-1]-1;
  if selseg>0:
    pos=exonstart[selseg]+(tpos-cumlen[selseg-1]);
  else:
    pos=exonstart[selseg]+tpos;
  return pos;
  
def writeBedline(fid,lineid,chromosome,direction,startrange,lenrange):
  """
  Write one line in .bed file.
  Need to provide information of chromosome, id, direction, segment starts and segment lengths
  """
  # skip if startrange is malformed
  if not startrange:
    return None
  bedrange=(startrange[0],startrange[-1]+lenrange[-1]);
  startrange=[i-startrange[0] for i in startrange];
  # directions
  if direction==1:
    direction='+';
  elif direction==-1:
    direction='-';
  #write line
  fid.write(chromosome + '\t' # 0th, chromosome
      + str(bedrange[0])+'\t'+str(bedrange[1])+'\t' # 1-2th, start and and
      + lineid + '\t' # 3th, id
      + '0\t'+ direction +'\t' # 4th, 5th, 0 and direction
      + str(bedrange[0])+'\t'+str(bedrange[1])+'\t' # 6-7th, same as 1-2
      + '0\t'+str(len(startrange))+'\t' # 8th, 0; 9th, number of segments
      + ''.join([str(i)+',' for i in lenrange]) + '\t' # 10th, length
      + ''.join([str(i)+',' for i in startrange]) +'\t' # 11th, start position
      +'\n');


