import sys
import re
filename=sys.argv[1]
outname=sys.argv[2]
fn = open(filename, 'r')
h = fn.readlines()
for i in range(0,len(h),2):
  seqrange=i+1
  header=h[i]
  line=h[seqrange]
  header=header.strip()
  line=line.strip()
  print(header)
  #print(line)
  #print(len(line))
  if line.startswith("ATG") and len(line)%3==0 and line.endswith(("TAA","TAG","TGA")) and re.match('^[ATGCN]*$', line):
    line=line[:-3]
    n=3
    codonlist=[]
    for index in range(0,len(line),n):
      codon=line[index:index+n]
      codonlist.append(codon)
    if 'TAA' not in codonlist and 'TAG' not in codonlist and 'TGA' not in codonlist:
      print(header, file=open(str(outname)+"_orf.txt", "a"))
      print(line, file=open(str(outname)+"_orf.txt", "a"))
