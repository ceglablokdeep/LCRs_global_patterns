##run the script as follows: Rscript pruning.r treename listofspeciestokeep.txt outputname
args=commandArgs(trailingOnly=TRUE)
t=args[1]
l=args[2]
o=args[3]
library(phytools)
a<-read.tree(t)
k<-read.table(l)
keep<-as.character(k$V1)
b<-keep.tip(a,keep)
write.tree(b,file=paste(o,".nwk",sep=''))
