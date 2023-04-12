####Repeat summary figure cladewise####
args = commandArgs(trailingOnly=TRUE)
if (length(args)<1) {
  stop("This script is to plot repeat density cladewise\n\nPlease provide arguments like Rscript repeatdensity.r cladename\n\n", call.=FALSE)
}
library(stringr)
cl=as.character(args[1])
cname=str_to_title(cl)
fname=paste(cname,'_all_flps_combined',sep="")
a<-read.table(fname,header=F)
colnames(a)=c("xp","species","glength","rptstart","rptend","rptlen","rpt")
colpal=data.frame(aa = c("A","B","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","U","V","W","X","Y","Z","LM"), col=c("#771155", "#AA4488", "#CC99BB", "#114477", "#4477AA", "#77AADD", "#117777", "#44AAAA", "#77CCCC", "#117744", "#44AA77", "#88CCAA", "#777711", "#AAAA44", "#DDDD77", "#774411", "#AA7744", "#DDAA77", "#771122", "#AA4455", "#DD7788","orange","cyan","violet","black"))
fin=a
fin$newstart=((fin$rptend+fin$rptstart)/(2*fin$glength))*100
fin$purity=(fin$rptlen/((fin$rptend-fin$rptstart)+1))*100
for(pur in seq(from=0,to=100,by=10)){
final=fin[fin$purity>=pur,c("rpt","newstart")]
coun=dim(final)[1]
char <- strsplit(as.character(final$rpt), '')
final=data.frame(rpt=unlist(char), value=rep(final$newstart, sapply(char, FUN=length)))
dfp=data.frame()
for(rp in unique(final$rpt)){
b1=final[final$rpt==rp,]
if(dim(b1)[1] > 30) {
h=hist(b1$value,plot=F,breaks=seq(from=0,to=100,by=5))
h1=data.frame(rpt=rp,mid=h$mids,val=h$counts)
h1$normval=(h1$val/sum(h1$val))*100
dfp=rbind(dfp,h1)} }
ptitle=paste('LPS with >= ',pur,'% purity in ',cname,' clade across all genes',sep="")
mdt=hist(final$value,plot=F,breaks=seq(from=0,to=100,by=5))
rdf=data.frame(mid=mdt$mids,val=mdt$counts)
sldf=rdf
sldf$rpt="LM"
sldf$normval=(sldf$val/sum(sldf$val))*100
jpeg(paste('Repeat_position_',pur,'_purity_in_gene_in_',cname,'.jpeg',sep=''),width=28,height=18,units="in",res=300)
plot(sldf$mid,sldf$normval,ylim=c(0,max(dfp$normval)+10),xlim=c(0,100),type='l',lwd=4,col="black",main=ptitle,xlab="Normalized position in gene",ylab="",cex.main=2,cex.axis=1,cex.lab=1.5)
title(ylab="Frequency", line=2, cex.lab=1.5)
for(rp in unique(dfp$rpt)){
b1=dfp[dfp$rpt==rp,]
lines(b1$mid,b1$normval,type='l',lwd=2,col=colpal[colpal$aa==rp,2])}
tx=paste('Total count ',coun,' for the plot',sep="")
legend("topleft",legend=tx,bty='n',cex=1.5)
legend("topright",legend=colpal$aa,col=colpal$col,lty=1,bty='n',cex=1.2,seg.len=1.2,lwd=3.5,title="Repeats")
rect(0,2,100,8,border="transparent",col=rgb(0.5,0.5,0.5,alpha=0.5))
dev.off()
}
