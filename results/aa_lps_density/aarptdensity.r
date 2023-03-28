library(stringr)
a1=read.table("all_clades_flps_out_combined",header=F)
colnames(a1)=c("xp","species","glength","rptstart","rptend","rptlen","rpt","clade")
a1$newstart=((a1$rptend+a1$rptstart)/(2*a1$glength))*100
a1$purity=(a1$rptlen/((a1$rptend-a1$rptstart)+1))*100
colpal=data.frame(clade = c("Afrotheria","Amphibia","Artiodactyla","Aves","Carnivora","Chiroptera","Lagomorpha","Marsupialia","Perissodactyla","Primates","Rodentia","Squamata","Testudines","LM"), col=c("#771155", "#AA4488", "#CC99BB", "#114477", "#4477AA", "#77AADD", "#117777", "#44AAAA", "#77CCCC", "#117744", "#44AA77", "#88CCAA", "#777711","black"))
char <- strsplit(as.character(a1$rpt), '')
final=data.frame(rpt=unlist(char), nstart=rep(a1$newstart, sapply(char, FUN=length)),cl=rep(a1$clade, sapply(char, FUN=length)),pur=rep(a1$purity, sapply(char, FUN=length)))
####First loop: for purity
for(pur in seq(from=0,to=100,by=10)){
a=final[final$pur>=pur,]
####Second loop: for repeats
for(aa in unique(a$rpt)){
spa=a[a$rpt==aa,]
coun=dim(spa)[1]
##firstly making the general line of all clades together for the AA type
mdt=hist(spa$nstart,plot=F,breaks=seq(from=0,to=100,by=5))
rdf=data.frame(mid=mdt$mids,val=mdt$counts)
sldf=rdf
sldf$cl="LM"
sldf$normval=(sldf$val/sum(sldf$val))*100
####Third loop: for clades
for(c in unique(spa$cl)){
spacl=spa[spa$cl==c,]
htem=hist(spacl$nstart,plot=F,breaks=seq(from=0,to=100,by=5))
cldf=data.frame(mid=htem$mids,val=htem$counts)
acldf=cldf
acldf$cl=c
acldf$normval=(acldf$val/sum(acldf$val))*100
sldf=rbind(sldf,acldf)}
ptitle=paste('Amino acid ',aa,' LPS distribution for >=',pur,'% purity in all clades across all genes',sep="")
jpeg(paste('AA_repeat_position_',aa,'_for_',pur,'_purity_in_all_genes.jpeg',sep=''),width=28,height=18,units="in",res=300)
plot(1,type="n",ylim=c(0,max(sldf$normval)+3),xlim=c(0,100),main=ptitle,xlab="Normalized position in gene", ylab="", cex.main=4, cex.axis=1, cex.lab=1.5)
title(ylab="Frequency", line=2, cex.lab=1.5)
for(c in unique(sldf$cl)){
j1=sldf[sldf$cl==c,]
lines(j1$mid,j1$normval,type='l',lwd=1.5,col=colpal[colpal$clade==c,2])}
j1=sldf[sldf$cl=="LM",]
lines(j1$mid,j1$normval,type='l',lwd=2.5,col=colpal[colpal$clade=="LM",2])
tx=paste('Total count ',coun,' for the plot',sep="")
legend("topleft",legend=tx,bty='n',cex=1.2)
legend("topright",legend=colpal$clade,col=colpal$col,lty=1,bty='n',cex=1.2,seg.len=1.2,lwd=2.5,title="Clades")
rect(0,2,100,8,border="transparent",col=rgb(0.5,0.5,0.5,alpha=0.5))
dev.off()
}
}
