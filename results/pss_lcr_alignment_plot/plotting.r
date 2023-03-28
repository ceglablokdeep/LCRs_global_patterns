s=read.table("sitenorm.txt")
r=read.table("rpt_norm.txt")
a=read.table("final_alndensity.txt")
t=cbind(a,hist(s$V2,breaks=100,plot=F)$density)
t=cbind(t,hist(r$V2,breaks=100,plot=F)$density)
colnames(t)=c("breaks","Alignment_coverage","Positive_sites","Repeats_location")
jpeg("Positive_repeat_coverage_distribution.jpeg",width=28,height=18,units="in",res=300)
par(mar=c(5.1,5.1,5,2.1))
plot(t$Alignment_coverage/20,type='l',col="darkorchid4",lwd=3,ylim=c(0,0.05),xlab="Normalized gene position",ylab="Density",main="Density distribution on normalized gene position",cex.main=2.5,cex.lab=2,cex.axis=2)
lines(t$Positive_sites,lwd=3,col="hotpink2")
lines(t$Repeats_location,lwd=3,col="salmon2")
legend(80,0.03,box.col="white",bg ="white", box.lwd =0,legend=c("Alignment coverage/20", "Positive sites","Repeats location"),fill = c("darkorchid4","hotpink2","salmon2"),cex=2)
dev.off()
library(corrplot)
b3=t[,-1]
m=cor(b3,method="kendall")
testRes = cor.mtest(b3, conf.level = 0.95)
jpeg("Positive_repeat_coverage_corrplot.jpeg",width=10,height=11,units="in",res=300)
par(mar=c(1,3,10,2.1))
corrplot(m,p.mat = testRes$p, method="shade",title="",order ="hclust",hclust.method="complete",insig="blank",addCoef.col ='white',sig.level=0.05,tl.col = "black",tl.cex=2)
dev.off()
