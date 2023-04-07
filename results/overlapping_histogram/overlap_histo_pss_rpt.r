r=read.table("all_flps_norm_pos")
s=read.table("sitenorm.txt")
r$V3="Repeat"
s$V3="Positive-site"
a=s[,c(3,2)]
b=r[,c(3,2)]
nps=dim(a)[1]
nrp=dim(b)[1]
p1 <- hist(a$V2,breaks=19,freq=F) 
p2 <- hist(b$V2,breaks=19,freq=F)
ps=data.frame(V1=p1$mids,V2=p1$density)
pr=data.frame(V1=p2$mids,V2=p2$density)
xlim <- range(0,100)
ylim <- range(0,max(ps$V2,pr$V2)+0.01)
jpeg(paste('overlap_of_lcr_and_pssites_in_all_genes.jpeg',sep=''),width=28,height=18,units="in",res=300)
par(mar=c(6,7,5,2))
plot(1, type = "n", xlab = "",ylab = "", xlim = xlim, ylim = ylim,main="Distribution of low-complexity regions and positively selected sites on normalized gene position",cex.axis=2,cex.main=3,xaxt='n')
for (i in seq(0,100,5)){abline(v=i,lty=2,lwd=0.6,col="gray50")}
axis(side=1, at=seq(2.5,97.5, 5), labels=c("(0,5]","(5,10]","(10,15]","(15,20]","(20,25]","(25,30]","(30,35]","(35,40]","(40,45]","(45,50]","(50,55]","(55,60]","(60,65]","(65,70]","(70,75]","(75,80]","(80,85]","(85,90]","(90,95]","(95,100]"),cex.axis=1.3)
for(i in 1:dim(ps)[1]){
r1=ps[i,]
x1=as.numeric(r1[1,1])-2.5
x2=as.numeric(r1[1,1])
y1=0
y2=as.numeric(r1[1,2])
rect(x1,y1,x2,y2,border="transparent",col=rgb(0,0,1,1/3))
}
for(i in 1:dim(pr)[1]){
r1=pr[i,]
x1=as.numeric(r1[1,1])
x2=as.numeric(r1[1,1])+2.5
y1=0
y2=as.numeric(r1[1,2])
rect(x1,y1,x2,y2,border="transparent",col=rgb(1,0,0,1/3))
}
lines(density(b$V2),col=rgb(1,0,0,1/2),lwd=3)
lines(density(a$V2),col=rgb(0,0,1,1/2),lwd=3)
title(ylab="Proportion", line=3.5, cex.lab=2.3)
title(xlab="Normalized position in gene", line=3, cex.lab=2.3)
pt1=paste('Positively-selected sites ( n =',nps,')')
pr1=paste('Low-complexity regions ( n =',nrp,')')
legend(x=50,y=0.04,c(pt1,pr1),fill = c(rgb(0,0,1,1/3), rgb(1,0,0,1/3)),border = NA,cex=2.5,box.col = "white",bg = "white")
dev.off()
