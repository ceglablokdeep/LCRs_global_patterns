args = commandArgs(trailingOnly=TRUE)
if (length(args)<2) {
  stop("Please provide arguments like Rscript overlap_histo.r afrotheria_all_flps_combined cladename [the values in file should look like this: XP_006877394.1_Chrysochloris_asiatica 519 395 423 9 {P}]", call.=FALSE)
}
library(stringr)
r=read.table(args[1],header=F)
cname=args[2]
clade=str_to_title(cname)
s=read.table("sitenorm_clades.txt")
r$rpnorm=((r$V4+r$V3)*100)/(2*(r$V2))
r=r[,c('V1','rpnorm')]
s=s[s$V1==cname,]
r$V3="Repeat"
s$V4="Positive-site"
a=s[,c(4,3)]
b=r[,c(3,2)]
nps=dim(a)[1]
nrp=dim(b)[1]
p1 <- hist(a$V3,breaks=19,freq=F) 
p2 <- hist(b$rpnorm,breaks=19,freq=F)
ps=data.frame(V1=p1$mids,V2=p1$density)
pr=data.frame(V1=p2$mids,V2=p2$density)
xlim <- range(0,100)
ylim <- range(0,max(ps$V2,pr$V2)+0.01)
jpeg(paste('overlap_of_lcr_and_pssites_in_',cname,'_in_all_genes.jpeg',sep=''),width=28,height=18,units="in",res=300)
par(mar=c(6,7,5,2))
mtitle=paste("Distribution of low-complexity regions and positively selected sites in",clade,"on normalized gene position",sep=" ")
plot(1, type = "n", xlab = "",ylab = "", xlim = xlim, ylim = ylim,main=mtitle,cex.axis=2,cex.main=2.5,xaxt='n')
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
lines(density(b$rpnorm),col=rgb(1,0,0,1/2),lwd=3)
lines(density(a$V3),col=rgb(0,0,1,1/2),lwd=3)
title(ylab="Proportion", line=3.5, cex.lab=2.3)
title(xlab="Normalized position in gene", line=3, cex.lab=2.3)
pt1=paste('Positively-selected sites ( n =',nps,')')
pr1=paste('Low-complexity regions ( n =',nrp,')')
legend(x=50,y=0.04,c(pt1,pr1),fill = c(rgb(0,0,1,1/3), rgb(1,0,0,1/3)),border = NA,cex=2.5,box.col = "white",bg = "white")
dev.off()
