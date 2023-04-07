r=read.table("all_flps_output_combined",header=F)
r$rpnorm=((r$V4+r$V3)*100)/(2*(r$V2))
r70=r[(r$V5/((r$V4-r$V3)+1))>0.7,]
rpure=r[(r$V5==((r$V4-r$V3)+1)),]
r=r[,c('V1','rpnorm')]
r$V3="Repeat"
b=r[,c(3,2)]
nrp=dim(b)[1]
p2 <- hist(b$rpnorm,breaks=19,freq=F)
pr=data.frame(V1=p2$mids,V2=p2$density)
#####more than 70 %###########
r70=r70[,c('V1','rpnorm')]
r70$V3="Repeat"
b70=r70[,c(3,2)]
nrp70=dim(b70)[1]
p270 <- hist(b70$rpnorm,breaks=19,freq=F)
pr70=data.frame(V1=p270$mids,V2=p270$density)
########pure repeats########
rpure=rpure[,c('V1','rpnorm')]
rpure$V3="Repeat"
bpure=rpure[,c(3,2)]
nrppure=dim(bpure)[1]
p2pure <- hist(bpure$rpnorm,breaks=19,freq=F)
prpure=data.frame(V1=p2pure$mids,V2=p2pure$density)
xlim <- range(0,100)
ylim <- range(0,max(pr$V2,pr70$V2,prpure$V2)+0.01)
jpeg(paste('overlap_of_repeatgrads_in_all_clades_in_all_genes.jpeg',sep=''),width=28,height=18,units="in",res=300)
par(mar=c(6,7,5,2))
mtitle=paste("Distribution of simple sequences in all clades on normalized gene position",sep=" ")
plot(1, type = "n", xlab = "",ylab = "", xlim = xlim, ylim = ylim,main=mtitle,cex.axis=2,cex.main=2.5,xaxt='n')
for (i in seq(0,100,5)){abline(v=i,lty=2,lwd=0.6,col="gray50")}
axis(side=1, at=seq(2.5,97.5, 5), labels=c("(0,5]","(5,10]","(10,15]","(15,20]","(20,25]","(25,30]","(30,35]","(35,40]","(40,45]","(45,50]","(50,55]","(55,60]","(60,65]","(65,70]","(70,75]","(75,80]","(80,85]","(85,90]","(90,95]","(95,100]"),cex.axis=1.3)
for(i in 1:dim(pr)[1]){
r1=pr[i,]
x1=as.numeric(r1[1,1])-2.5
x2=as.numeric(r1[1,1])-0.833333333
y1=0
y2=as.numeric(r1[1,2])
rect(x1,y1,x2,y2,border="transparent",col=rgb(0,0,1,1/3))
}
for(i in 1:dim(pr70)[1]){
r1=pr70[i,]
x1=as.numeric(r1[1,1])-0.833333333
x2=as.numeric(r1[1,1])+0.833333333
y1=0
y2=as.numeric(r1[1,2])
rect(x1,y1,x2,y2,border="transparent",col=rgb(1,0,0,1/3))
}
for(i in 1:dim(prpure)[1]){
r1=prpure[i,]
x1=as.numeric(r1[1,1])+0.833333333
x2=as.numeric(r1[1,1])+2.5
y1=0
y2=as.numeric(r1[1,2])
rect(x1,y1,x2,y2,border="transparent",col=rgb(0,1,0,1/3))
}
lines(density(b$rpnorm),col=rgb(0,0,1,1/2),lwd=3)
lines(density(b70$rpnorm),col=rgb(1,0,0,1/2),lwd=3)
lines(density(bpure$rpnorm),col=rgb(0,1,0,1/2),lwd=3)
title(ylab="Proportion", line=3.5, cex.lab=2.3)
title(xlab="Normalized position in gene", line=3.5, cex.lab=2.3)
pr1=paste('Compositionally-biased regions ( n =',nrp,')')
pr70t=paste('Compositionally-biased regions of > 70% purity ( n =',nrp70,')')
prp=paste('Pure amino acid repeats ( n =',nrppure,')')
legend(x=50,y=0.032,c(pr1,pr70t,prp),fill = c(rgb(0,0,1,1/3), rgb(1,0,0,1/3),rgb(0,1,0,1/3)),border = NA,cex=2.5,box.col = "white",bg = "white")
dev.off()
