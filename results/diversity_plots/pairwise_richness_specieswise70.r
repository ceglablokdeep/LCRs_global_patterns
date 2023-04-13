a=read.table("morethanseventy_diversity_specieswise",header=F)
colnames(a)=c("clades","species","repeat","count")
library(stringr)
a$clades=str_to_title(a$clades)
finaldfr<-data.frame(first=c(1:10))
for(c in unique(a$clades)){
finaldf<-data.frame(clades=character(), shd=numeric())
s1=a[a$clades==c,]
for(sp in unique(s1$species)){
s=s1[s1$species==sp,]
simpd=dim(s)[1]
sd=data.frame(simpd)
finaldf=rbind(finaldf,sd)
}
colnames(finaldf)=c
finaldfr=qpcR:::cbind.na(finaldfr,finaldf)
}
d=finaldfr[,-1]
d1=data.frame(finaldfr[,1])
for(clade in c("Afrotheria","Primates","Rodentia","Lagomorpha","Chiroptera","Artiodactyla","Perissodactyla","Carnivora","Marsupialia","Testudines","Aves","Squamata","Amphibia")){
d2=data.frame(d[,colnames(d)==clade])
colnames(d2)=as.character(clade)
d1=cbind(d1,d2)}
d=d1[,-1]
library(stringr)
str_to_title(colnames(d))->colnames(d)
species<-colnames(d)
len<-length(species)
totalcol=dim(d)[2]
ymin=c(min(d[,1],na.rm=T))
for ( i in 2:length(d)){
ymin[i]<-c(min(d[,i],na.rm=T))}
maxquan=0
P<- data.frame(id = character(0), p.value = numeric(0))
for ( i in seq(1,(totalcol-1),1)){
for ( j in seq(i+1,totalcol,1)){
if (species[i] == species[j]){
next
}
else
{
print(paste(species[i],species[j],sep=" "))
fg1<-d[,i]
bg1<-d[,j]
if (maxquan < quantile(fg1,0.75,na.rm=T)){maxquan=quantile(fg1,0.75,na.rm=T) }
wilcox.test(fg1,bg1,alternative = "two.sided")->K
P1<-data.frame(names=paste(species[i],"vs",species[j],sep="_"),p.value=K$p.value)
seg2=paste(species[i],"vs",species[j],sep="_")
P<-rbind(P,P1) }
}
}
options(scipen=0,"digits"=2)
P$qvalue<-p.adjust(P$p.value, method = "BH")
P$V1<-format(P$qvalue,scientific=T)
jpeg("Pairwise_boxplot_specieswise_richness_morethanseventy.jpeg",height=2500,width=3600)
par(mai=c(5,5,2,1),mgp=c(18,0.5,0.1))
boxplot(d,outline=F,main="Richness of LCRs across clades", xlab="Clades",ylab="Richness of LCRs", col="Lightsalmon", border="black", las=2, cex=4, cex.main=4, cex.axis=3, cex.lab=4, lwd=2, ylim=c(min(ymin),as.numeric(maxquan)+1+(2*(length(P[P$qvalue<0.05,3])))))
start=maxquan
for ( i in seq(1,len,1)){
for (j in seq(2,len,1)){
tryCatch({
if (species[i] == species[j]){
next
}
else
{
fg=as.character(species[i])
bg=as.character(species[j])
seg=paste(fg,"vs",bg,sep="_")
P[P$names==seg,3]-> qval
P[P$names==seg,4] -> qvalchar
if(qval > 0.05){
next
}
else
{
segments(i,start,j,start,lty=3,col="purple",lwd=2)
text(((i+j)/2),start+1,labels=qvalchar,cex=2)
}
}
start = start+2
}
, error=function(e){})}
}
dev.off()
jpeg("Boxplot_richness_specieswise_morethanseventy.jpeg",height=2500,width=3600)
par(mai=c(5,5,2,1),mgp=c(18,0.5,0.1))
boxplot(d,outline=F,main="Richness of LCRs across clades", xlab="Clades",ylab="Richness of LCRs", col="Lightsalmon", border="black", las=2, cex=4, cex.main=4, cex.axis=3, cex.lab=4, lwd=2)
dev.off()
####### for slope (normalized) #########
finaldfr <- data.frame(first=c(1:10))
for(c in unique(a$clades)){
finaldf <- data.frame(clades=character(), shd=numeric())
b1=a[a$clades==c,]
for(sp in unique(b1$species)){
b=b1[b1$species==sp,]
rich1=length(a[a$species==sp,4])
numrep=sum(a[a$species==sp,4])
slope=log(rich1)/log(numrep)
rich=data.frame(slope)
finaldf=rbind(finaldf,rich)
}
colnames(finaldf)=c
finaldfr=qpcR:::cbind.na(finaldfr,finaldf)
}
d=finaldfr[,-1]
d1=data.frame(finaldfr[,1])
for(clade in c("Afrotheria","Primates","Rodentia","Lagomorpha","Chiroptera","Artiodactyla","Perissodactyla","Carnivora","Marsupialia","Testudines","Aves","Squamata","Amphibia")){
d2=data.frame(d[,colnames(d)==clade])
colnames(d2)=as.character(clade)
d1=cbind(d1,d2)}
d=d1[,-1]
library(stringr)
str_to_title(colnames(d))->colnames(d)
species<-colnames(d)
len<-length(species)
totalcol=dim(d)[2]
ymin=c(min(d[,1],na.rm=T))
for ( i in 2:length(d)){
ymin[i]<-c(min(d[,i],na.rm=T))}
maxquan=0
P<- data.frame(id = character(0), p.value = numeric(0))
for ( i in seq(1,(totalcol-1),1)){
for ( j in seq(i+1,totalcol,1)){
if (species[i] == species[j]){
next
}
else
{
print(paste(species[i],species[j],sep=" "))
fg1<-d[,i]
bg1<-d[,j]
if (maxquan < quantile(fg1,0.75,na.rm=T)){maxquan=quantile(fg1,0.75,na.rm=T) }
wilcox.test(fg1,bg1,alternative = "two.sided")->K
P1<-data.frame(names=paste(species[i],"vs",species[j],sep="_"),p.value=K$p.value)
seg2=paste(species[i],"vs",species[j],sep="_")
P<-rbind(P,P1) }
}
}
options(scipen=0,"digits"=2)
P$qvalue<-p.adjust(P$p.value, method = "BH")
P$V1<-format(P$qvalue,scientific=T)
jpeg("Pairwise_boxplot_specieswise_richness_slope_morethanseventy.jpeg",height=2500,width=3600)
par(mai=c(5,5,2,1),mgp=c(18,0.5,0.1))
boxplot(d,outline=F,main="Slope of richness of LCRs across clades", xlab="Clades",ylab="Slope of richness of LCRs", col="Lightsalmon", border="black", las=2, cex=4, cex.main=4, cex.axis=3, cex.lab=4, lwd=2, ylim=c(min(ymin),as.numeric(maxquan)+1+(2*(length(P[P$qvalue<0.05,3])))))
start=maxquan
for ( i in seq(1,len,1)){
for (j in seq(2,len,1)){
tryCatch({
if (species[i] == species[j]){
next
}
else
{
fg=as.character(species[i])
bg=as.character(species[j])
seg=paste(fg,"vs",bg,sep="_")
P[P$names==seg,3]-> qval
P[P$names==seg,4] -> qvalchar
if(qval > 0.05){
next
}
else
{
segments(i,start,j,start,lty=3,col="purple",lwd=2)
text(((i+j)/2),start+1,labels=qvalchar,cex=2)
}
}
start = start+2
}
, error=function(e){})}
}
dev.off()
jpeg("Boxplot_specieswise_richness_slope_morethanseventy.jpeg",height=2500,width=3600)
par(mai=c(5,5,2,1),mgp=c(18,0.5,0.1))
boxplot(d,outline=F,main="Slope of richness of LCRs across clades", xlab="Clades",ylab="Slope of richness of LCRs", col="Lightsalmon", border="black", las=2, cex=4, cex.main=4, cex.axis=3, cex.lab=4, lwd=2)
dev.off()
########### slope of richness ###log-transformed richness and number of sequences #############
finaldf <- data.frame(clades=character(), species=character(), rich=numeric(), numrep=numeric())
for(c in unique(a$clades)){
b1=a[a$clades==c,]
for(sp in unique(b1$species)){
b=b1[b1$species==sp,]
rich1=length(a[a$species==sp,4])
numrep=sum(a[a$species==sp,4])
slope=log(rich1)/log(numrep)
rich=data.frame(c,sp,rich1,numrep)
finaldf=rbind(finaldf,rich)
}
}
col13rainbow= c("#771155", "#AA4488", "#CC99BB", "#114477", "#4477AA", "#77AADD", "#117777", "#44AAAA", "#77CCCC", "#117744", "#44AA77", "#88CCAA", "#777711")
groups=as.character(unique(finaldf$c))
## all clades in one-log plots
clade=groups[1]
jpeg("log_transformed_richness_specieswise_morethanseventy.jpeg",width=23.2,height=13,units="in",res=300)
b=lm(formula = log(finaldf$rich[finaldf$c==clade]) ~ log(finaldf$numrep[finaldf$c==clade]), data = finaldf)
plot(log(finaldf$numrep[finaldf$c==clade]),log(finaldf$rich[finaldf$c==clade]),xlim=c(5,9),ylim=c(4,6),ylab="log(richness)",xlab="log(number of LCRs)",col=col13rainbow[1],pch=20,main="log transformed richness of LCRs and number of LCRs across different clades")
abline(b,col=col13rainbow[1],lty=2)
for(clade in 2:length(groups)){
b=lm(formula = log(finaldf$rich[finaldf$c==groups[clade]]) ~ log(finaldf$numrep[finaldf$c==groups[clade]]), data = finaldf)
points(log(finaldf$numrep[finaldf$c==groups[clade]]),log(finaldf$rich[finaldf$c==groups[clade]]),ylab="log(richness)",xlab="log(number of LCRs)",col=col13rainbow[clade],pch=20,main="log transformed richness of LCRs and number of LCRs across different clades",cex=2)
abline(b,col=col13rainbow[clade],lty=2,lwd=2)
}
legend("topleft",legend=groups,fill=col13rainbow,col=col13rainbow,bty="n")
dev.off()
## cladewise log plots
for(clade in 1:length(groups)){
jpeg(paste(groups[clade],"_log_transformed_richness_specieswise_morethanseventy.jpeg",sep=""),width=23.2,height=13,units="in",res=300)
b=lm(formula = log(finaldf$rich[finaldf$c==groups[clade]]) ~ log(finaldf$numrep[finaldf$c==groups[clade]]), data = finaldf)
plot(log(finaldf$numrep[finaldf$c==groups[clade]]),log(finaldf$rich[finaldf$c==groups[clade]]),ylab="log(richness)",xlab="log(number of LCRs)",col=col13rainbow[clade],pch=20,main=paste("log transformed richness of LCRs and number of LCRs in ", groups[clade],sep=""),cex=2)
abline(b,col=col13rainbow[clade],lty=2,lwd=2)
legend("topleft",legend=groups[clade],fill=col13rainbow[clade],col=col13rainbow[clade],bty="n")
dev.off()
}
for(clade in 1:length(groups)){
jpeg(paste(groups[clade],"_log_transformed_richness_specieswise_morethanseventy.jpeg",sep=""),width=23.2,height=13,units="in",res=300)
b=lm(formula = log(finaldf$rich[finaldf$c==groups[clade]]) ~ log(finaldf$numrep[finaldf$c==groups[clade]]), data = finaldf)
plot(log(finaldf$numrep[finaldf$c==groups[clade]]),log(finaldf$rich[finaldf$c==groups[clade]]),ylab="log(richness)",xlab="log(number of LCRs)",col=col13rainbow[clade],pch=20,main=paste("log transformed richness of LCRs and number of LCRs in ", groups[clade],sep=""),cex=2)
abline(b,col=col13rainbow[clade],lty=2,lwd=2)
legend("topleft",legend=groups[clade],fill=col13rainbow[clade],col=col13rainbow[clade],bty="n")
dev.off()
}
