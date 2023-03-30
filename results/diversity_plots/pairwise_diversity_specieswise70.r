a=read.table("morethanseventy_diversity_specieswise",header=F)
colnames(a)=c("clades","species","repeat","count")
library(stringr)
a$clades=str_to_title(a$clades)
#install.packages("qpcR")
library("qpcR")
####for Shannon diversity####
finaldfr <- data.frame(first=c(1:10))
for(c in unique(a$clades)){
finaldf <- data.frame(clades=character(), shd=numeric())
b1=a[a$clades==c,]
for(sp in unique(b1$species)){
b=b1[b1$species==sp,]
b$V4=(b$count/sum(b$count))*(-log(b$count/sum(b$count)))
sd1=sum(b$V4)
sd=data.frame(sd1)
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
if (maxquan < quantile(fg1,0.75,na.rm=T)){maxquan=quantile(fg1,0.75,na.rm=T)}
wilcox.test(fg1,bg1,alternative = "two.sided")->K
P1<-data.frame(names=paste(species[i],"vs",species[j],sep="_"),p.value=K$p.value)
seg2=paste(species[i],"vs",species[j],sep="_")
P<-rbind(P,P1) }
}
}
options(scipen=0,"digits"=2)
P$qvalue<-p.adjust(P$p.value, method = "BH")
P$V1<-format(P$qvalue,scientific=T)
jpeg("Pairwise_boxplot_specieswise_Shannon_diversity_morethanseventy.jpeg",height=2500,width=3600)
par(mai=c(5,5,2,1),mgp=c(18,0.5,0.1))
boxplot(d,outline=F,main="Shannon diversity of repeats across different clades", xlab="Clades",ylab="Shannon diversity", col="Lightsalmon", border="black", las=2, cex=4, cex.main=4, cex.axis=3, cex.lab=4, lwd=2, ylim=c(min(ymin),as.numeric(maxquan)+1+(2*(length(P[P$qvalue<0.05,3])))))
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
jpeg("Shannon_diversity_specieswise_boxplot_morethanseventy.jpeg",height=2500,width=3600)
par(mai=c(5,5,2,1),mgp=c(18,0.5,0.1))
boxplot(d,outline=F,main="Shannon diversity of repeats across different clades", xlab="Clades",ylab="Shannon diversity", col="Lightsalmon", border="black", las=2, cex=4, cex.main=4, cex.axis=3, cex.lab=4, lwd=2)
dev.off()
####for Simpson diversity####
finaldfr<-data.frame(first=c(1:10))
for(c in unique(a$clades)){
finaldf<-data.frame(clades=character(), shd=numeric())
s1=a[a$clades==c,]
for(sp in unique(s1$species)){
s=s1[s1$species==sp,]
s$V4=s$count*(s$count - 1)
simpd=round(1-(sum(s$V4)/(sum(s$count)*(sum(s$count)-1))),4)
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
jpeg("Pairwise_boxplot_specieswise_Simpson_diversity_morethanseventy.jpeg",height=2500,width=3600)
par(mai=c(5,5,2,1),mgp=c(18,0.5,0.1))
boxplot(d,outline=F,main="Simpson diversity of repeats across different clades", xlab="Clades",ylab="Simpson diversity", col="Lightsalmon", border="black", las=2, cex=4, cex.main=4, cex.axis=3, cex.lab=4, lwd=2, ylim=c(min(ymin),as.numeric(maxquan)+1+(2*(length(P[P$qvalue<0.05,3])))))
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
jpeg("Simpson_diversity_specieswise_boxplot_morethanseventy.jpeg",height=2500,width=3600)
par(mai=c(5,5,2,1),mgp=c(18,0.5,0.1))
boxplot(d,outline=F,main="Simpson diversity of repeats across different clades", xlab="Clades",ylab="Simpson diversity", col="Lightsalmon", border="black", las=2, cex=4, cex.main=4, cex.axis=3, cex.lab=4, lwd=2)
dev.off()
