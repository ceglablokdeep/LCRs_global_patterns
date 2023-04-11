a=read.table("exactstretch_diversity",header=F)
colnames(a)=c("clades","repeat","count")
library(stringr)
a$clades=str_to_title(a$clades)
finaldf <- data.frame(clades=character(), shd=numeric())
for(c in unique(a$clades)){
b=a[a$clades==c,]
b$V4=(b$count/sum(b$count))*(-log(b$count/sum(b$count)))
sd=sum(b$V4)
d=data.frame(c,sd)
finaldf=rbind(finaldf,d)
}
mf=data.frame()
for(clade in c("Afrotheria","Primates","Rodentia","Lagomorpha","Chiroptera","Artiodactyla","Perissodactyla","Carnivora","Marsupialia","Testudines","Aves","Squamata","Amphibia")){
m1=finaldf[finaldf$c==clade,]
mf=rbind(mf,m1)
}
finaldf=mf
jpeg("shannon_diversity_exactstretch.jpeg",width=23.2,height=13,units="in",res=300)
barplot(finaldf$sd,names=finaldf$c,ylab="Shannon diversity",ylim=c(0,max(finaldf$sd)),xlab="Clades",col="lightcyan",main="Shannon diversity of exact LPSs across different clades")
dev.off()
simpsondf <- data.frame(clades=character(), simd=numeric())
for(c in unique(a$clades)){
s=a[a$clades==c,]
s$V4=s$count*(s$count - 1)
simpd=round(1-(sum(s$V4)/(sum(s$count)*(sum(s$count)-1))),4)
simp=data.frame(c,simpd)
simpsondf=rbind(simpsondf,simp)}
sf=data.frame()
for(clade in c("Afrotheria","Primates","Rodentia","Lagomorpha","Chiroptera","Artiodactyla","Perissodactyla","Carnivora","Marsupialia","Testudines","Aves","Squamata","Amphibia")){
s1=simpsondf[simpsondf$c==clade,]
sf=rbind(sf,s1)
}
simpsondf=sf
jpeg("simpson_diversity_exactstretch.jpeg",width=23.2,height=13,units="in",res=300)
barplot(simpsondf$simpd,names=simpsondf$c,ylab="Simpson diversity",ylim=c(0,1),xlab="Clades",col="lightcyan",main="Simpson diversity of exact LPSs across different clades")
dev.off()
unnorrich <- data.frame(clades=character(), richness=numeric())
for(c in unique(a$clades)){
b=a[a$clades==c,]
numrep=dim(b)[1]
ric=data.frame(c,numrep)
unnorrich=rbind(unnorrich,ric)}
unno=data.frame()
for(clade in c("Afrotheria","Primates","Rodentia","Lagomorpha","Chiroptera","Artiodactyla","Perissodactyla","Carnivora","Marsupialia","Testudines","Aves","Squamata","Amphibia")){
s1=unnorrich[unnorrich$c==clade,]
unno=rbind(unno,s1)
}
unnorrich=unno
jpeg("unnormalized_LCR_richness_exactstretch.jpeg",width=23.2,height=13,units="in",res=300)
barplot(unnorrich$numrep,names=unnorrich$c,ylab="LCR richness",xlab="Clades",ylim=c(0,max(unnorrich$numrep)+2),col="#69b3a2",main="LCRs richness for exact LPSs across different clades")
dev.off()
numlps=read.table("numberoflps_exact.txt",header=F)
numlps$V1=str_to_title(numlps$V1)
norrich <- data.frame(clades=character(), richness=numeric())
for(c in unique(a$clades)){
rich=unnorrich[unnorrich$c==c,2]
numrep=numlps[numlps$V1==c,2]
slope=log(rich)/log(numrep)
rich=data.frame(c,slope)
norrich=rbind(norrich,rich)}
nono=data.frame()
for(clade in c("Afrotheria","Primates","Rodentia","Lagomorpha","Chiroptera","Artiodactyla","Perissodactyla","Carnivora","Marsupialia","Testudines","Aves","Squamata","Amphibia")){
s1=norrich[norrich$c==clade,]
nono=rbind(nono,s1)
}
norrich=nono
jpeg("normalized_LCR_richness_exactstretch.jpeg",width=23.2,height=13,units="in",res=300)
barplot(norrich$slope,names=norrich$c,ylab="LCR richness",xlab="Clades",ylim=c(0,max(norrich$slope)+0.2),col="#69b3a2",main="Normalized LCRs richness for exact LPSs across different clades")
dev.off()
norrich <- data.frame(clades=character(), richness=numeric(),numrepeats=numeric())
for(c in unique(a$clades)){
rich=unnorrich[unnorrich$c==c,2]
numrep=numlps[numlps$V1==c,2]
rich=data.frame(c,rich,numrep)
norrich=rbind(norrich,rich)}
nono=data.frame()
for(clade in c("Afrotheria","Primates","Rodentia","Lagomorpha","Chiroptera","Artiodactyla","Perissodactyla","Carnivora","Marsupialia","Testudines","Aves","Squamata","Amphibia")){
s1=norrich[norrich$c==clade,]
nono=rbind(nono,s1)
}
norrich=nono
jpeg("log_transformed_richness_exactstretch.jpeg",width=23.2,height=13,units="in",res=300)
b=lm(formula = log(norrich$rich) ~ log(norrich$numrep), data = norrich)
plot(log(norrich$numrep),log(norrich$rich),ylab="log(richness)",xlab="log(number of sequences)",col="blue3",pch=20,main="log transformed richness of exact LCRs and number of sequences across different clades")
abline(b,col="red",lty=2)
dev.off()
