##in R
a=read.table("flps_all_output_with_xm_gene_sp_clade_phylo_arranged")
b=a[,c(1,6,7,8)]
b$rplen=(b$V7-b$V6)+1
b$rppur=b$V8/b$rplen
library(vioplot)
library(stringr)
colpal=c("#771155", "#AA4488", "#CC99BB", "#114477", "#4477AA", "#77AADD", "#117777", "#44AAAA", "#77CCCC", "#117744", "#44AA77", "#88CCAA", "#777711", "#AAAA44", "#DDDD77", "#774411", "#AA7744", "#DDAA77", "#771122", "#AA4455", "#DD7788","orange","cyan","violet")
colpal1=c("#777711", "#AAAA44", "#DDDD77", "#774411", "#AA7744", "#DDAA77", "#771122", "#AA4455", "#DD7788","orange","cyan","violet", "#AA7744")
cnames=c("Afrotheria","Primates","Rodentia","Lagomorpha","Chiroptera","Artiodactyla","Perissodactyla","Carnivora","Marsupialia","Testudines","Aves","Squamata","Amphibia")
for(pur in seq(from=0,to=1,by=0.1)){
print(pur)
pp=pur*100
c=b[b$rppur>=pur,c(1,5)]
c$V1=str_to_title(c$V1)
jpeg(paste('lps_length_distribution_',pur,'_in_all_genes.jpeg',sep=''),width=32,height=18,units="in",res=300)
par(mar=c(7,7,5,2))
mt=paste("LPS length with more than or equal to ",pp,"% purity across different clades",sep="")
vioplot(c$rplen[c$V1=="Afrotheria"],c$rplen[c$V1=="Primates"],c$rplen[c$V1=="Rodentia"],c$rplen[c$V1=="Lagomorpha"],c$rplen[c$V1=="Chiroptera"],c$rplen[c$V1=="Artiodactyla"],c$rplen[c$V1=="Perissodactyla"],c$rplen[c$V1=="Carnivora"],c$rplen[c$V1=="Marsupialia"],c$rplen[c$V1=="Testudines"],c$rplen[c$V1=="Aves"],c$rplen[c$V1=="Squamata"],c$rplen[c$V1=="Amphibia"],col = colpal,xlab = "",ylab = "",border=colpal,main=mt,cex.names=1.8,cex.main=3,cex.axis=1.8,names=cnames)
pt1=paste('Afrotheria ( n =',length(c$rplen[c$V1=="Afrotheria"]),')')
pt2=paste('Primates ( n =',length(c$rplen[c$V1=="Primates"]),')')
pt3=paste('Rodentia ( n =',length(c$rplen[c$V1=="Rodentia"]),')')
pt4=paste('Lagomorpha ( n =',length(c$rplen[c$V1=="Lagomorpha"]),')')
pt5=paste('Chiroptera ( n =',length(c$rplen[c$V1=="Chiroptera"]),')')
pt6=paste('Artiodactyla ( n =',length(c$rplen[c$V1=="Artiodactyla"]),')')
pt7=paste('Perissodactyla ( n =',length(c$rplen[c$V1=="Perissodactyla"]),')')
pt8=paste('Carnivora ( n =',length(c$rplen[c$V1=="Carnivora"]),')')
pt9=paste('Marsupialia ( n =',length(c$rplen[c$V1=="Marsupialia"]),')')
pt10=paste('Testudines ( n =',length(c$rplen[c$V1=="Testudines"]),')')
pt11=paste('Aves ( n =',length(c$rplen[c$V1=="Aves"]),')')
pt12=paste('Squamata ( n =',length(c$rplen[c$V1=="Squamata"]),')')
pt13=paste('Amphibia ( n =',length(c$rplen[c$V1=="Amphibia"]),')')
legend("topright",c(pt1,pt2,pt3,pt4,pt5,pt6,pt7,pt8,pt9,pt10,pt11,pt12,pt13),fill = colpal,border = NA,cex=2.5,box.col = "white",bg = "white")
title(xlab="Clades",ylab="Length of LPS (bp)",cex.lab=2.2,line=3.5)
dev.off()
print(pur)
}
