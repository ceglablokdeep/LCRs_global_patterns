####in R
for(clade in c("Afrotheria","Primates","Rodentia","Chiroptera","Artiodactyla","Perissodactyla","Carnivora","Marsupialia","Testudines","Aves","Squamata","Amphibia")){
print(clade)
fname=paste(clade,'_contingency',sep="")
a=read.table(fname,header=F)
colnames(a)=c("#pss_genes","#no_pss")
row.names(a)=c("#rpt_genes", "#no_rpt_genes")
imgname=paste('mosaicplot_of_pss_rpt_for_',clade,'.jpeg',sep="")
jpeg(imgname,width=9,height=6,units="in",res=300)
mosaicplot(a,col=c("slategray2","skyblue2"),main=paste('Occurrence of LCR containing genes and PSS containing genes in ',clade,sep=""))
dev.off()
fg=fisher.test(a, alternative = "g")
fl=fisher.test(a, alternative = "l")
ft=fisher.test(a, alternative = "t")
fisherout=paste(clade,as.numeric(fg$estimate),as.numeric(fg$p.value),as.numeric(fl$p.value),as.numeric(ft$p.value),sep="\t")
write.table(fisherout,file="fisher_contingency_out.txt",append=T,quote=F,sep="\t",col.names=F,row.names=F)
}
