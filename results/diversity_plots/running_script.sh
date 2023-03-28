for clades in afrotheria amphibia artiodactyla aves carnivore chiroptera lagomorpha marsupials perissodactyla primates rodents squamata testudines
do
cut -f1 "$clades"_morethanseventy_combined|sort|uniq -c|awk '{print $2,$1}' OFS='\t' > "$clades"_numberoflps.txt
done
###############################
##calculating diversity cladewise for >70% LPS
awk '{print $1,$11}' OFS='\t' morethanseventy_combined > morethanseventy_diversity
for clades in afrotheria amphibia artiodactyla aves carnivore chiroptera lagomorpha marsupials perissodactyla primates rodents squamata testudines
do
grep "$clades" morethanseventy_diversity|sort -k2|uniq -c|awk '{print $2,$3,$1}' OFS='\t'|sort -k3nr >> temp.txt
done
mv temp.txt morethanseventy_diversity
################Plotting in R#############
echo 'a=read.table("morethanseventy_diversity",header=F)' > diversity_richness70.r
echo 'colnames(a)=c("clades","repeat","count")' >> diversity_richness70.r
echo 'library(stringr)' >> diversity_richness70.r
echo 'a$clades=str_to_title(a$clades)' >> diversity_richness70.r
####for Shannon diversity####
echo 'finaldf <- data.frame(clades=character(), shd=numeric())' >> diversity_richness70.r
echo 'for(c in unique(a$clades)){' >> diversity_richness70.r
echo 'b=a[a$clades==c,]' >> diversity_richness70.r
echo 'b$V4=(b$count/sum(b$count))*(-log(b$count/sum(b$count)))' >> diversity_richness70.r
echo 'sd=sum(b$V4)' >> diversity_richness70.r
echo 'd=data.frame(c,sd)' >> diversity_richness70.r
echo 'finaldf=rbind(finaldf,d)' >> diversity_richness70.r
echo '}' >> diversity_richness70.r
echo 'jpeg("shannon_diversity_morethanseventy.jpeg",width=23.2,height=13,units="in",res=300)' >> diversity_richness70.r
echo 'barplot(finaldf$sd,names=finaldf$c,ylab="Shannon diversity",ylim=c(0,max(finaldf$sd)),xlab="Clades",col="lightcyan",main="Shannon diversity of >70% LPSs across different clades")' >> diversity_richness70.r
echo 'dev.off()' >> diversity_richness70.r
####for Simpson diversity####
echo 'simpsondf <- data.frame(clades=character(), simd=numeric())' >> diversity_richness70.r
echo 'for(c in unique(a$clades)){' >> diversity_richness70.r
echo 's=a[a$clades==c,]' >> diversity_richness70.r
echo 's$V4=s$count*(s$count - 1)' >> diversity_richness70.r
echo 'simpd=round(1-(sum(s$V4)/(sum(s$count)*(sum(s$count)-1))),4)' >> diversity_richness70.r
echo 'simp=data.frame(c,simpd)' >> diversity_richness70.r
echo 'simpsondf=rbind(simpsondf,simp)}' >> diversity_richness70.r
echo 'jpeg("simpson_diversity_morethanseventy.jpeg",width=23.2,height=13,units="in",res=300)' >> diversity_richness70.r
echo 'barplot(simpsondf$simpd,names=simpsondf$c,ylab="Simpson diversity",ylim=c(0,1),xlab="Clades",col="lightcyan",main="Simpson diversity of >70% LPSs across different clades")' >> diversity_richness70.r
echo 'dev.off()' >> diversity_richness70.r
####Richness and slope of richness
echo 'unnorrich <- data.frame(clades=character(), richness=numeric())' >> diversity_richness70.r
echo 'for(c in unique(a$clades)){' >> diversity_richness70.r
echo 'b=a[a$clades==c,]' >> diversity_richness70.r
echo 'numrep=dim(b)[1]' >> diversity_richness70.r
echo 'ric=data.frame(c,numrep)' >> diversity_richness70.r
echo 'unnorrich=rbind(unnorrich,ric)}' >> diversity_richness70.r
echo 'jpeg("unnormalized_repeat_richness_morethanseventy.jpeg",width=23.2,height=13,units="in",res=300)' >> diversity_richness70.r
echo 'barplot(unnorrich$numrep,names=unnorrich$c,ylab="Repeat richness",xlab="Clades",ylim=c(0,max(unnorrich$numrep)+2),col="#69b3a2",main="Repeats richness for >70% LPSs across different clades")' >> diversity_richness70.r
echo 'dev.off()' >> diversity_richness70.r
##normalizing the richness
echo 'numlps=read.table("numberoflps.txt",header=F)' >> diversity_richness70.r
echo 'numlps$V1=str_to_title(numlps$V1)' >> diversity_richness70.r
echo 'norrich <- data.frame(clades=character(), richness=numeric())' >> diversity_richness70.r
echo 'for(c in unique(a$clades)){' >> diversity_richness70.r
echo 'rich=unnorrich[unnorrich$c==c,2]' >> diversity_richness70.r
echo 'numrep=numlps[numlps$V1==c,2]' >> diversity_richness70.r
echo 'slope=log(rich)/log(numrep)' >> diversity_richness70.r
echo 'rich=data.frame(c,slope)' >> diversity_richness70.r
echo 'norrich=rbind(norrich,rich)}' >> diversity_richness70.r
echo 'jpeg("normalized_repeat_richness_morethanseventy.jpeg",width=23.2,height=13,units="in",res=300)' >> diversity_richness70.r
echo 'barplot(norrich$slope,names=norrich$c,ylab="Repeat richness",xlab="Clades",ylim=c(0,max(norrich$slope)+0.2),col="#69b3a2",main="Normalized repeats richness for >70% LPSs across different clades")' >> diversity_richness70.r
echo 'dev.off()' >> diversity_richness70.r
##slope of richness
###log-transformed richness and number of sequences
echo 'norrich <- data.frame(clades=character(), richness=numeric(),numrepeats=numeric())' >> diversity_richness70.r
echo 'for(c in unique(a$clades)){' >> diversity_richness70.r
echo 'rich=unnorrich[unnorrich$c==c,2]' >> diversity_richness70.r
echo 'numrep=numlps[numlps$V1==c,2]' >> diversity_richness70.r
echo 'rich=data.frame(c,rich,numrep)' >> diversity_richness70.r
echo 'norrich=rbind(norrich,rich)}' >> diversity_richness70.r
echo 'jpeg("log_transformed_richness_morethanseventy.jpeg",width=23.2,height=13,units="in",res=300)' >> diversity_richness70.r
echo 'b=lm(formula = log(norrich$rich) ~ log(norrich$numrep), data = norrich)' >> diversity_richness70.r
echo 'plot(log(norrich$numrep),log(norrich$rich),ylab="log(richness)",xlab="log(number of sequences)",col="blue3",pch=20,main="log transformed richness of >70% repeats and number of sequences across different clades")' >> diversity_richness70.r
echo 'abline(b,col="red",lty=2)' >> diversity_richness70.r
echo 'dev.off()' >> diversity_richness70.r
Rscript diversity_richness70.r
####R script end here####


#####################For exact stretch now###################
cut -f1 exactstretch_combined|sort|uniq -c|awk '{print $2,$1}' OFS='\t' > numberoflps_exact.txt
###############################
##calculating diversity cladewise for exact repeats
awk '{print $1,$11}' OFS='\t' exactstretch_combined > exactstretch_diversity
for clades in afrotheria amphibia artiodactyla aves carnivore chiroptera lagomorpha marsupials perissodactyla primates rodents squamata testudines
do
grep "$clades" exactstretch_diversity|sort -k2|uniq -c|awk '{print $2,$3,$1}' OFS='\t'|sort -k3nr >> temp.txt
done
mv temp.txt exactstretch_diversity
################Plotting in R#############
##R
echo 'a=read.table("exactstretch_diversity",header=F)' > diversity_richness_exact.r
echo 'colnames(a)=c("clades","repeat","count")' >> diversity_richness_exact.r
echo 'library(stringr)' >> diversity_richness_exact.r
echo 'a$clades=str_to_title(a$clades)' >> diversity_richness_exact.r
####for Shannon diversity####
echo 'finaldf <- data.frame(clades=character(), shd=numeric())' >> diversity_richness_exact.r
echo 'for(c in unique(a$clades)){' >> diversity_richness_exact.r
echo 'b=a[a$clades==c,]' >> diversity_richness_exact.r
echo 'b$V4=(b$count/sum(b$count))*(-log(b$count/sum(b$count)))' >> diversity_richness_exact.r
echo 'sd=sum(b$V4)' >> diversity_richness_exact.r
echo 'd=data.frame(c,sd)' >> diversity_richness_exact.r
echo 'finaldf=rbind(finaldf,d)' >> diversity_richness_exact.r
echo '}' >> diversity_richness_exact.r
echo 'jpeg("shannon_diversity_exactstretch.jpeg",width=23.2,height=13,units="in",res=300)' >> diversity_richness_exact.r
echo 'barplot(finaldf$sd,names=finaldf$c,ylab="Shannon diversity",ylim=c(0,max(finaldf$sd)),xlab="Clades",col="lightcyan",main="Shannon diversity of exact LPSs across different clades")' >> diversity_richness_exact.r
echo 'dev.off()' >> diversity_richness_exact.r
####for Simpson diversity####
echo 'simpsondf <- data.frame(clades=character(), simd=numeric())' >> diversity_richness_exact.r
echo 'for(c in unique(a$clades)){' >> diversity_richness_exact.r
echo 's=a[a$clades==c,]' >> diversity_richness_exact.r
echo 's$V4=s$count*(s$count - 1)' >> diversity_richness_exact.r
echo 'simpd=round(1-(sum(s$V4)/(sum(s$count)*(sum(s$count)-1))),4)' >> diversity_richness_exact.r
echo 'simp=data.frame(c,simpd)' >> diversity_richness_exact.r
echo 'simpsondf=rbind(simpsondf,simp)}' >> diversity_richness_exact.r
echo 'jpeg("simpson_diversity_exactstretch.jpeg",width=23.2,height=13,units="in",res=300)' >> diversity_richness_exact.r
echo 'barplot(simpsondf$simpd,names=simpsondf$c,ylab="Simpson diversity",ylim=c(0,1),xlab="Clades",col="lightcyan",main="Simpson diversity of exact LPSs across different clades")' >> diversity_richness_exact.r
echo 'dev.off()' >> diversity_richness_exact.r
####Richness and slope of richness
echo 'unnorrich <- data.frame(clades=character(), richness=numeric())' >> diversity_richness_exact.r
echo 'for(c in unique(a$clades)){' >> diversity_richness_exact.r
echo 'b=a[a$clades==c,]' >> diversity_richness_exact.r
echo 'numrep=dim(b)[1]' >> diversity_richness_exact.r
echo 'ric=data.frame(c,numrep)' >> diversity_richness_exact.r
echo 'unnorrich=rbind(unnorrich,ric)}' >> diversity_richness_exact.r
echo 'jpeg("unnormalized_repeat_richness_exactstretch.jpeg",width=23.2,height=13,units="in",res=300)' >> diversity_richness_exact.r
echo 'barplot(unnorrich$numrep,names=unnorrich$c,ylab="Repeat richness",xlab="Clades",ylim=c(0,max(unnorrich$numrep)+2),col="#69b3a2",main="Repeats richness for exact LPSs across different clades")' >> diversity_richness_exact.r
echo 'dev.off()' >> diversity_richness_exact.r
##normalizing the richness
echo 'numlps=read.table("numberoflps_exact.txt",header=F)' >> diversity_richness_exact.r
echo 'numlps$V1=str_to_title(numlps$V1)' >> diversity_richness_exact.r
echo 'norrich <- data.frame(clades=character(), richness=numeric())' >> diversity_richness_exact.r
echo 'for(c in unique(a$clades)){' >> diversity_richness_exact.r
echo 'rich=unnorrich[unnorrich$c==c,2]' >> diversity_richness_exact.r
echo 'numrep=numlps[numlps$V1==c,2]' >> diversity_richness_exact.r
echo 'slope=log(rich)/log(numrep)' >> diversity_richness_exact.r
echo 'rich=data.frame(c,slope)' >> diversity_richness_exact.r
echo 'norrich=rbind(norrich,rich)}' >> diversity_richness_exact.r
echo 'jpeg("normalized_repeat_richness_exactstretch.jpeg",width=23.2,height=13,units="in",res=300)' >> diversity_richness_exact.r
echo 'barplot(norrich$slope,names=norrich$c,ylab="Repeat richness",xlab="Clades",ylim=c(0,max(norrich$slope)+0.2),col="#69b3a2",main="Normalized repeats richness for exact LPSs across different clades")' >> diversity_richness_exact.r
echo 'dev.off()' >> diversity_richness_exact.r
##slope of richness
###log-transformed richness and number of sequences
echo 'norrich <- data.frame(clades=character(), richness=numeric(),numrepeats=numeric())' >> diversity_richness_exact.r
echo 'for(c in unique(a$clades)){' >> diversity_richness_exact.r
echo 'rich=unnorrich[unnorrich$c==c,2]' >> diversity_richness_exact.r
echo 'numrep=numlps[numlps$V1==c,2]' >> diversity_richness_exact.r
echo 'rich=data.frame(c,rich,numrep)' >> diversity_richness_exact.r
echo 'norrich=rbind(norrich,rich)}' >> diversity_richness_exact.r
echo 'jpeg("log_transformed_richness_exactstretch.jpeg",width=23.2,height=13,units="in",res=300)' >> diversity_richness_exact.r
echo 'b=lm(formula = log(norrich$rich) ~ log(norrich$numrep), data = norrich)' >> diversity_richness_exact.r
echo 'plot(log(norrich$numrep),log(norrich$rich),ylab="log(richness)",xlab="log(number of sequences)",col="blue3",pch=20,main="log transformed richness of exact repeats and number of sequences across different clades")' >> diversity_richness_exact.r
echo 'abline(b,col="red",lty=2)' >> diversity_richness_exact.r
echo 'dev.off()' >> diversity_richness_exact.r
Rscript diversity_richness_exact.r


######################
#################### For morethanseventy repeats where species is also included ###########################
## making boxplots for richness and density
###### making the diversity input file
awk '{print $1,$3,$11}' OFS='\t' morethanseventy_combined|awk '{ gsub(/.*.[0-9]_/,"", $2); print }' OFS="\t" > morethanseventy_diversity_specieswise
for clades in afrotheria amphibia artiodactyla aves carnivore chiroptera lagomorpha marsupials perissodactyla primates rodents squamata testudines
do
grep "$clades" morethanseventy_diversity_specieswise |sort -k2|uniq -c|awk '{print $2,$3,$4,$1}' OFS='\t'|sort -k3nr >> temp.txt
done
mv temp.txt morethanseventy_diversity_specieswise
######### for diversity ###################
### making the diversity input file ###
##R
echo 'a=read.table("morethanseventy_diversity_specieswise",header=F)' > pairwise_diversity_specieswise70.r
echo 'colnames(a)=c("clades","species","repeat","count")' >> pairwise_diversity_specieswise70.r
echo 'library(stringr)' >> pairwise_diversity_specieswise70.r
echo 'a$clades=str_to_title(a$clades)' >> pairwise_diversity_specieswise70.r
echo '#install.packages("qpcR")' >> pairwise_diversity_specieswise70.r
echo 'library("qpcR")' >> pairwise_diversity_specieswise70.r
echo '####for Shannon diversity####' >> pairwise_diversity_specieswise70.r
echo 'finaldfr <- data.frame(first=c(1:10))' >> pairwise_diversity_specieswise70.r
echo 'for(c in unique(a$clades)){' >> pairwise_diversity_specieswise70.r
echo 'finaldf <- data.frame(clades=character(), shd=numeric())' >> pairwise_diversity_specieswise70.r
echo 'b1=a[a$clades==c,]' >> pairwise_diversity_specieswise70.r
echo 'for(sp in unique(b1$species)){' >> pairwise_diversity_specieswise70.r
echo 'b=b1[b1$species==sp,]' >> pairwise_diversity_specieswise70.r
echo 'b$V4=(b$count/sum(b$count))*(-log(b$count/sum(b$count)))' >> pairwise_diversity_specieswise70.r
echo 'sd1=sum(b$V4)' >> pairwise_diversity_specieswise70.r
echo 'sd=data.frame(sd1)' >> pairwise_diversity_specieswise70.r
echo 'finaldf=rbind(finaldf,sd)' >> pairwise_diversity_specieswise70.r
echo '}' >> pairwise_diversity_specieswise70.r
echo 'colnames(finaldf)=c' >> pairwise_diversity_specieswise70.r
echo 'finaldfr=qpcR:::cbind.na(finaldfr,finaldf)' >> pairwise_diversity_specieswise70.r
echo '}' >> pairwise_diversity_specieswise70.r
echo 'd=finaldfr[,-1]' >> pairwise_diversity_specieswise70.r
echo 'library(stringr)' >> pairwise_diversity_specieswise70.r
echo 'str_to_title(colnames(d))->colnames(d)' >> pairwise_diversity_specieswise70.r
echo 'species<-colnames(d)' >> pairwise_diversity_specieswise70.r
echo 'len<-length(species)' >> pairwise_diversity_specieswise70.r
echo 'totalcol=dim(d)[2]' >> pairwise_diversity_specieswise70.r
echo 'ymin=c(min(d[,1],na.rm=T))' >> pairwise_diversity_specieswise70.r
echo 'for ( i in 2:length(d)){' >> pairwise_diversity_specieswise70.r
echo 'ymin[i]<-c(min(d[,i],na.rm=T))}' >> pairwise_diversity_specieswise70.r
echo 'maxquan=0' >> pairwise_diversity_specieswise70.r
echo 'P<- data.frame(id = character(0), p.value = numeric(0))' >> pairwise_diversity_specieswise70.r
echo 'for ( i in seq(1,(totalcol-1),1)){' >> pairwise_diversity_specieswise70.r
echo 'for ( j in seq(i+1,totalcol,1)){' >> pairwise_diversity_specieswise70.r
echo 'if (species[i] == species[j]){' >> pairwise_diversity_specieswise70.r
echo 'next' >> pairwise_diversity_specieswise70.r
echo '}' >> pairwise_diversity_specieswise70.r
echo 'else' >> pairwise_diversity_specieswise70.r
echo '{' >> pairwise_diversity_specieswise70.r
echo 'print(paste(species[i],species[j],sep=" "))' >> pairwise_diversity_specieswise70.r
echo 'fg1<-d[,i]' >> pairwise_diversity_specieswise70.r
echo 'bg1<-d[,j]' >> pairwise_diversity_specieswise70.r
echo 'if (maxquan < quantile(fg1,0.75,na.rm=T)){maxquan=quantile(fg1,0.75,na.rm=T)}' >> pairwise_diversity_specieswise70.r
echo 'wilcox.test(fg1,bg1,alternative = "two.sided")->K' >> pairwise_diversity_specieswise70.r
echo 'P1<-data.frame(names=paste(species[i],"vs",species[j],sep="_"),p.value=K$p.value)' >> pairwise_diversity_specieswise70.r
echo 'seg2=paste(species[i],"vs",species[j],sep="_")' >> pairwise_diversity_specieswise70.r
echo 'P<-rbind(P,P1) }' >> pairwise_diversity_specieswise70.r
echo '}' >> pairwise_diversity_specieswise70.r
echo '}' >> pairwise_diversity_specieswise70.r
echo 'options(scipen=0,"digits"=2)' >> pairwise_diversity_specieswise70.r
echo 'P$qvalue<-p.adjust(P$p.value, method = "BH")' >> pairwise_diversity_specieswise70.r
echo 'P$V1<-format(P$qvalue,scientific=T)' >> pairwise_diversity_specieswise70.r
echo 'jpeg("Pairwise_boxplot_specieswise_Shannon_diversity_morethanseventy.jpeg",height=2500,width=3600)' >> pairwise_diversity_specieswise70.r
echo 'par(mai=c(5,5,2,1),mgp=c(18,0.5,0.1))' >> pairwise_diversity_specieswise70.r
echo 'boxplot(d,outline=F,main="Shannon diversity of >70% LPSs across different clades", xlab="Clades",ylab="Shannon diversity", col="Lightsalmon", border="black", las=2, cex=4, cex.main=4, cex.axis=3, cex.lab=4, lwd=2, ylim=c(min(ymin),as.numeric(maxquan)+1+(2*(length(P[P$qvalue<0.05,3])))))' >> pairwise_diversity_specieswise70.r
echo 'start=maxquan' >> pairwise_diversity_specieswise70.r
echo 'for ( i in seq(1,len,1)){' >> pairwise_diversity_specieswise70.r
echo 'for (j in seq(2,len,1)){' >> pairwise_diversity_specieswise70.r
echo 'tryCatch({' >> pairwise_diversity_specieswise70.r
echo 'if (species[i] == species[j]){' >> pairwise_diversity_specieswise70.r
echo 'next' >> pairwise_diversity_specieswise70.r
echo '}' >> pairwise_diversity_specieswise70.r
echo 'else' >> pairwise_diversity_specieswise70.r
echo '{' >> pairwise_diversity_specieswise70.r
echo 'fg=as.character(species[i])' >> pairwise_diversity_specieswise70.r
echo 'bg=as.character(species[j])' >> pairwise_diversity_specieswise70.r
echo 'seg=paste(fg,"vs",bg,sep="_")' >> pairwise_diversity_specieswise70.r
echo 'P[P$names==seg,3]-> qval' >> pairwise_diversity_specieswise70.r
echo 'P[P$names==seg,4] -> qvalchar' >> pairwise_diversity_specieswise70.r
echo 'if(qval > 0.05){' >> pairwise_diversity_specieswise70.r
echo 'next' >> pairwise_diversity_specieswise70.r
echo '}' >> pairwise_diversity_specieswise70.r
echo 'else' >> pairwise_diversity_specieswise70.r
echo '{' >> pairwise_diversity_specieswise70.r
echo 'segments(i,start,j,start,lty=3,col="purple",lwd=2)' >> pairwise_diversity_specieswise70.r
echo 'text(((i+j)/2),start+1,labels=qvalchar,cex=2)' >> pairwise_diversity_specieswise70.r
echo '}' >> pairwise_diversity_specieswise70.r
echo '}' >> pairwise_diversity_specieswise70.r
echo 'start = start+2' >> pairwise_diversity_specieswise70.r
echo '}' >> pairwise_diversity_specieswise70.r
echo ', error=function(e){})}' >> pairwise_diversity_specieswise70.r
echo '}' >> pairwise_diversity_specieswise70.r
echo 'dev.off()' >> pairwise_diversity_specieswise70.r
echo 'jpeg("Shannon_diversity_specieswise_boxplot_morethanseventy.jpeg",height=2500,width=3600)' >> pairwise_diversity_specieswise70.r
echo 'par(mai=c(5,5,2,1),mgp=c(18,0.5,0.1))' >> pairwise_diversity_specieswise70.r
echo 'boxplot(d,outline=F,main="Shannon diversity of >70% LPSs across different clades", xlab="Clades",ylab="Shannon diversity", col="Lightsalmon", border="black", las=2, cex=4, cex.main=4, cex.axis=3, cex.lab=4, lwd=2)' >> pairwise_diversity_specieswise70.r
echo 'dev.off()' >> pairwise_diversity_specieswise70.r
echo '####for Simpson diversity####' >> pairwise_diversity_specieswise70.r
echo 'finaldfr<-data.frame(first=c(1:10))' >> pairwise_diversity_specieswise70.r
echo 'for(c in unique(a$clades)){' >> pairwise_diversity_specieswise70.r
echo 'finaldf<-data.frame(clades=character(), shd=numeric())' >> pairwise_diversity_specieswise70.r
echo 's1=a[a$clades==c,]' >> pairwise_diversity_specieswise70.r
echo 'for(sp in unique(s1$species)){' >> pairwise_diversity_specieswise70.r
echo 's=s1[s1$species==sp,]' >> pairwise_diversity_specieswise70.r
echo 's$V4=s$count*(s$count - 1)' >> pairwise_diversity_specieswise70.r
echo 'simpd=round(1-(sum(s$V4)/(sum(s$count)*(sum(s$count)-1))),4)' >> pairwise_diversity_specieswise70.r
echo 'sd=data.frame(simpd)' >> pairwise_diversity_specieswise70.r
echo 'finaldf=rbind(finaldf,sd)' >> pairwise_diversity_specieswise70.r
echo '}' >> pairwise_diversity_specieswise70.r
echo 'colnames(finaldf)=c' >> pairwise_diversity_specieswise70.r
echo 'finaldfr=qpcR:::cbind.na(finaldfr,finaldf)' >> pairwise_diversity_specieswise70.r
echo '}' >> pairwise_diversity_specieswise70.r
echo 'd=finaldfr[,-1]' >> pairwise_diversity_specieswise70.r
echo 'library(stringr)' >> pairwise_diversity_specieswise70.r
echo 'str_to_title(colnames(d))->colnames(d)' >> pairwise_diversity_specieswise70.r
echo 'species<-colnames(d)' >> pairwise_diversity_specieswise70.r
echo 'len<-length(species)' >> pairwise_diversity_specieswise70.r
echo 'totalcol=dim(d)[2]' >> pairwise_diversity_specieswise70.r
echo 'ymin=c(min(d[,1],na.rm=T))' >> pairwise_diversity_specieswise70.r
echo 'for ( i in 2:length(d)){' >> pairwise_diversity_specieswise70.r
echo 'ymin[i]<-c(min(d[,i],na.rm=T))}' >> pairwise_diversity_specieswise70.r
echo 'maxquan=0' >> pairwise_diversity_specieswise70.r
echo 'P<- data.frame(id = character(0), p.value = numeric(0))' >> pairwise_diversity_specieswise70.r
echo 'for ( i in seq(1,(totalcol-1),1)){' >> pairwise_diversity_specieswise70.r
echo 'for ( j in seq(i+1,totalcol,1)){' >> pairwise_diversity_specieswise70.r
echo 'if (species[i] == species[j]){' >> pairwise_diversity_specieswise70.r
echo 'next' >> pairwise_diversity_specieswise70.r
echo '}' >> pairwise_diversity_specieswise70.r
echo 'else' >> pairwise_diversity_specieswise70.r
echo '{' >> pairwise_diversity_specieswise70.r
echo 'print(paste(species[i],species[j],sep=" "))' >> pairwise_diversity_specieswise70.r
echo 'fg1<-d[,i]' >> pairwise_diversity_specieswise70.r
echo 'bg1<-d[,j]' >> pairwise_diversity_specieswise70.r
echo 'if (maxquan < quantile(fg1,0.75,na.rm=T)){maxquan=quantile(fg1,0.75,na.rm=T) }' >> pairwise_diversity_specieswise70.r
echo 'wilcox.test(fg1,bg1,alternative = "two.sided")->K' >> pairwise_diversity_specieswise70.r
echo 'P1<-data.frame(names=paste(species[i],"vs",species[j],sep="_"),p.value=K$p.value)' >> pairwise_diversity_specieswise70.r
echo 'seg2=paste(species[i],"vs",species[j],sep="_")' >> pairwise_diversity_specieswise70.r
echo 'P<-rbind(P,P1) }' >> pairwise_diversity_specieswise70.r
echo '}' >> pairwise_diversity_specieswise70.r
echo '}' >> pairwise_diversity_specieswise70.r
echo 'options(scipen=0,"digits"=2)' >> pairwise_diversity_specieswise70.r
echo 'P$qvalue<-p.adjust(P$p.value, method = "BH")' >> pairwise_diversity_specieswise70.r
echo 'P$V1<-format(P$qvalue,scientific=T)' >> pairwise_diversity_specieswise70.r
echo 'jpeg("Pairwise_boxplot_specieswise_Simpson_diversity_morethanseventy.jpeg",height=2500,width=3600)' >> pairwise_diversity_specieswise70.r
echo 'par(mai=c(5,5,2,1),mgp=c(18,0.5,0.1))' >> pairwise_diversity_specieswise70.r
echo 'boxplot(d,outline=F,main="Simpson diversity of >70% LPSs across different clades", xlab="Clades",ylab="Simpson diversity", col="Lightsalmon", border="black", las=2, cex=4, cex.main=4, cex.axis=3, cex.lab=4, lwd=2, ylim=c(min(ymin),as.numeric(maxquan)+1+(2*(length(P[P$qvalue<0.05,3])))))' >> pairwise_diversity_specieswise70.r
echo 'start=maxquan' >> pairwise_diversity_specieswise70.r
echo 'for ( i in seq(1,len,1)){' >> pairwise_diversity_specieswise70.r
echo 'for (j in seq(2,len,1)){' >> pairwise_diversity_specieswise70.r
echo 'tryCatch({' >> pairwise_diversity_specieswise70.r
echo 'if (species[i] == species[j]){' >> pairwise_diversity_specieswise70.r
echo 'next' >> pairwise_diversity_specieswise70.r
echo '}' >> pairwise_diversity_specieswise70.r
echo 'else' >> pairwise_diversity_specieswise70.r
echo '{' >> pairwise_diversity_specieswise70.r
echo 'fg=as.character(species[i])' >> pairwise_diversity_specieswise70.r
echo 'bg=as.character(species[j])' >> pairwise_diversity_specieswise70.r
echo 'seg=paste(fg,"vs",bg,sep="_")' >> pairwise_diversity_specieswise70.r
echo 'P[P$names==seg,3]-> qval' >> pairwise_diversity_specieswise70.r
echo 'P[P$names==seg,4] -> qvalchar' >> pairwise_diversity_specieswise70.r
echo 'if(qval > 0.05){' >> pairwise_diversity_specieswise70.r
echo 'next' >> pairwise_diversity_specieswise70.r
echo '}' >> pairwise_diversity_specieswise70.r
echo 'else' >> pairwise_diversity_specieswise70.r
echo '{' >> pairwise_diversity_specieswise70.r
echo 'segments(i,start,j,start,lty=3,col="purple",lwd=2)' >> pairwise_diversity_specieswise70.r
echo 'text(((i+j)/2),start+1,labels=qvalchar,cex=2)' >> pairwise_diversity_specieswise70.r
echo '}' >> pairwise_diversity_specieswise70.r
echo '}' >> pairwise_diversity_specieswise70.r
echo 'start = start+2' >> pairwise_diversity_specieswise70.r
echo '}' >> pairwise_diversity_specieswise70.r
echo ', error=function(e){})}' >> pairwise_diversity_specieswise70.r
echo '}' >> pairwise_diversity_specieswise70.r
echo 'dev.off()' >> pairwise_diversity_specieswise70.r
echo 'jpeg("Simpson_diversity_specieswise_boxplot_morethanseventy.jpeg",height=2500,width=3600)' >> pairwise_diversity_specieswise70.r
echo 'par(mai=c(5,5,2,1),mgp=c(18,0.5,0.1))' >> pairwise_diversity_specieswise70.r
echo 'boxplot(d,outline=F,main="Simpson diversity of >70% LPSs across different clades", xlab="Clades",ylab="Simpson diversity", col="Lightsalmon", border="black", las=2, cex=4, cex.main=4, cex.axis=3, cex.lab=4, lwd=2)' >> pairwise_diversity_specieswise70.r
echo 'dev.off()' >> pairwise_diversity_specieswise70.r
Rscript pairwise_diversity_specieswise70.r
############################################ calculating richness cladewise for >70% LPS ###############################
####Richness and slope of richness
##R
echo 'a=read.table("morethanseventy_diversity_specieswise",header=F)' > pairwise_richness_specieswise70.r
echo 'colnames(a)=c("clades","species","repeat","count")' >> pairwise_richness_specieswise70.r
echo 'library(stringr)' >> pairwise_richness_specieswise70.r
echo 'a$clades=str_to_title(a$clades)' >> pairwise_richness_specieswise70.r
echo 'finaldfr<-data.frame(first=c(1:10))' >> pairwise_richness_specieswise70.r
echo 'for(c in unique(a$clades)){' >> pairwise_richness_specieswise70.r
echo 'finaldf<-data.frame(clades=character(), shd=numeric())' >> pairwise_richness_specieswise70.r
echo 's1=a[a$clades==c,]' >> pairwise_richness_specieswise70.r
echo 'for(sp in unique(s1$species)){' >> pairwise_richness_specieswise70.r
echo 's=s1[s1$species==sp,]' >> pairwise_richness_specieswise70.r
echo 'simpd=dim(s)[1]' >> pairwise_richness_specieswise70.r
echo 'sd=data.frame(simpd)' >> pairwise_richness_specieswise70.r
echo 'finaldf=rbind(finaldf,sd)' >> pairwise_richness_specieswise70.r
echo '}' >> pairwise_richness_specieswise70.r
echo 'colnames(finaldf)=c' >> pairwise_richness_specieswise70.r
echo 'finaldfr=qpcR:::cbind.na(finaldfr,finaldf)' >> pairwise_richness_specieswise70.r
echo '}' >> pairwise_richness_specieswise70.r
echo 'd=finaldfr[,-1]' >> pairwise_richness_specieswise70.r
echo 'library(stringr)' >> pairwise_richness_specieswise70.r
echo 'str_to_title(colnames(d))->colnames(d)' >> pairwise_richness_specieswise70.r
echo 'species<-colnames(d)' >> pairwise_richness_specieswise70.r
echo 'len<-length(species)' >> pairwise_richness_specieswise70.r
echo 'totalcol=dim(d)[2]' >> pairwise_richness_specieswise70.r
echo 'ymin=c(min(d[,1],na.rm=T))' >> pairwise_richness_specieswise70.r
echo 'for ( i in 2:length(d)){' >> pairwise_richness_specieswise70.r
echo 'ymin[i]<-c(min(d[,i],na.rm=T))}' >> pairwise_richness_specieswise70.r
echo 'maxquan=0' >> pairwise_richness_specieswise70.r
echo 'P<- data.frame(id = character(0), p.value = numeric(0))' >> pairwise_richness_specieswise70.r
echo 'for ( i in seq(1,(totalcol-1),1)){' >> pairwise_richness_specieswise70.r
echo 'for ( j in seq(i+1,totalcol,1)){' >> pairwise_richness_specieswise70.r
echo 'if (species[i] == species[j]){' >> pairwise_richness_specieswise70.r
echo 'next' >> pairwise_richness_specieswise70.r
echo '}' >> pairwise_richness_specieswise70.r
echo 'else' >> pairwise_richness_specieswise70.r
echo '{' >> pairwise_richness_specieswise70.r
echo 'print(paste(species[i],species[j],sep=" "))' >> pairwise_richness_specieswise70.r
echo 'fg1<-d[,i]' >> pairwise_richness_specieswise70.r
echo 'bg1<-d[,j]' >> pairwise_richness_specieswise70.r
echo 'if (maxquan < quantile(fg1,0.75,na.rm=T)){maxquan=quantile(fg1,0.75,na.rm=T) }' >> pairwise_richness_specieswise70.r
echo 'wilcox.test(fg1,bg1,alternative = "two.sided")->K' >> pairwise_richness_specieswise70.r
echo 'P1<-data.frame(names=paste(species[i],"vs",species[j],sep="_"),p.value=K$p.value)' >> pairwise_richness_specieswise70.r
echo 'seg2=paste(species[i],"vs",species[j],sep="_")' >> pairwise_richness_specieswise70.r
echo 'P<-rbind(P,P1) }' >> pairwise_richness_specieswise70.r
echo '}' >> pairwise_richness_specieswise70.r
echo '}' >> pairwise_richness_specieswise70.r
echo 'options(scipen=0,"digits"=2)' >> pairwise_richness_specieswise70.r
echo 'P$qvalue<-p.adjust(P$p.value, method = "BH")' >> pairwise_richness_specieswise70.r
echo 'P$V1<-format(P$qvalue,scientific=T)' >> pairwise_richness_specieswise70.r
echo 'jpeg("Pairwise_boxplot_specieswise_richness_morethanseventy.jpeg",height=2500,width=3600)' >> pairwise_richness_specieswise70.r
echo 'par(mai=c(5,5,2,1),mgp=c(18,0.5,0.1))' >> pairwise_richness_specieswise70.r
echo 'boxplot(d,outline=F,main="Richness of >70% LPSs across clades", xlab="Clades",ylab="Richness of LPSs", col="Lightsalmon", border="black", las=2, cex=4, cex.main=4, cex.axis=3, cex.lab=4, lwd=2, ylim=c(min(ymin),as.numeric(maxquan)+1+(2*(length(P[P$qvalue<0.05,3])))))' >> pairwise_richness_specieswise70.r
echo 'start=maxquan' >> pairwise_richness_specieswise70.r
echo 'for ( i in seq(1,len,1)){' >> pairwise_richness_specieswise70.r
echo 'for (j in seq(2,len,1)){' >> pairwise_richness_specieswise70.r
echo 'tryCatch({' >> pairwise_richness_specieswise70.r
echo 'if (species[i] == species[j]){' >> pairwise_richness_specieswise70.r
echo 'next' >> pairwise_richness_specieswise70.r
echo '}' >> pairwise_richness_specieswise70.r
echo 'else' >> pairwise_richness_specieswise70.r
echo '{' >> pairwise_richness_specieswise70.r
echo 'fg=as.character(species[i])' >> pairwise_richness_specieswise70.r
echo 'bg=as.character(species[j])' >> pairwise_richness_specieswise70.r
echo 'seg=paste(fg,"vs",bg,sep="_")' >> pairwise_richness_specieswise70.r
echo 'P[P$names==seg,3]-> qval' >> pairwise_richness_specieswise70.r
echo 'P[P$names==seg,4] -> qvalchar' >> pairwise_richness_specieswise70.r
echo 'if(qval > 0.05){' >> pairwise_richness_specieswise70.r
echo 'next' >> pairwise_richness_specieswise70.r
echo '}' >> pairwise_richness_specieswise70.r
echo 'else' >> pairwise_richness_specieswise70.r
echo '{' >> pairwise_richness_specieswise70.r
echo 'segments(i,start,j,start,lty=3,col="purple",lwd=2)' >> pairwise_richness_specieswise70.r
echo 'text(((i+j)/2),start+1,labels=qvalchar,cex=2)' >> pairwise_richness_specieswise70.r
echo '}' >> pairwise_richness_specieswise70.r
echo '}' >> pairwise_richness_specieswise70.r
echo 'start = start+2' >> pairwise_richness_specieswise70.r
echo '}' >> pairwise_richness_specieswise70.r
echo ', error=function(e){})}' >> pairwise_richness_specieswise70.r
echo '}' >> pairwise_richness_specieswise70.r
echo 'dev.off()' >> pairwise_richness_specieswise70.r
echo 'jpeg("Boxplot_richness_specieswise_morethanseventy.jpeg",height=2500,width=3600)' >> pairwise_richness_specieswise70.r
echo 'par(mai=c(5,5,2,1),mgp=c(18,0.5,0.1))' >> pairwise_richness_specieswise70.r
echo 'boxplot(d,outline=F,main="Richness of >70% LPSs across clades", xlab="Clades",ylab="Richness of LPSs", col="Lightsalmon", border="black", las=2, cex=4, cex.main=4, cex.axis=3, cex.lab=4, lwd=2)' >> pairwise_richness_specieswise70.r
echo 'dev.off()' >> pairwise_richness_specieswise70.r
echo '####### for slope (normalized) #########' >> pairwise_richness_specieswise70.r
echo 'finaldfr <- data.frame(first=c(1:10))' >> pairwise_richness_specieswise70.r
echo 'for(c in unique(a$clades)){' >> pairwise_richness_specieswise70.r
echo 'finaldf <- data.frame(clades=character(), shd=numeric())' >> pairwise_richness_specieswise70.r
echo 'b1=a[a$clades==c,]' >> pairwise_richness_specieswise70.r
echo 'for(sp in unique(b1$species)){' >> pairwise_richness_specieswise70.r
echo 'b=b1[b1$species==sp,]' >> pairwise_richness_specieswise70.r
echo 'rich1=length(a[a$species==sp,4])' >> pairwise_richness_specieswise70.r
echo 'numrep=sum(a[a$species==sp,4])' >> pairwise_richness_specieswise70.r
echo 'slope=log(rich1)/log(numrep)' >> pairwise_richness_specieswise70.r
echo 'rich=data.frame(slope)' >> pairwise_richness_specieswise70.r
echo 'finaldf=rbind(finaldf,rich)' >> pairwise_richness_specieswise70.r
echo '}' >> pairwise_richness_specieswise70.r
echo 'colnames(finaldf)=c' >> pairwise_richness_specieswise70.r
echo 'finaldfr=qpcR:::cbind.na(finaldfr,finaldf)' >> pairwise_richness_specieswise70.r
echo '}' >> pairwise_richness_specieswise70.r
echo 'd=finaldfr[,-1]' >> pairwise_richness_specieswise70.r
echo 'library(stringr)' >> pairwise_richness_specieswise70.r
echo 'str_to_title(colnames(d))->colnames(d)' >> pairwise_richness_specieswise70.r
echo 'species<-colnames(d)' >> pairwise_richness_specieswise70.r
echo 'len<-length(species)' >> pairwise_richness_specieswise70.r
echo 'totalcol=dim(d)[2]' >> pairwise_richness_specieswise70.r
echo 'ymin=c(min(d[,1],na.rm=T))' >> pairwise_richness_specieswise70.r
echo 'for ( i in 2:length(d)){' >> pairwise_richness_specieswise70.r
echo 'ymin[i]<-c(min(d[,i],na.rm=T))}' >> pairwise_richness_specieswise70.r
echo 'maxquan=0' >> pairwise_richness_specieswise70.r
echo 'P<- data.frame(id = character(0), p.value = numeric(0))' >> pairwise_richness_specieswise70.r
echo 'for ( i in seq(1,(totalcol-1),1)){' >> pairwise_richness_specieswise70.r
echo 'for ( j in seq(i+1,totalcol,1)){' >> pairwise_richness_specieswise70.r
echo 'if (species[i] == species[j]){' >> pairwise_richness_specieswise70.r
echo 'next' >> pairwise_richness_specieswise70.r
echo '}' >> pairwise_richness_specieswise70.r
echo 'else' >> pairwise_richness_specieswise70.r
echo '{' >> pairwise_richness_specieswise70.r
echo 'print(paste(species[i],species[j],sep=" "))' >> pairwise_richness_specieswise70.r
echo 'fg1<-d[,i]' >> pairwise_richness_specieswise70.r
echo 'bg1<-d[,j]' >> pairwise_richness_specieswise70.r
echo 'if (maxquan < quantile(fg1,0.75,na.rm=T)){maxquan=quantile(fg1,0.75,na.rm=T) }' >> pairwise_richness_specieswise70.r
echo 'wilcox.test(fg1,bg1,alternative = "two.sided")->K' >> pairwise_richness_specieswise70.r
echo 'P1<-data.frame(names=paste(species[i],"vs",species[j],sep="_"),p.value=K$p.value)' >> pairwise_richness_specieswise70.r
echo 'seg2=paste(species[i],"vs",species[j],sep="_")' >> pairwise_richness_specieswise70.r
echo 'P<-rbind(P,P1) }' >> pairwise_richness_specieswise70.r
echo '}' >> pairwise_richness_specieswise70.r
echo '}' >> pairwise_richness_specieswise70.r
echo 'options(scipen=0,"digits"=2)' >> pairwise_richness_specieswise70.r
echo 'P$qvalue<-p.adjust(P$p.value, method = "BH")' >> pairwise_richness_specieswise70.r
echo 'P$V1<-format(P$qvalue,scientific=T)' >> pairwise_richness_specieswise70.r
echo 'jpeg("Pairwise_boxplot_specieswise_richness_slope_morethanseventy.jpeg",height=2500,width=3600)' >> pairwise_richness_specieswise70.r
echo 'par(mai=c(5,5,2,1),mgp=c(18,0.5,0.1))' >> pairwise_richness_specieswise70.r
echo 'boxplot(d,outline=F,main="Slope of richness of >70% LPSs across clades", xlab="Clades",ylab="Slope of richness of LPSs", col="Lightsalmon", border="black", las=2, cex=4, cex.main=4, cex.axis=3, cex.lab=4, lwd=2, ylim=c(min(ymin),as.numeric(maxquan)+1+(2*(length(P[P$qvalue<0.05,3])))))' >> pairwise_richness_specieswise70.r
echo 'start=maxquan' >> pairwise_richness_specieswise70.r
echo 'for ( i in seq(1,len,1)){' >> pairwise_richness_specieswise70.r
echo 'for (j in seq(2,len,1)){' >> pairwise_richness_specieswise70.r
echo 'tryCatch({' >> pairwise_richness_specieswise70.r
echo 'if (species[i] == species[j]){' >> pairwise_richness_specieswise70.r
echo 'next' >> pairwise_richness_specieswise70.r
echo '}' >> pairwise_richness_specieswise70.r
echo 'else' >> pairwise_richness_specieswise70.r
echo '{' >> pairwise_richness_specieswise70.r
echo 'fg=as.character(species[i])' >> pairwise_richness_specieswise70.r
echo 'bg=as.character(species[j])' >> pairwise_richness_specieswise70.r
echo 'seg=paste(fg,"vs",bg,sep="_")' >> pairwise_richness_specieswise70.r
echo 'P[P$names==seg,3]-> qval' >> pairwise_richness_specieswise70.r
echo 'P[P$names==seg,4] -> qvalchar' >> pairwise_richness_specieswise70.r
echo 'if(qval > 0.05){' >> pairwise_richness_specieswise70.r
echo 'next' >> pairwise_richness_specieswise70.r
echo '}' >> pairwise_richness_specieswise70.r
echo 'else' >> pairwise_richness_specieswise70.r
echo '{' >> pairwise_richness_specieswise70.r
echo 'segments(i,start,j,start,lty=3,col="purple",lwd=2)' >> pairwise_richness_specieswise70.r
echo 'text(((i+j)/2),start+1,labels=qvalchar,cex=2)' >> pairwise_richness_specieswise70.r
echo '}' >> pairwise_richness_specieswise70.r
echo '}' >> pairwise_richness_specieswise70.r
echo 'start = start+2' >> pairwise_richness_specieswise70.r
echo '}' >> pairwise_richness_specieswise70.r
echo ', error=function(e){})}' >> pairwise_richness_specieswise70.r
echo '}' >> pairwise_richness_specieswise70.r
echo 'dev.off()' >> pairwise_richness_specieswise70.r
echo 'jpeg("Boxplot_specieswise_richness_slope_morethanseventy.jpeg",height=2500,width=3600)' >> pairwise_richness_specieswise70.r
echo 'par(mai=c(5,5,2,1),mgp=c(18,0.5,0.1))' >> pairwise_richness_specieswise70.r
echo 'boxplot(d,outline=F,main="Slope of richness of >70% LPSs across clades", xlab="Clades",ylab="Slope of richness of LPSs", col="Lightsalmon", border="black", las=2, cex=4, cex.main=4, cex.axis=3, cex.lab=4, lwd=2)' >> pairwise_richness_specieswise70.r
echo 'dev.off()' >> pairwise_richness_specieswise70.r
echo '########### slope of richness ###log-transformed richness and number of sequences #############' >> pairwise_richness_specieswise70.r
echo 'finaldf <- data.frame(clades=character(), species=character(), rich=numeric(), numrep=numeric())' >> pairwise_richness_specieswise70.r
echo 'for(c in unique(a$clades)){' >> pairwise_richness_specieswise70.r
echo 'b1=a[a$clades==c,]' >> pairwise_richness_specieswise70.r
echo 'for(sp in unique(b1$species)){' >> pairwise_richness_specieswise70.r
echo 'b=b1[b1$species==sp,]' >> pairwise_richness_specieswise70.r
echo 'rich1=length(a[a$species==sp,4])' >> pairwise_richness_specieswise70.r
echo 'numrep=sum(a[a$species==sp,4])' >> pairwise_richness_specieswise70.r
echo 'slope=log(rich1)/log(numrep)' >> pairwise_richness_specieswise70.r
echo 'rich=data.frame(c,sp,rich1,numrep)' >> pairwise_richness_specieswise70.r
echo 'finaldf=rbind(finaldf,rich)' >> pairwise_richness_specieswise70.r
echo '}' >> pairwise_richness_specieswise70.r
echo '}' >> pairwise_richness_specieswise70.r
echo 'col13rainbow= c("#771155", "#AA4488", "#CC99BB", "#114477", "#4477AA", "#77AADD", "#117777", "#44AAAA", "#77CCCC", "#117744", "#44AA77", "#88CCAA", "#777711")' >> pairwise_richness_specieswise70.r
echo 'groups=as.character(unique(finaldf$c))' >> pairwise_richness_specieswise70.r
echo '## all clades in one-log plots' >> pairwise_richness_specieswise70.r
echo 'clade=groups[1]' >> pairwise_richness_specieswise70.r
echo 'jpeg("log_transformed_richness_specieswise_morethanseventy.jpeg",width=23.2,height=13,units="in",res=300)' >> pairwise_richness_specieswise70.r
echo 'b=lm(formula = log(finaldf$rich[finaldf$c==clade]) ~ log(finaldf$numrep[finaldf$c==clade]), data = finaldf)' >> pairwise_richness_specieswise70.r
echo 'plot(log(finaldf$numrep[finaldf$c==clade]),log(finaldf$rich[finaldf$c==clade]),ylab="log(richness)",xlab="log(number of sequences)",col=col13rainbow[1],pch=20,main="log transformed richness of >70% repeats and number of sequences across different clades")' >> pairwise_richness_specieswise70.r
echo 'abline(b,col=col13rainbow[1],lty=2)' >> pairwise_richness_specieswise70.r
echo 'for(clade in 2:length(groups)){' >> pairwise_richness_specieswise70.r
echo 'b=lm(formula = log(finaldf$rich[finaldf$c==groups[clade]]) ~ log(finaldf$numrep[finaldf$c==groups[clade]]), data = finaldf)' >> pairwise_richness_specieswise70.r
echo 'points(log(finaldf$numrep[finaldf$c==groups[clade]]),log(finaldf$rich[finaldf$c==groups[clade]]),ylab="log(richness)",xlab="log(number of sequences)",col=col13rainbow[clade],pch=20,main="log transformed richness of >70% repeats and number of sequences across different clades",cex=2)' >> pairwise_richness_specieswise70.r
echo 'abline(b,col=col13rainbow[clade],lty=2,lwd=2)' >> pairwise_richness_specieswise70.r
echo '}' >> pairwise_richness_specieswise70.r
echo 'legend("topleft",legend=groups,fill=col13rainbow,col=col13rainbow,bty="n")' >> pairwise_richness_specieswise70.r
echo 'dev.off()' >> pairwise_richness_specieswise70.r
echo '## cladewise log plots' >> pairwise_richness_specieswise70.r
echo 'for(clade in 1:length(groups)){' >> pairwise_richness_specieswise70.r
echo 'jpeg(paste(groups[clade],"_log_transformed_richness_specieswise_morethanseventy.jpeg",sep=""),width=23.2,height=13,units="in",res=300)' >> pairwise_richness_specieswise70.r
echo 'b=lm(formula = log(finaldf$rich[finaldf$c==groups[clade]]) ~ log(finaldf$numrep[finaldf$c==groups[clade]]), data = finaldf)' >> pairwise_richness_specieswise70.r
echo 'plot(log(finaldf$numrep[finaldf$c==groups[clade]]),log(finaldf$rich[finaldf$c==groups[clade]]),ylab="log(richness)",xlab="log(number of sequences)",col=col13rainbow[clade],pch=20,main=paste("log transformed richness of >70% repeats and number of sequences in ", groups[clade],sep=""),cex=2)' >> pairwise_richness_specieswise70.r
echo 'abline(b,col=col13rainbow[clade],lty=2,lwd=2)' >> pairwise_richness_specieswise70.r
echo 'legend("topleft",legend=groups[clade],fill=col13rainbow,col=col13rainbow,bty="n")' >> pairwise_richness_specieswise70.r
echo 'dev.off()' >> pairwise_richness_specieswise70.r
echo '}' >> pairwise_richness_specieswise70.r
## cladewise log plots
echo 'for(clade in 1:length(groups)){' >> pairwise_richness_specieswise70.r
echo 'jpeg(paste(groups[clade],"_log_transformed_richness_specieswise_morethanseventy.jpeg",sep=""),width=23.2,height=13,units="in",res=300)' >> pairwise_richness_specieswise70.r
echo 'b=lm(formula = log(finaldf$rich[finaldf$c==groups[clade]]) ~ log(finaldf$numrep[finaldf$c==groups[clade]]), data = finaldf)' >> pairwise_richness_specieswise70.r
echo 'plot(log(finaldf$numrep[finaldf$c==groups[clade]]),log(finaldf$rich[finaldf$c==groups[clade]]),ylab="log(richness)",xlab="log(number of sequences)",col=col13rainbow[clade],pch=20,main=paste("log transformed richness of >70% repeats and number of sequences in ", groups[clade],sep=""),cex=2)' >> pairwise_richness_specieswise70.r
echo 'abline(b,col=col13rainbow[clade],lty=2,lwd=2)' >> pairwise_richness_specieswise70.r
echo 'legend("topleft",legend=groups[clade],fill=col13rainbow,col=col13rainbow,bty="n")' >> pairwise_richness_specieswise70.r
echo 'dev.off()' >> pairwise_richness_specieswise70.r
echo '}' >> pairwise_richness_specieswise70.r
Rscript pairwise_richness_specieswise70.r


#################### For exact repeats ##########################################
## making boxplots for richness and density
###### making the diversity input file
awk '(($8-$7)+1)==$9{print $1,$3,$11}' OFS='\t' morethanseventy_combined|awk '{ gsub(/.*.[0-9]_/,"", $2); print }' OFS="\t" > exact_diversity_specieswise
for clades in afrotheria amphibia artiodactyla aves carnivore chiroptera lagomorpha marsupials perissodactyla primates rodents squamata testudines
do
grep "$clades" exact_diversity_specieswise |sort -k2|uniq -c|awk '{print $2,$3,$4,$1}' OFS='\t'|sort -k3nr >> temp.txt
done
mv temp.txt exact_diversity_specieswise
######### for diversity ###################
### making the diversity input file ###
##R
echo 'a=read.table("exact_diversity_specieswise",header=F)' > pairwise_species_diversity_richness_exact.r
echo 'colnames(a)=c("clades","species","repeat","count")' >> pairwise_species_diversity_richness_exact.r
echo 'library(stringr)' >> pairwise_species_diversity_richness_exact.r
echo 'a$clades=str_to_title(a$clades)' >> pairwise_species_diversity_richness_exact.r
echo '#install.packages("qpcR")' >> pairwise_species_diversity_richness_exact.r
echo 'library("qpcR")' >> pairwise_species_diversity_richness_exact.r
echo '####for Shannon diversity####' >> pairwise_species_diversity_richness_exact.r
echo 'finaldfr <- data.frame(first=c(1:10))' >> pairwise_species_diversity_richness_exact.r
echo 'for(c in unique(a$clades)){' >> pairwise_species_diversity_richness_exact.r
echo 'finaldf <- data.frame(clades=character(), shd=numeric())' >> pairwise_species_diversity_richness_exact.r
echo 'b1=a[a$clades==c,]' >> pairwise_species_diversity_richness_exact.r
echo 'for(sp in unique(b1$species)){' >> pairwise_species_diversity_richness_exact.r
echo 'b=b1[b1$species==sp,]' >> pairwise_species_diversity_richness_exact.r
echo 'b$V4=(b$count/sum(b$count))*(-log(b$count/sum(b$count)))' >> pairwise_species_diversity_richness_exact.r
echo 'sd1=sum(b$V4)' >> pairwise_species_diversity_richness_exact.r
echo 'sd=data.frame(sd1)' >> pairwise_species_diversity_richness_exact.r
echo 'finaldf=rbind(finaldf,sd)' >> pairwise_species_diversity_richness_exact.r
echo '}' >> pairwise_species_diversity_richness_exact.r
echo 'colnames(finaldf)=c' >> pairwise_species_diversity_richness_exact.r
echo 'finaldfr=qpcR:::cbind.na(finaldfr,finaldf)' >> pairwise_species_diversity_richness_exact.r
echo '}' >> pairwise_species_diversity_richness_exact.r
echo 'd=finaldfr[,-1]' >> pairwise_species_diversity_richness_exact.r
echo 'library(stringr)' >> pairwise_species_diversity_richness_exact.r
echo 'str_to_title(colnames(d))->colnames(d)' >> pairwise_species_diversity_richness_exact.r
echo 'species<-colnames(d)' >> pairwise_species_diversity_richness_exact.r
echo 'len<-length(species)' >> pairwise_species_diversity_richness_exact.r
echo 'totalcol=dim(d)[2]' >> pairwise_species_diversity_richness_exact.r
echo 'ymin=c(min(d[,1],na.rm=T))' >> pairwise_species_diversity_richness_exact.r
echo 'for ( i in 2:length(d)){' >> pairwise_species_diversity_richness_exact.r
echo 'ymin[i]<-c(min(d[,i],na.rm=T))}' >> pairwise_species_diversity_richness_exact.r
echo 'maxquan=0' >> pairwise_species_diversity_richness_exact.r
echo 'P<- data.frame(id = character(0), p.value = numeric(0))' >> pairwise_species_diversity_richness_exact.r
echo 'for ( i in seq(1,(totalcol-1),1)){' >> pairwise_species_diversity_richness_exact.r
echo 'for ( j in seq(i+1,totalcol,1)){' >> pairwise_species_diversity_richness_exact.r
echo 'if (species[i] == species[j]){' >> pairwise_species_diversity_richness_exact.r
echo 'next' >> pairwise_species_diversity_richness_exact.r
echo '}' >> pairwise_species_diversity_richness_exact.r
echo 'else' >> pairwise_species_diversity_richness_exact.r
echo '{' >> pairwise_species_diversity_richness_exact.r
echo 'print(paste(species[i],species[j],sep=" "))' >> pairwise_species_diversity_richness_exact.r
echo 'fg1<-d[,i]' >> pairwise_species_diversity_richness_exact.r
echo 'bg1<-d[,j]' >> pairwise_species_diversity_richness_exact.r
echo 'if (maxquan < quantile(fg1,0.75,na.rm=T)){maxquan=quantile(fg1,0.75,na.rm=T)}' >> pairwise_species_diversity_richness_exact.r
echo 'wilcox.test(fg1,bg1,alternative = "two.sided")->K' >> pairwise_species_diversity_richness_exact.r
echo 'P1<-data.frame(names=paste(species[i],"vs",species[j],sep="_"),p.value=K$p.value)' >> pairwise_species_diversity_richness_exact.r
echo 'seg2=paste(species[i],"vs",species[j],sep="_")' >> pairwise_species_diversity_richness_exact.r
echo 'P<-rbind(P,P1) }' >> pairwise_species_diversity_richness_exact.r
echo '}' >> pairwise_species_diversity_richness_exact.r
echo '}' >> pairwise_species_diversity_richness_exact.r
echo 'options(scipen=0,"digits"=2)' >> pairwise_species_diversity_richness_exact.r
echo 'P$qvalue<-p.adjust(P$p.value, method = "BH")' >> pairwise_species_diversity_richness_exact.r
echo 'P$V1<-format(P$qvalue,scientific=T)' >> pairwise_species_diversity_richness_exact.r
echo 'jpeg("Pairwise_boxplot_specieswise_Shannon_diversity_exactrepeats.jpeg",height=2500,width=3600)' >> pairwise_species_diversity_richness_exact.r
echo 'par(mai=c(5,5,2,1),mgp=c(18,0.5,0.1))' >> pairwise_species_diversity_richness_exact.r
echo 'boxplot(d,outline=F,main="Shannon diversity of exact repeats across different clades", xlab="Clades",ylab="Shannon diversity", col="Lightsalmon", border="black", las=2, cex=4, cex.main=4, cex.axis=3, cex.lab=4, lwd=2, ylim=c(min(ymin),as.numeric(maxquan)+1+(2*(length(P[P$qvalue<0.05,3])))))' >> pairwise_species_diversity_richness_exact.r
echo 'start=maxquan' >> pairwise_species_diversity_richness_exact.r
echo 'for ( i in seq(1,len,1)){' >> pairwise_species_diversity_richness_exact.r
echo 'for (j in seq(2,len,1)){' >> pairwise_species_diversity_richness_exact.r
echo 'tryCatch({' >> pairwise_species_diversity_richness_exact.r
echo 'if (species[i] == species[j]){' >> pairwise_species_diversity_richness_exact.r
echo 'next' >> pairwise_species_diversity_richness_exact.r
echo '}' >> pairwise_species_diversity_richness_exact.r
echo 'else' >> pairwise_species_diversity_richness_exact.r
echo '{' >> pairwise_species_diversity_richness_exact.r
echo 'fg=as.character(species[i])' >> pairwise_species_diversity_richness_exact.r
echo 'bg=as.character(species[j])' >> pairwise_species_diversity_richness_exact.r
echo 'seg=paste(fg,"vs",bg,sep="_")' >> pairwise_species_diversity_richness_exact.r
echo 'P[P$names==seg,3]-> qval' >> pairwise_species_diversity_richness_exact.r
echo 'P[P$names==seg,4] -> qvalchar' >> pairwise_species_diversity_richness_exact.r
echo 'if(qval > 0.05){' >> pairwise_species_diversity_richness_exact.r
echo 'next' >> pairwise_species_diversity_richness_exact.r
echo '}' >> pairwise_species_diversity_richness_exact.r
echo 'else' >> pairwise_species_diversity_richness_exact.r
echo '{' >> pairwise_species_diversity_richness_exact.r
echo 'segments(i,start,j,start,lty=3,col="purple",lwd=2)' >> pairwise_species_diversity_richness_exact.r
echo 'text(((i+j)/2),start+1,labels=qvalchar,cex=2)' >> pairwise_species_diversity_richness_exact.r
echo '}' >> pairwise_species_diversity_richness_exact.r
echo '}' >> pairwise_species_diversity_richness_exact.r
echo 'start = start+2' >> pairwise_species_diversity_richness_exact.r
echo '}' >> pairwise_species_diversity_richness_exact.r
echo ', error=function(e){})}' >> pairwise_species_diversity_richness_exact.r
echo '}' >> pairwise_species_diversity_richness_exact.r
echo 'dev.off()' >> pairwise_species_diversity_richness_exact.r
echo 'jpeg("Shannon_diversity_boxplot_specieswise_exactrepeats.jpeg",height=2500,width=3600)' >> pairwise_species_diversity_richness_exact.r
echo 'par(mai=c(5,5,2,1),mgp=c(18,0.5,0.1))' >> pairwise_species_diversity_richness_exact.r
echo 'boxplot(d,outline=F,main="Shannon diversity of exact repeats across different clades", xlab="Clades",ylab="Shannon diversity", col="Lightsalmon", border="black", las=2, cex=4, cex.main=4, cex.axis=3, cex.lab=4, lwd=2)' >> pairwise_species_diversity_richness_exact.r
echo 'dev.off()' >> pairwise_species_diversity_richness_exact.r
echo '####for Simpson diversity####' >> pairwise_species_diversity_richness_exact.r
echo 'finaldfr<-data.frame(first=c(1:10))' >> pairwise_species_diversity_richness_exact.r
echo 'for(c in unique(a$clades)){' >> pairwise_species_diversity_richness_exact.r
echo 'finaldf<-data.frame(clades=character(), shd=numeric())' >> pairwise_species_diversity_richness_exact.r
echo 's1=a[a$clades==c,]' >> pairwise_species_diversity_richness_exact.r
echo 'for(sp in unique(s1$species)){' >> pairwise_species_diversity_richness_exact.r
echo 's=s1[s1$species==sp,]' >> pairwise_species_diversity_richness_exact.r
echo 's$V4=s$count*(s$count - 1)' >> pairwise_species_diversity_richness_exact.r
echo 'simpd=round(1-(sum(s$V4)/(sum(s$count)*(sum(s$count)-1))),4)' >> pairwise_species_diversity_richness_exact.r
echo 'sd=data.frame(simpd)' >> pairwise_species_diversity_richness_exact.r
echo 'finaldf=rbind(finaldf,sd)' >> pairwise_species_diversity_richness_exact.r
echo '}' >> pairwise_species_diversity_richness_exact.r
echo 'colnames(finaldf)=c' >> pairwise_species_diversity_richness_exact.r
echo 'finaldfr=qpcR:::cbind.na(finaldfr,finaldf)' >> pairwise_species_diversity_richness_exact.r
echo '}' >> pairwise_species_diversity_richness_exact.r
echo 'd=finaldfr[,-1]' >> pairwise_species_diversity_richness_exact.r
echo 'library(stringr)' >> pairwise_species_diversity_richness_exact.r
echo 'str_to_title(colnames(d))->colnames(d)' >> pairwise_species_diversity_richness_exact.r
echo 'species<-colnames(d)' >> pairwise_species_diversity_richness_exact.r
echo 'len<-length(species)' >> pairwise_species_diversity_richness_exact.r
echo 'totalcol=dim(d)[2]' >> pairwise_species_diversity_richness_exact.r
echo 'ymin=c(min(d[,1],na.rm=T))' >> pairwise_species_diversity_richness_exact.r
echo 'for ( i in 2:length(d)){' >> pairwise_species_diversity_richness_exact.r
echo 'ymin[i]<-c(min(d[,i],na.rm=T))}' >> pairwise_species_diversity_richness_exact.r
echo 'maxquan=0' >> pairwise_species_diversity_richness_exact.r
echo 'P<- data.frame(id = character(0), p.value = numeric(0))' >> pairwise_species_diversity_richness_exact.r
echo 'for ( i in seq(1,(totalcol-1),1)){' >> pairwise_species_diversity_richness_exact.r
echo 'for ( j in seq(i+1,totalcol,1)){' >> pairwise_species_diversity_richness_exact.r
echo 'if (species[i] == species[j]){' >> pairwise_species_diversity_richness_exact.r
echo 'next' >> pairwise_species_diversity_richness_exact.r
echo '}' >> pairwise_species_diversity_richness_exact.r
echo 'else' >> pairwise_species_diversity_richness_exact.r
echo '{' >> pairwise_species_diversity_richness_exact.r
echo 'print(paste(species[i],species[j],sep=" "))' >> pairwise_species_diversity_richness_exact.r
echo 'fg1<-d[,i]' >> pairwise_species_diversity_richness_exact.r
echo 'bg1<-d[,j]' >> pairwise_species_diversity_richness_exact.r
echo 'if (maxquan < quantile(fg1,0.75,na.rm=T)){maxquan=quantile(fg1,0.75,na.rm=T) }' >> pairwise_species_diversity_richness_exact.r
echo 'wilcox.test(fg1,bg1,alternative = "two.sided")->K' >> pairwise_species_diversity_richness_exact.r
echo 'P1<-data.frame(names=paste(species[i],"vs",species[j],sep="_"),p.value=K$p.value)' >> pairwise_species_diversity_richness_exact.r
echo 'seg2=paste(species[i],"vs",species[j],sep="_")' >> pairwise_species_diversity_richness_exact.r
echo 'P<-rbind(P,P1) }' >> pairwise_species_diversity_richness_exact.r
echo '}' >> pairwise_species_diversity_richness_exact.r
echo '}' >> pairwise_species_diversity_richness_exact.r
echo 'options(scipen=0,"digits"=2)' >> pairwise_species_diversity_richness_exact.r
echo 'P$qvalue<-p.adjust(P$p.value, method = "BH")' >> pairwise_species_diversity_richness_exact.r
echo 'P$V1<-format(P$qvalue,scientific=T)' >> pairwise_species_diversity_richness_exact.r
echo 'jpeg("Pairwise_boxplot_Simpson_diversity_specieswise_exactrepeats.jpeg",height=2500,width=3600)' >> pairwise_species_diversity_richness_exact.r
echo 'par(mai=c(5,5,2,1),mgp=c(18,0.5,0.1))' >> pairwise_species_diversity_richness_exact.r
echo 'boxplot(d,outline=F,main="Simpson diversity of exact repeats across different clades", xlab="Clades",ylab="Simpson diversity", col="Lightsalmon", border="black", las=2, cex=4, cex.main=4, cex.axis=3, cex.lab=4, lwd=2, ylim=c(min(ymin),as.numeric(maxquan)+1+(2*(length(P[P$qvalue<0.05,3])))))' >> pairwise_species_diversity_richness_exact.r
echo 'start=maxquan' >> pairwise_species_diversity_richness_exact.r
echo 'for ( i in seq(1,len,1)){' >> pairwise_species_diversity_richness_exact.r
echo 'for (j in seq(2,len,1)){' >> pairwise_species_diversity_richness_exact.r
echo 'tryCatch({' >> pairwise_species_diversity_richness_exact.r
echo 'if (species[i] == species[j]){' >> pairwise_species_diversity_richness_exact.r
echo 'next' >> pairwise_species_diversity_richness_exact.r
echo '}' >> pairwise_species_diversity_richness_exact.r
echo 'else' >> pairwise_species_diversity_richness_exact.r
echo '{' >> pairwise_species_diversity_richness_exact.r
echo 'fg=as.character(species[i])' >> pairwise_species_diversity_richness_exact.r
echo 'bg=as.character(species[j])' >> pairwise_species_diversity_richness_exact.r
echo 'seg=paste(fg,"vs",bg,sep="_")' >> pairwise_species_diversity_richness_exact.r
echo 'P[P$names==seg,3]-> qval' >> pairwise_species_diversity_richness_exact.r
echo 'P[P$names==seg,4] -> qvalchar' >> pairwise_species_diversity_richness_exact.r
echo 'if(qval > 0.05){' >> pairwise_species_diversity_richness_exact.r
echo 'next' >> pairwise_species_diversity_richness_exact.r
echo '}' >> pairwise_species_diversity_richness_exact.r
echo 'else' >> pairwise_species_diversity_richness_exact.r
echo '{' >> pairwise_species_diversity_richness_exact.r
echo 'segments(i,start,j,start,lty=3,col="purple",lwd=2)' >> pairwise_species_diversity_richness_exact.r
echo 'text(((i+j)/2),start+1,labels=qvalchar,cex=2)' >> pairwise_species_diversity_richness_exact.r
echo '}' >> pairwise_species_diversity_richness_exact.r
echo '}' >> pairwise_species_diversity_richness_exact.r
echo 'start = start+2' >> pairwise_species_diversity_richness_exact.r
echo '}' >> pairwise_species_diversity_richness_exact.r
echo ', error=function(e){})}' >> pairwise_species_diversity_richness_exact.r
echo '}' >> pairwise_species_diversity_richness_exact.r
echo 'dev.off()' >> pairwise_species_diversity_richness_exact.r
echo 'jpeg("Simpson_diversity_boxplot_specieswise_exactrepeats.jpeg",height=2500,width=3600)' >> pairwise_species_diversity_richness_exact.r
echo 'par(mai=c(5,5,2,1),mgp=c(18,0.5,0.1))' >> pairwise_species_diversity_richness_exact.r
echo 'boxplot(d,outline=F,main="Simpson diversity of exact repeats across different clades", xlab="Clades",ylab="Simpson diversity", col="Lightsalmon", border="black", las=2, cex=4, cex.main=4, cex.axis=3, cex.lab=4, lwd=2)' >> pairwise_species_diversity_richness_exact.r
echo 'dev.off()' >> pairwise_species_diversity_richness_exact.r
echo '############################################ calculating richness cladewise for exact repeats ###############################' >> pairwise_species_diversity_richness_exact.r
echo '############## R #################' >> pairwise_species_diversity_richness_exact.r
echo '####Richness and slope of richness' >> pairwise_species_diversity_richness_exact.r
echo 'a=read.table("exact_diversity_specieswise",header=F)' >> pairwise_species_diversity_richness_exact.r
echo 'colnames(a)=c("clades","species","repeat","count")' >> pairwise_species_diversity_richness_exact.r
echo 'library(stringr)' >> pairwise_species_diversity_richness_exact.r
echo 'a$clades=str_to_title(a$clades)' >> pairwise_species_diversity_richness_exact.r
echo 'finaldfr<-data.frame(first=c(1:10))' >> pairwise_species_diversity_richness_exact.r
echo 'for(c in unique(a$clades)){' >> pairwise_species_diversity_richness_exact.r
echo 'finaldf<-data.frame(clades=character(), shd=numeric())' >> pairwise_species_diversity_richness_exact.r
echo 's1=a[a$clades==c,]' >> pairwise_species_diversity_richness_exact.r
echo 'for(sp in unique(s1$species)){' >> pairwise_species_diversity_richness_exact.r
echo 's=s1[s1$species==sp,]' >> pairwise_species_diversity_richness_exact.r
echo 'simpd=dim(s)[1]' >> pairwise_species_diversity_richness_exact.r
echo 'sd=data.frame(simpd)' >> pairwise_species_diversity_richness_exact.r
echo 'finaldf=rbind(finaldf,sd)' >> pairwise_species_diversity_richness_exact.r
echo '}' >> pairwise_species_diversity_richness_exact.r
echo 'colnames(finaldf)=c' >> pairwise_species_diversity_richness_exact.r
echo 'finaldfr=qpcR:::cbind.na(finaldfr,finaldf)' >> pairwise_species_diversity_richness_exact.r
echo '}' >> pairwise_species_diversity_richness_exact.r
echo 'd=finaldfr[,-1]' >> pairwise_species_diversity_richness_exact.r
echo 'library(stringr)' >> pairwise_species_diversity_richness_exact.r
echo 'str_to_title(colnames(d))->colnames(d)' >> pairwise_species_diversity_richness_exact.r
echo 'species<-colnames(d)' >> pairwise_species_diversity_richness_exact.r
echo 'len<-length(species)' >> pairwise_species_diversity_richness_exact.r
echo 'totalcol=dim(d)[2]' >> pairwise_species_diversity_richness_exact.r
echo 'ymin=c(min(d[,1],na.rm=T))' >> pairwise_species_diversity_richness_exact.r
echo 'for ( i in 2:length(d)){' >> pairwise_species_diversity_richness_exact.r
echo 'ymin[i]<-c(min(d[,i],na.rm=T))}' >> pairwise_species_diversity_richness_exact.r
echo 'maxquan=0' >> pairwise_species_diversity_richness_exact.r
echo 'P<- data.frame(id = character(0), p.value = numeric(0))' >> pairwise_species_diversity_richness_exact.r
echo 'for ( i in seq(1,(totalcol-1),1)){' >> pairwise_species_diversity_richness_exact.r
echo 'for ( j in seq(i+1,totalcol,1)){' >> pairwise_species_diversity_richness_exact.r
echo 'if (species[i] == species[j]){' >> pairwise_species_diversity_richness_exact.r
echo 'next' >> pairwise_species_diversity_richness_exact.r
echo '}' >> pairwise_species_diversity_richness_exact.r
echo 'else' >> pairwise_species_diversity_richness_exact.r
echo '{' >> pairwise_species_diversity_richness_exact.r
echo 'print(paste(species[i],species[j],sep=" "))' >> pairwise_species_diversity_richness_exact.r
echo 'fg1<-d[,i]' >> pairwise_species_diversity_richness_exact.r
echo 'bg1<-d[,j]' >> pairwise_species_diversity_richness_exact.r
echo 'if (maxquan < quantile(fg1,0.75,na.rm=T)){maxquan=quantile(fg1,0.75,na.rm=T) }' >> pairwise_species_diversity_richness_exact.r
echo 'wilcox.test(fg1,bg1,alternative = "two.sided")->K' >> pairwise_species_diversity_richness_exact.r
echo 'P1<-data.frame(names=paste(species[i],"vs",species[j],sep="_"),p.value=K$p.value)' >> pairwise_species_diversity_richness_exact.r
echo 'seg2=paste(species[i],"vs",species[j],sep="_")' >> pairwise_species_diversity_richness_exact.r
echo 'P<-rbind(P,P1) }' >> pairwise_species_diversity_richness_exact.r
echo '}' >> pairwise_species_diversity_richness_exact.r
echo '}' >> pairwise_species_diversity_richness_exact.r
echo 'options(scipen=0,"digits"=2)' >> pairwise_species_diversity_richness_exact.r
echo 'P$qvalue<-p.adjust(P$p.value, method = "BH")' >> pairwise_species_diversity_richness_exact.r
echo 'P$V1<-format(P$qvalue,scientific=T)' >> pairwise_species_diversity_richness_exact.r
echo 'jpeg("Pairwise_boxplot_specieswise_richness_exactrepeats.jpeg",height=2500,width=3600)' >> pairwise_species_diversity_richness_exact.r
echo 'par(mai=c(5,5,2,1),mgp=c(18,0.5,0.1))' >> pairwise_species_diversity_richness_exact.r
echo 'boxplot(d,outline=F,main="Richness of exact repeats across clades", xlab="Clades",ylab="Richness of exact repeats", col="Lightsalmon", border="black", las=2, cex=4, cex.main=4, cex.axis=3, cex.lab=4, lwd=2, ylim=c(min(ymin),as.numeric(maxquan)+1+(2*(length(P[P$qvalue<0.05,3])))))' >> pairwise_species_diversity_richness_exact.r
echo 'start=maxquan' >> pairwise_species_diversity_richness_exact.r
echo 'for ( i in seq(1,len,1)){' >> pairwise_species_diversity_richness_exact.r
echo 'for (j in seq(2,len,1)){' >> pairwise_species_diversity_richness_exact.r
echo 'tryCatch({' >> pairwise_species_diversity_richness_exact.r
echo 'if (species[i] == species[j]){' >> pairwise_species_diversity_richness_exact.r
echo 'next' >> pairwise_species_diversity_richness_exact.r
echo '}' >> pairwise_species_diversity_richness_exact.r
echo 'else' >> pairwise_species_diversity_richness_exact.r
echo '{' >> pairwise_species_diversity_richness_exact.r
echo 'fg=as.character(species[i])' >> pairwise_species_diversity_richness_exact.r
echo 'bg=as.character(species[j])' >> pairwise_species_diversity_richness_exact.r
echo 'seg=paste(fg,"vs",bg,sep="_")' >> pairwise_species_diversity_richness_exact.r
echo 'P[P$names==seg,3]-> qval' >> pairwise_species_diversity_richness_exact.r
echo 'P[P$names==seg,4] -> qvalchar' >> pairwise_species_diversity_richness_exact.r
echo 'if(qval > 0.05){' >> pairwise_species_diversity_richness_exact.r
echo 'next' >> pairwise_species_diversity_richness_exact.r
echo '}' >> pairwise_species_diversity_richness_exact.r
echo 'else' >> pairwise_species_diversity_richness_exact.r
echo '{' >> pairwise_species_diversity_richness_exact.r
echo 'segments(i,start,j,start,lty=3,col="purple",lwd=2)' >> pairwise_species_diversity_richness_exact.r
echo 'text(((i+j)/2),start+1,labels=qvalchar,cex=2)' >> pairwise_species_diversity_richness_exact.r
echo '}' >> pairwise_species_diversity_richness_exact.r
echo '}' >> pairwise_species_diversity_richness_exact.r
echo 'start = start+2' >> pairwise_species_diversity_richness_exact.r
echo '}' >> pairwise_species_diversity_richness_exact.r
echo ', error=function(e){})}' >> pairwise_species_diversity_richness_exact.r
echo '}' >> pairwise_species_diversity_richness_exact.r
echo 'dev.off()' >> pairwise_species_diversity_richness_exact.r
echo 'jpeg("Boxplot_specieswise_richness_exactrepeats.jpeg",height=2500,width=3600)' >> pairwise_species_diversity_richness_exact.r
echo 'par(mai=c(5,5,2,1),mgp=c(18,0.5,0.1))' >> pairwise_species_diversity_richness_exact.r
echo 'boxplot(d,outline=F,main="Richness of exact repeats across clades", xlab="Clades",ylab="Richness of exact repeats", col="Lightsalmon", border="black", las=2, cex=4, cex.main=4, cex.axis=3, cex.lab=4, lwd=2)' >> pairwise_species_diversity_richness_exact.r
echo 'dev.off()' >> pairwise_species_diversity_richness_exact.r
echo '####### for slope (normalized) #########' >> pairwise_species_diversity_richness_exact.r
echo 'a=read.table("exact_diversity_specieswise",header=F)' >> pairwise_species_diversity_richness_exact.r
echo 'colnames(a)=c("clades","species","repeat","count")' >> pairwise_species_diversity_richness_exact.r
echo 'library(stringr)' >> pairwise_species_diversity_richness_exact.r
echo 'a$clades=str_to_title(a$clades)' >> pairwise_species_diversity_richness_exact.r
echo 'finaldfr <- data.frame(first=c(1:10))' >> pairwise_species_diversity_richness_exact.r
echo 'for(c in unique(a$clades)){' >> pairwise_species_diversity_richness_exact.r
echo 'finaldf <- data.frame(clades=character(), shd=numeric())' >> pairwise_species_diversity_richness_exact.r
echo 'b1=a[a$clades==c,]' >> pairwise_species_diversity_richness_exact.r
echo 'for(sp in unique(b1$species)){' >> pairwise_species_diversity_richness_exact.r
echo 'b=b1[b1$species==sp,]' >> pairwise_species_diversity_richness_exact.r
echo 'rich1=length(a[a$species==sp,4])' >> pairwise_species_diversity_richness_exact.r
echo 'numrep=sum(a[a$species==sp,4])' >> pairwise_species_diversity_richness_exact.r
echo 'slope=log(rich1)/log(numrep)' >> pairwise_species_diversity_richness_exact.r
echo 'rich=data.frame(slope)' >> pairwise_species_diversity_richness_exact.r
echo 'finaldf=rbind(finaldf,rich)' >> pairwise_species_diversity_richness_exact.r
echo '}' >> pairwise_species_diversity_richness_exact.r
echo 'colnames(finaldf)=c' >> pairwise_species_diversity_richness_exact.r
echo 'finaldfr=qpcR:::cbind.na(finaldfr,finaldf)' >> pairwise_species_diversity_richness_exact.r
echo '}' >> pairwise_species_diversity_richness_exact.r
echo 'd=finaldfr[,-1]' >> pairwise_species_diversity_richness_exact.r
echo 'library(stringr)' >> pairwise_species_diversity_richness_exact.r
echo 'str_to_title(colnames(d))->colnames(d)' >> pairwise_species_diversity_richness_exact.r
echo 'species<-colnames(d)' >> pairwise_species_diversity_richness_exact.r
echo 'len<-length(species)' >> pairwise_species_diversity_richness_exact.r
echo 'totalcol=dim(d)[2]' >> pairwise_species_diversity_richness_exact.r
echo 'ymin=c(min(d[,1],na.rm=T))' >> pairwise_species_diversity_richness_exact.r
echo 'for ( i in 2:length(d)){' >> pairwise_species_diversity_richness_exact.r
echo 'ymin[i]<-c(min(d[,i],na.rm=T))}' >> pairwise_species_diversity_richness_exact.r
echo 'maxquan=0' >> pairwise_species_diversity_richness_exact.r
echo 'P<- data.frame(id = character(0), p.value = numeric(0))' >> pairwise_species_diversity_richness_exact.r
echo 'for ( i in seq(1,(totalcol-1),1)){' >> pairwise_species_diversity_richness_exact.r
echo 'for ( j in seq(i+1,totalcol,1)){' >> pairwise_species_diversity_richness_exact.r
echo 'if (species[i] == species[j]){' >> pairwise_species_diversity_richness_exact.r
echo 'next' >> pairwise_species_diversity_richness_exact.r
echo '}' >> pairwise_species_diversity_richness_exact.r
echo 'else' >> pairwise_species_diversity_richness_exact.r
echo '{' >> pairwise_species_diversity_richness_exact.r
echo 'print(paste(species[i],species[j],sep=" "))' >> pairwise_species_diversity_richness_exact.r
echo 'fg1<-d[,i]' >> pairwise_species_diversity_richness_exact.r
echo 'bg1<-d[,j]' >> pairwise_species_diversity_richness_exact.r
echo 'if (maxquan < quantile(fg1,0.75,na.rm=T)){maxquan=quantile(fg1,0.75,na.rm=T) }' >> pairwise_species_diversity_richness_exact.r
echo 'wilcox.test(fg1,bg1,alternative = "two.sided")->K' >> pairwise_species_diversity_richness_exact.r
echo 'P1<-data.frame(names=paste(species[i],"vs",species[j],sep="_"),p.value=K$p.value)' >> pairwise_species_diversity_richness_exact.r
echo 'seg2=paste(species[i],"vs",species[j],sep="_")' >> pairwise_species_diversity_richness_exact.r
echo 'P<-rbind(P,P1) }' >> pairwise_species_diversity_richness_exact.r
echo '}' >> pairwise_species_diversity_richness_exact.r
echo '}' >> pairwise_species_diversity_richness_exact.r
echo 'options(scipen=0,"digits"=2)' >> pairwise_species_diversity_richness_exact.r
echo 'P$qvalue<-p.adjust(P$p.value, method = "BH")' >> pairwise_species_diversity_richness_exact.r
echo 'P$V1<-format(P$qvalue,scientific=T)' >> pairwise_species_diversity_richness_exact.r
echo 'jpeg("Pairwise_boxplot_specieswise_richness_slope_exactrepeats.jpeg",height=2500,width=3600)' >> pairwise_species_diversity_richness_exact.r
echo 'par(mai=c(5,5,2,1),mgp=c(18,0.5,0.1))' >> pairwise_species_diversity_richness_exact.r
echo 'boxplot(d,outline=F,main="Slope of richness of exact repeats across clades", xlab="Clades",ylab="Slope of richness of exact repeats", col="Lightsalmon", border="black", las=2, cex=4, cex.main=4, cex.axis=3, cex.lab=4, lwd=2, ylim=c(min(ymin),as.numeric(maxquan)+1+(2*(length(P[P$qvalue<0.05,3])))))' >> pairwise_species_diversity_richness_exact.r
echo 'start=maxquan' >> pairwise_species_diversity_richness_exact.r
echo 'for ( i in seq(1,len,1)){' >> pairwise_species_diversity_richness_exact.r
echo 'for (j in seq(2,len,1)){' >> pairwise_species_diversity_richness_exact.r
echo 'tryCatch({' >> pairwise_species_diversity_richness_exact.r
echo 'if (species[i] == species[j]){' >> pairwise_species_diversity_richness_exact.r
echo 'next' >> pairwise_species_diversity_richness_exact.r
echo '}' >> pairwise_species_diversity_richness_exact.r
echo 'else' >> pairwise_species_diversity_richness_exact.r
echo '{' >> pairwise_species_diversity_richness_exact.r
echo 'fg=as.character(species[i])' >> pairwise_species_diversity_richness_exact.r
echo 'bg=as.character(species[j])' >> pairwise_species_diversity_richness_exact.r
echo 'seg=paste(fg,"vs",bg,sep="_")' >> pairwise_species_diversity_richness_exact.r
echo 'P[P$names==seg,3]-> qval' >> pairwise_species_diversity_richness_exact.r
echo 'P[P$names==seg,4] -> qvalchar' >> pairwise_species_diversity_richness_exact.r
echo 'if(qval > 0.05){' >> pairwise_species_diversity_richness_exact.r
echo 'next' >> pairwise_species_diversity_richness_exact.r
echo '}' >> pairwise_species_diversity_richness_exact.r
echo 'else' >> pairwise_species_diversity_richness_exact.r
echo '{' >> pairwise_species_diversity_richness_exact.r
echo 'segments(i,start,j,start,lty=3,col="purple",lwd=2)' >> pairwise_species_diversity_richness_exact.r
echo 'text(((i+j)/2),start+1,labels=qvalchar,cex=2)' >> pairwise_species_diversity_richness_exact.r
echo '}' >> pairwise_species_diversity_richness_exact.r
echo '}' >> pairwise_species_diversity_richness_exact.r
echo 'start = start+2' >> pairwise_species_diversity_richness_exact.r
echo '}' >> pairwise_species_diversity_richness_exact.r
echo ', error=function(e){})}' >> pairwise_species_diversity_richness_exact.r
echo '}' >> pairwise_species_diversity_richness_exact.r
echo 'dev.off()' >> pairwise_species_diversity_richness_exact.r
echo 'jpeg("Boxplot_specieswise_richness_slope_exactrepeats.jpeg",height=2500,width=3600)' >> pairwise_species_diversity_richness_exact.r
echo 'par(mai=c(5,5,2,1),mgp=c(18,0.5,0.1))' >> pairwise_species_diversity_richness_exact.r
echo 'boxplot(d,outline=F,main="Slope of richness of exact repeats across clades", xlab="Clades",ylab="Slope of richness of exact repeats", col="Lightsalmon", border="black", las=2, cex=4, cex.main=4, cex.axis=3, cex.lab=4, lwd=2)' >> pairwise_species_diversity_richness_exact.r
echo 'dev.off()' >> pairwise_species_diversity_richness_exact.r
echo '########### slope of richness ###log-transformed richness and number of sequences #############' >> pairwise_species_diversity_richness_exact.r
echo 'finaldf <- data.frame(clades=character(), species=character(), rich=numeric(), numrep=numeric())' >> pairwise_species_diversity_richness_exact.r
echo 'for(c in unique(a$clades)){' >> pairwise_species_diversity_richness_exact.r
echo 'b1=a[a$clades==c,]' >> pairwise_species_diversity_richness_exact.r
echo 'for(sp in unique(b1$species)){' >> pairwise_species_diversity_richness_exact.r
echo 'b=b1[b1$species==sp,]' >> pairwise_species_diversity_richness_exact.r
echo 'rich1=length(a[a$species==sp,4])' >> pairwise_species_diversity_richness_exact.r
echo 'numrep=sum(a[a$species==sp,4])' >> pairwise_species_diversity_richness_exact.r
echo 'slope=log(rich1)/log(numrep)' >> pairwise_species_diversity_richness_exact.r
echo 'rich=data.frame(c,sp,rich1,numrep)' >> pairwise_species_diversity_richness_exact.r
echo 'finaldf=rbind(finaldf,rich)' >> pairwise_species_diversity_richness_exact.r
echo '}' >> pairwise_species_diversity_richness_exact.r
echo '}' >> pairwise_species_diversity_richness_exact.r
echo 'col13rainbow= c("#771155", "#AA4488", "#CC99BB", "#114477", "#4477AA", "#77AADD", "#117777", "#44AAAA", "#77CCCC", "#117744", "#44AA77", "#88CCAA", "#777711")' >> pairwise_species_diversity_richness_exact.r
echo 'groups=as.character(unique(finaldf$c))' >> pairwise_species_diversity_richness_exact.r
echo '## all clades in one-log plots' >> pairwise_species_diversity_richness_exact.r
echo 'clade=groups[1]' >> pairwise_species_diversity_richness_exact.r
echo 'jpeg("log_transformed_specieswise_richness_exactrepeats.jpeg",width=23.2,height=13,units="in",res=300)' >> pairwise_species_diversity_richness_exact.r
echo 'b=lm(formula = log(finaldf$rich[finaldf$c==clade]) ~ log(finaldf$numrep[finaldf$c==clade]), data = finaldf)' >> pairwise_species_diversity_richness_exact.r
echo 'plot(log(finaldf$numrep[finaldf$c==clade]),log(finaldf$rich[finaldf$c==clade]),ylab="log(richness)",xlab="log(number of sequences)",col=col13rainbow[1],pch=20,main="log transformed richness of exact repeats and number of sequences across different clades",cex=2,xlim=c(2.5,6),ylim=c(2,4))' >> pairwise_species_diversity_richness_exact.r
echo 'abline(b,col=col13rainbow[1],lty=2,lwd=2)' >> pairwise_species_diversity_richness_exact.r
echo 'for(clade in 2:length(groups)){' >> pairwise_species_diversity_richness_exact.r
echo 'b=lm(formula = log(finaldf$rich[finaldf$c==groups[clade]]) ~ log(finaldf$numrep[finaldf$c==groups[clade]]), data = finaldf)' >> pairwise_species_diversity_richness_exact.r
echo 'points(log(finaldf$numrep[finaldf$c==groups[clade]]),log(finaldf$rich[finaldf$c==groups[clade]]),ylab="log(richness)",xlab="log(number of sequences)",col=col13rainbow[clade],pch=20,cex=2)' >> pairwise_species_diversity_richness_exact.r
echo 'abline(b,col=col13rainbow[clade],lty=2,lwd=2)' >> pairwise_species_diversity_richness_exact.r
echo '}' >> pairwise_species_diversity_richness_exact.r
echo 'legend("topleft",legend=groups,fill=col13rainbow,col=col13rainbow,bty="n")' >> pairwise_species_diversity_richness_exact.r
echo 'dev.off()' >> pairwise_species_diversity_richness_exact.r
echo '## cladewise log plots' >> pairwise_species_diversity_richness_exact.r
echo 'for(clade in 1:length(groups)){' >> pairwise_species_diversity_richness_exact.r
echo 'jpeg(paste(groups[clade],"_log_transformed_specieswise_richness_exactrepeats.jpeg",sep=""),width=23.2,height=13,units="in",res=300)' >> pairwise_species_diversity_richness_exact.r
echo 'b=lm(formula = log(finaldf$rich[finaldf$c==groups[clade]]) ~ log(finaldf$numrep[finaldf$c==groups[clade]]), data = finaldf)' >> pairwise_species_diversity_richness_exact.r
echo 'plot(log(finaldf$numrep[finaldf$c==groups[clade]]),log(finaldf$rich[finaldf$c==groups[clade]]),ylab="log(richness)",xlab="log(number of sequences)",col=col13rainbow[clade],pch=20,main=paste("log transformed richness of exact repeats and number of sequences in ", groups[clade],sep=""),cex=2)' >> pairwise_species_diversity_richness_exact.r
echo 'abline(b,col=col13rainbow[clade],lty=2,lwd=2)' >> pairwise_species_diversity_richness_exact.r
echo 'legend("topleft",legend=groups[clade],fill=col13rainbow[clade],col=col13rainbow[clade],bty="n")' >> pairwise_species_diversity_richness_exact.r
echo 'dev.off()' >> pairwise_species_diversity_richness_exact.r
echo '}' >> pairwise_species_diversity_richness_exact.r
Rscript pairwise_species_diversity_richness_exact.r
