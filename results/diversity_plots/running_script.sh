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
sed -i 's/afrotheria/Afrotheria/g' morethanseventy_diversity
sed -i 's/amphibia/Amphibia/g' morethanseventy_diversity
sed -i 's/artiodactyla/Artiodactyla/g' morethanseventy_diversity
sed -i 's/aves/Aves/g' morethanseventy_diversity
sed -i 's/carnivore/Carnivora/g' morethanseventy_diversity
sed -i 's/chiroptera/Chiroptera/g' morethanseventy_diversity
sed -i 's/lagomorpha/Lagomorpha/g' morethanseventy_diversity
sed -i 's/marsupials/Marsupialia/g' morethanseventy_diversity
sed -i 's/perissodactyla/Perissodactyla/g' morethanseventy_diversity
sed -i 's/primates/Primates/g' morethanseventy_diversity
sed -i 's/rodents/Rodentia/g' morethanseventy_diversity
sed -i 's/squamata/Squamata/g' morethanseventy_diversity
sed -i 's/testudines/Testudines/g' morethanseventy_diversity

sed -i 's/afrotheria/Afrotheria/g' numberoflps.txt
sed -i 's/amphibia/Amphibia/g' numberoflps.txt
sed -i 's/artiodactyla/Artiodactyla/g' numberoflps.txt
sed -i 's/aves/Aves/g' numberoflps.txt
sed -i 's/carnivore/Carnivora/g' numberoflps.txt
sed -i 's/chiroptera/Chiroptera/g' numberoflps.txt
sed -i 's/lagomorpha/Lagomorpha/g' numberoflps.txt
sed -i 's/marsupials/Marsupialia/g' numberoflps.txt
sed -i 's/perissodactyla/Perissodactyla/g' numberoflps.txt
sed -i 's/primates/Primates/g' numberoflps.txt
sed -i 's/rodents/Rodentia/g' numberoflps.txt
sed -i 's/squamata/Squamata/g' numberoflps.txt
sed -i 's/testudines/Testudines/g' numberoflps.txt

################Plotting in R#############
echo 'a=read.table("morethanseventy_diversity",header=F)' > diversity_richness70.r
echo 'colnames(a)=c("clades","repeat","count")' >> diversity_richness70.r
echo 'library(stringr)' >> diversity_richness70.r
echo 'a$clades=str_to_title(a$clades)' >> diversity_richness70.r
echo 'a$clades <- factor(a$clades,)' >> diversity_richness70.r
echo 'finaldf <- data.frame(clades=character(), shd=numeric())' >> diversity_richness70.r
echo 'for(c in unique(a$clades)){' >> diversity_richness70.r
echo 'b=a[a$clades==c,]' >> diversity_richness70.r
echo 'b$V4=(b$count/sum(b$count))*(-log(b$count/sum(b$count)))' >> diversity_richness70.r
echo 'sd=sum(b$V4)' >> diversity_richness70.r
echo 'd=data.frame(c,sd)' >> diversity_richness70.r
echo 'finaldf=rbind(finaldf,d)' >> diversity_richness70.r
echo '}' >> diversity_richness70.r
echo 'mf=data.frame()' >> diversity_richness70.r
echo 'for(clade in c("Afrotheria","Primates","Rodentia","Lagomorpha","Chiroptera","Artiodactyla","Perissodactyla","Carnivora","Marsupialia","Testudines","Aves","Squamata","Amphibia")){' >> diversity_richness70.r
echo 'm1=finaldf[finaldf$c==clade,]' >> diversity_richness70.r
echo 'mf=rbind(mf,m1)' >> diversity_richness70.r
echo '}' >> diversity_richness70.r
echo 'finaldf=mf' >> diversity_richness70.r
echo 'jpeg("shannon_diversity_morethanseventy.jpeg",width=23.2,height=13,units="in",res=300)' >> diversity_richness70.r
echo 'barplot(finaldf$sd,names=finaldf$c,ylab="Shannon diversity",ylim=c(0,max(finaldf$sd)),xlab="Clades",col="lightcyan",main="Shannon diversity of >70% LPSs across different clades")' >> diversity_richness70.r
echo 'dev.off()' >> diversity_richness70.r
echo 'simpsondf <- data.frame(clades=character(), simd=numeric())' >> diversity_richness70.r
echo 'for(c in unique(a$clades)){' >> diversity_richness70.r
echo 's=a[a$clades==c,]' >> diversity_richness70.r
echo 's$V4=s$count*(s$count - 1)' >> diversity_richness70.r
echo 'simpd=round(1-(sum(s$V4)/(sum(s$count)*(sum(s$count)-1))),4)' >> diversity_richness70.r
echo 'simp=data.frame(c,simpd)' >> diversity_richness70.r
echo 'simpsondf=rbind(simpsondf,simp)}' >> diversity_richness70.r
echo 'sf=data.frame()' >> diversity_richness70.r
echo 'for(clade in c("Afrotheria","Primates","Rodentia","Lagomorpha","Chiroptera","Artiodactyla","Perissodactyla","Carnivora","Marsupialia","Testudines","Aves","Squamata","Amphibia")){' >> diversity_richness70.r
echo 's1=simpsondf[simpsondf$c==clade,]' >> diversity_richness70.r
echo 'sf=rbind(sf,s1)' >> diversity_richness70.r
echo '}' >> diversity_richness70.r
echo 'simpsondf=sf' >> diversity_richness70.r
echo 'jpeg("simpson_diversity_morethanseventy.jpeg",width=23.2,height=13,units="in",res=300)' >> diversity_richness70.r
echo 'barplot(simpsondf$simpd,names=simpsondf$c,ylab="Simpson diversity",ylim=c(0,1),xlab="Clades",col="lightcyan",main="Simpson diversity of >70% LPSs across different clades")' >> diversity_richness70.r
echo 'dev.off()' >> diversity_richness70.r
echo 'unnorrich <- data.frame(clades=character(), richness=numeric())' >> diversity_richness70.r
echo 'for(c in unique(a$clades)){' >> diversity_richness70.r
echo 'b=a[a$clades==c,]' >> diversity_richness70.r
echo 'numrep=dim(b)[1]' >> diversity_richness70.r
echo 'ric=data.frame(c,numrep)' >> diversity_richness70.r
echo 'unnorrich=rbind(unnorrich,ric)}' >> diversity_richness70.r
echo 'unno=data.frame()' >> diversity_richness70.r
echo 'for(clade in c("Afrotheria","Primates","Rodentia","Lagomorpha","Chiroptera","Artiodactyla","Perissodactyla","Carnivora","Marsupialia","Testudines","Aves","Squamata","Amphibia")){' >> diversity_richness70.r
echo 's1=unnorrich[unnorrich$c==clade,]' >> diversity_richness70.r
echo 'unno=rbind(unno,s1)' >> diversity_richness70.r
echo '}' >> diversity_richness70.r
echo 'unnorrich=unno' >> diversity_richness70.r
echo 'jpeg("unnormalized_repeat_richness_morethanseventy.jpeg",width=23.2,height=13,units="in",res=300)' >> diversity_richness70.r
echo 'barplot(unnorrich$numrep,names=unnorrich$c,ylab="Repeat richness",xlab="Clades",ylim=c(0,max(unnorrich$numrep)+2),col="#69b3a2",main="Repeats richness for >70% LPSs across different clades")' >> diversity_richness70.r
echo 'dev.off()' >> diversity_richness70.r
echo 'numlps=read.table("numberoflps.txt",header=F)' >> diversity_richness70.r
echo 'numlps$V1=str_to_title(numlps$V1)' >> diversity_richness70.r
echo 'norrich <- data.frame(clades=character(), richness=numeric())' >> diversity_richness70.r
echo 'for(c in unique(a$clades)){' >> diversity_richness70.r
echo 'rich=unnorrich[unnorrich$c==c,2]' >> diversity_richness70.r
echo 'numrep=numlps[numlps$V1==c,2]' >> diversity_richness70.r
echo 'slope=log(rich)/log(numrep)' >> diversity_richness70.r
echo 'rich=data.frame(c,slope)' >> diversity_richness70.r
echo 'norrich=rbind(norrich,rich)}' >> diversity_richness70.r
echo 'nono=data.frame()' >> diversity_richness70.r
echo 'for(clade in c("Afrotheria","Primates","Rodentia","Lagomorpha","Chiroptera","Artiodactyla","Perissodactyla","Carnivora","Marsupialia","Testudines","Aves","Squamata","Amphibia")){' >> diversity_richness70.r
echo 's1=norrich[norrich$c==clade,]' >> diversity_richness70.r
echo 'nono=rbind(nono,s1)' >> diversity_richness70.r
echo '}' >> diversity_richness70.r
echo 'norrich=nono' >> diversity_richness70.r
echo 'jpeg("normalized_repeat_richness_morethanseventy.jpeg",width=23.2,height=13,units="in",res=300)' >> diversity_richness70.r
echo 'barplot(norrich$slope,names=norrich$c,ylab="Repeat richness",xlab="Clades",ylim=c(0,max(norrich$slope)+0.2),col="#69b3a2",main="Normalized repeats richness for >70% LPSs across different clades")' >> diversity_richness70.r
echo 'dev.off()' >> diversity_richness70.r
norrich <- data.frame(clades=character(), richness=numeric(),numrepeats=numeric())
for(c in unique(a$clades)){
rich=unnorrich[unnorrich$c==c,2]
numrep=numlps[numlps$V1==c,2]
rich=data.frame(c,rich,numrep)
norrich=rbind(norrich,rich)}
echo 'nono=data.frame()' >> diversity_richness70.r
echo 'for(clade in c("Afrotheria","Primates","Rodentia","Lagomorpha","Chiroptera","Artiodactyla","Perissodactyla","Carnivora","Marsupialia","Testudines","Aves","Squamata","Amphibia")){' >> diversity_richness70.r
echo 's1=norrich[norrich$c==clade,]' >> diversity_richness70.r
echo 'nono=rbind(nono,s1)' >> diversity_richness70.r
echo '}' >> diversity_richness70.r
echo 'norrich=nono' >> diversity_richness70.r
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

sed -i 's/afrotheria/Afrotheria/g' exactstretch_diversity
sed -i 's/amphibia/Amphibia/g' exactstretch_diversity
sed -i 's/artiodactyla/Artiodactyla/g' exactstretch_diversity
sed -i 's/aves/Aves/g' exactstretch_diversity
sed -i 's/carnivore/Carnivora/g' exactstretch_diversity
sed -i 's/chiroptera/Chiroptera/g' exactstretch_diversity
sed -i 's/lagomorpha/Lagomorpha/g' exactstretch_diversity
sed -i 's/marsupials/Marsupialia/g' exactstretch_diversity
sed -i 's/perissodactyla/Perissodactyla/g' exactstretch_diversity
sed -i 's/primates/Primates/g' exactstretch_diversity
sed -i 's/rodents/Rodentia/g' exactstretch_diversity
sed -i 's/squamata/Squamata/g' exactstretch_diversity
sed -i 's/testudines/Testudines/g' exactstretch_diversity

sed -i 's/afrotheria/Afrotheria/g' numberoflps_exact.txt
sed -i 's/amphibia/Amphibia/g' numberoflps_exact.txt
sed -i 's/artiodactyla/Artiodactyla/g' numberoflps_exact.txt
sed -i 's/aves/Aves/g' numberoflps_exact.txt
sed -i 's/carnivore/Carnivora/g' numberoflps_exact.txt
sed -i 's/chiroptera/Chiroptera/g' numberoflps_exact.txt
sed -i 's/lagomorpha/Lagomorpha/g' numberoflps_exact.txt
sed -i 's/marsupials/Marsupialia/g' numberoflps_exact.txt
sed -i 's/perissodactyla/Perissodactyla/g' numberoflps_exact.txt
sed -i 's/primates/Primates/g' numberoflps_exact.txt
sed -i 's/rodents/Rodentia/g' numberoflps_exact.txt
sed -i 's/squamata/Squamata/g' numberoflps_exact.txt
sed -i 's/testudines/Testudines/g' numberoflps_exact.txt

################Plotting in R#############
##R
echo 'a=read.table("exactstretch_diversity",header=F)' >> diversity_richness_exact.r
echo 'colnames(a)=c("clades","repeat","count")' >> diversity_richness_exact.r
echo 'library(stringr)' >> diversity_richness_exact.r
echo 'a$clades=str_to_title(a$clades)' >> diversity_richness_exact.r
echo 'finaldf <- data.frame(clades=character(), shd=numeric())' >> diversity_richness_exact.r
echo 'for(c in unique(a$clades)){' >> diversity_richness_exact.r
echo 'b=a[a$clades==c,]' >> diversity_richness_exact.r
echo 'b$V4=(b$count/sum(b$count))*(-log(b$count/sum(b$count)))' >> diversity_richness_exact.r
echo 'sd=sum(b$V4)' >> diversity_richness_exact.r
echo 'd=data.frame(c,sd)' >> diversity_richness_exact.r
echo 'finaldf=rbind(finaldf,d)' >> diversity_richness_exact.r
echo '}' >> diversity_richness_exact.r
echo 'mf=data.frame()' >> diversity_richness_exact.r
echo 'for(clade in c("Afrotheria","Primates","Rodentia","Lagomorpha","Chiroptera","Artiodactyla","Perissodactyla","Carnivora","Marsupialia","Testudines","Aves","Squamata","Amphibia")){' >> diversity_richness_exact.r
echo 'm1=finaldf[finaldf$c==clade,]' >> diversity_richness_exact.r
echo 'mf=rbind(mf,m1)' >> diversity_richness_exact.r
echo '}' >> diversity_richness_exact.r
echo 'finaldf=mf' >> diversity_richness_exact.r
echo 'jpeg("shannon_diversity_exactstretch.jpeg",width=23.2,height=13,units="in",res=300)' >> diversity_richness_exact.r
echo 'barplot(finaldf$sd,names=finaldf$c,ylab="Shannon diversity",ylim=c(0,max(finaldf$sd)),xlab="Clades",col="lightcyan",main="Shannon diversity of exact LPSs across different clades")' >> diversity_richness_exact.r
echo 'dev.off()' >> diversity_richness_exact.r
echo 'simpsondf <- data.frame(clades=character(), simd=numeric())' >> diversity_richness_exact.r
echo 'for(c in unique(a$clades)){' >> diversity_richness_exact.r
echo 's=a[a$clades==c,]' >> diversity_richness_exact.r
echo 's$V4=s$count*(s$count - 1)' >> diversity_richness_exact.r
echo 'simpd=round(1-(sum(s$V4)/(sum(s$count)*(sum(s$count)-1))),4)' >> diversity_richness_exact.r
echo 'simp=data.frame(c,simpd)' >> diversity_richness_exact.r
echo 'simpsondf=rbind(simpsondf,simp)}' >> diversity_richness_exact.r
echo 'sf=data.frame()' >> diversity_richness_exact.r
echo 'for(clade in c("Afrotheria","Primates","Rodentia","Lagomorpha","Chiroptera","Artiodactyla","Perissodactyla","Carnivora","Marsupialia","Testudines","Aves","Squamata","Amphibia")){' >> diversity_richness_exact.r
echo 's1=simpsondf[simpsondf$c==clade,]' >> diversity_richness_exact.r
echo 'sf=rbind(sf,s1)' >> diversity_richness_exact.r
echo '}' >> diversity_richness_exact.r
echo 'simpsondf=sf' >> diversity_richness_exact.r
echo 'jpeg("simpson_diversity_exactstretch.jpeg",width=23.2,height=13,units="in",res=300)' >> diversity_richness_exact.r
echo 'barplot(simpsondf$simpd,names=simpsondf$c,ylab="Simpson diversity",ylim=c(0,1),xlab="Clades",col="lightcyan",main="Simpson diversity of exact LPSs across different clades")' >> diversity_richness_exact.r
echo 'dev.off()' >> diversity_richness_exact.r
echo 'unnorrich <- data.frame(clades=character(), richness=numeric())' >> diversity_richness_exact.r
echo 'for(c in unique(a$clades)){' >> diversity_richness_exact.r
echo 'b=a[a$clades==c,]' >> diversity_richness_exact.r
echo 'numrep=dim(b)[1]' >> diversity_richness_exact.r
echo 'ric=data.frame(c,numrep)' >> diversity_richness_exact.r
echo 'unnorrich=rbind(unnorrich,ric)}' >> diversity_richness_exact.r
echo 'unno=data.frame()' >> diversity_richness_exact.r
echo 'for(clade in c("Afrotheria","Primates","Rodentia","Lagomorpha","Chiroptera","Artiodactyla","Perissodactyla","Carnivora","Marsupialia","Testudines","Aves","Squamata","Amphibia")){' >> diversity_richness_exact.r
echo 's1=unnorrich[unnorrich$c==clade,]' >> diversity_richness_exact.r
echo 'unno=rbind(unno,s1)' >> diversity_richness_exact.r
echo '}' >> diversity_richness_exact.r
echo 'unnorrich=unno' >> diversity_richness_exact.r
echo 'jpeg("unnormalized_repeat_richness_exactstretch.jpeg",width=23.2,height=13,units="in",res=300)' >> diversity_richness_exact.r
echo 'barplot(unnorrich$numrep,names=unnorrich$c,ylab="Repeat richness",xlab="Clades",ylim=c(0,max(unnorrich$numrep)+2),col="#69b3a2",main="Repeats richness for exact LPSs across different clades")' >> diversity_richness_exact.r
echo 'dev.off()' >> diversity_richness_exact.r
echo 'numlps=read.table("numberoflps_exact.txt",header=F)' >> diversity_richness_exact.r
echo 'numlps$V1=str_to_title(numlps$V1)' >> diversity_richness_exact.r
echo 'norrich <- data.frame(clades=character(), richness=numeric())' >> diversity_richness_exact.r
echo 'for(c in unique(a$clades)){' >> diversity_richness_exact.r
echo 'rich=unnorrich[unnorrich$c==c,2]' >> diversity_richness_exact.r
echo 'numrep=numlps[numlps$V1==c,2]' >> diversity_richness_exact.r
echo 'slope=log(rich)/log(numrep)' >> diversity_richness_exact.r
echo 'rich=data.frame(c,slope)' >> diversity_richness_exact.r
echo 'norrich=rbind(norrich,rich)}' >> diversity_richness_exact.r
echo 'nono=data.frame()' >> diversity_richness_exact.r
echo 'for(clade in c("Afrotheria","Primates","Rodentia","Lagomorpha","Chiroptera","Artiodactyla","Perissodactyla","Carnivora","Marsupialia","Testudines","Aves","Squamata","Amphibia")){' >> diversity_richness_exact.r
echo 's1=norrich[norrich$c==clade,]' >> diversity_richness_exact.r
echo 'nono=rbind(nono,s1)' >> diversity_richness_exact.r
echo '}' >> diversity_richness_exact.r
echo 'norrich=nono' >> diversity_richness_exact.r
echo 'jpeg("normalized_repeat_richness_exactstretch.jpeg",width=23.2,height=13,units="in",res=300)' >> diversity_richness_exact.r
echo 'barplot(norrich$slope,names=norrich$c,ylab="Repeat richness",xlab="Clades",ylim=c(0,max(norrich$slope)+0.2),col="#69b3a2",main="Normalized repeats richness for exact LPSs across different clades")' >> diversity_richness_exact.r
echo 'dev.off()' >> diversity_richness_exact.r
echo 'norrich <- data.frame(clades=character(), richness=numeric(),numrepeats=numeric())' >> diversity_richness_exact.r
echo 'for(c in unique(a$clades)){' >> diversity_richness_exact.r
echo 'rich=unnorrich[unnorrich$c==c,2]' >> diversity_richness_exact.r
echo 'numrep=numlps[numlps$V1==c,2]' >> diversity_richness_exact.r
echo 'rich=data.frame(c,rich,numrep)' >> diversity_richness_exact.r
echo 'norrich=rbind(norrich,rich)}' >> diversity_richness_exact.r
echo 'nono=data.frame()' >> diversity_richness_exact.r
echo 'for(clade in c("Afrotheria","Primates","Rodentia","Lagomorpha","Chiroptera","Artiodactyla","Perissodactyla","Carnivora","Marsupialia","Testudines","Aves","Squamata","Amphibia")){' >> diversity_richness_exact.r
echo 's1=norrich[norrich$c==clade,]' >> diversity_richness_exact.r
echo 'nono=rbind(nono,s1)' >> diversity_richness_exact.r
echo '}' >> diversity_richness_exact.r
echo 'norrich=nono' >> diversity_richness_exact.r
echo 'jpeg("log_transformed_richness_exactstretch.jpeg",width=23.2,height=13,units="in",res=300)' >> diversity_richness_exact.r
echo 'b=lm(formula = log(norrich$rich) ~ log(norrich$numrep), data = norrich)' >> diversity_richness_exact.r
echo 'plot(log(norrich$numrep),log(norrich$rich),ylab="log(richness)",xlab="log(number of sequences)",col="blue3",pch=20,main="log transformed richness of exact repeats and number of sequences across different clades")' >> diversity_richness_exact.r
echo 'abline(b,col="red",lty=2)' >> diversity_richness_exact.r
echo 'dev.off()' >> diversity_richness_exact.r
Rscript diversity_richness_exact.r
