####making stacked barplot of top 20 repeat across clades####
##a. making proportion file
rm -f colone.txt
for clades in afrotheria amphibia artiodactyla aves carnivore chiroptera lagomorpha marsupials perissodactyla primates rodents squamata testudines
do
sumcol=`grep "\b$clades\b" morethanseventy_diversity|awk '{sum1+=$3} END{print sum1}'`
grep -i "\b$clades\b" morethanseventy_diversity|sort -k3nr|awk  -v s=$sumcol '{print $0,$3/s}' OFS='\t' > "$clades"_morethanseventy_proportion.txt
head -n20 "$clades"_morethanseventy_proportion.txt|awk '{print $2}' >> colone.txt
done
sort -u colone.txt > temp.txt
mv temp.txt colone.txt
##b. getting proportion of these aa LPSs from each clade
echo -e "repeats\tafrotheria\tamphibia\tartiodactyla\taves\tcarnivore\tchiroptera\tlagomorpha\tmarsupials\tperissodactyla\tprimates\trodents\tsquamata\ttestudines" > top_twenty_proportion_morethanseventy.txt
for repeats in `cat colone.txt`
do
echo "$repeats" > abc.txt
for clades in afrotheria amphibia artiodactyla aves carnivore chiroptera lagomorpha marsupials perissodactyla primates rodents squamata testudines
do
files=`ls "$clades"_morethanseventy_proportion.txt`
var1=`grep -c "\b$repeats\b" $files`
if [ $var1 -gt 0 ]
then
num=`grep "\b$repeats\b" $files|awk '{print $4}'`
else
num=0
fi
echo "$num" >> abc.txt
done
cat abc.txt|xargs|sed 's/ /\t/g' >> top_twenty_proportion_morethanseventy.txt
rm abc.txt
done
rm colone.txt
###c. adding remaining proportion to the file
other=`awk '{for (i=2;i<=NF;i++) sum[i]+=$i;}; END{for (i in sum) print 1-sum[i];}' top_twenty_proportion_morethanseventy.txt|xargs|sed 's/ /\t/g'`
echo -e "other\t$other">> top_twenty_proportion_morethanseventy.txt
##d. plotting in R
#R
echo 'm=read.table("top_twenty_proportion_morethanseventy.txt",header=T)' > top_twenty_proportion.r
echo 'library(stringr)' >> top_twenty_proportion.r
echo 'colnames(m)=str_to_title(colnames(m))' >> top_twenty_proportion.r
echo 'rownames(m)=m[,1]' >> top_twenty_proportion.r
echo 'm=m[,-1]' >> top_twenty_proportion.r
echo 'col21rainbow= c("#771155", "#AA4488", "#CC99BB", "#114477", "#4477AA", "#77AADD", "#117777", "#44AAAA", "#77CCCC", "#117744", "#44AA77", "#88CCAA", "#777711", "#AAAA44", "#DDDD77", "#774411", "#AA7744", "#DDAA77", "#771122", "#AA4455", "#DD7788","orange","cyan","black","violet","red","darkblue","yellow","thistle1","greenyellow","wheat1")' >> top_twenty_proportion.r
echo 'jpeg("top20_morethanseventy_lps_proportion.jpeg",width=23.2,height=13,units="in",res=300)' >> top_twenty_proportion.r
echo 'par(mai=c(1.25,3,1,2))' >> top_twenty_proportion.r
echo 'barplot(as.matrix(m),col=col21rainbow,main="Proportion of top 20 >70% LPSs across clades",cex.main=2,xlab="Proportion of LPSs",ylab="",horiz=T,las=2,cex.axis=0.5,cex.lab=2)' >> top_twenty_proportion.r
echo 'legend("topright",inset=c(-0.08,0.1),legend=row.names(m),col=col21rainbow,fill=col21rainbow,bty="n",xpd=T)' >> top_twenty_proportion.r
echo 'title(ylab="Clades",line=6.5,cex.lab=2)' >> top_twenty_proportion.r
echo 'dev.off()' >> top_twenty_proportion.r
Rscript top_twenty_proportion.r
