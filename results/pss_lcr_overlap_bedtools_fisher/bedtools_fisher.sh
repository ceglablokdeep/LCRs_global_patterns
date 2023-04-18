#####Bedtools fisher to look for genes that significantly favor overlap of PSS and repeats
#The directory contains coordinatecompiledwithalnlength sitemodel.outresults_together.txt files
#firstly, we convert Bos_indicus_x_Bos_taurus to Bos_hybrid otherwise it can effect our downstream analyses
sed -i 's/Bos_indicus_x_Bos_taurus/Bos_hybrid/g' coordinatecompiledwithalnlength
sed -i 's/Bos_indicus_x_Bos_taurus/Bos_hybrid/g' sitemodel.outresults_together.txt
###making file for bedtools merge
awk '{print $1"_"$2,$7,$8}' OFS='\t' coordinatecompiledwithalnlength|sort -u > repeatscoordinatesfiles
sort -u repeatscoordinatesfiles > sortedrepeatfile
###since bedtools merge cannot work for multiple chromosomes/scaffolds in one file, we will use it on one clade_gene at a time and concat them
touch mergedfile
for gc in `cut -f1 sortedrepeatfile|sort -u`
do
grep "\b$gc\b" sortedrepeatfile|sort -k2n > temp.bed
bedtools merge -i temp.bed > out.bed
cat mergedfile out.bed >> mergedfile
done
rm temp.bed 
####the output of PAML sitemodel.outresults_together.txt has amino acid coordinates, we need to convert it to codon aligned coordinates
awk '{print $1"_"$2,$3*3,($3*3)+2}' OFS='\t' sitemodel.outresults_together.txt > sitemodelcoordinates
sed -i 's/rodents/rodentia/g' sitemodelcoordinates
paste <(awk '{print $1}' sitemodelcoordinates|sort -u) <(awk '{print $1}' mergedfile|sort -u) -d "\n"|sort|uniq -c|awk '$1>1{print $2}'|sed '/^$/d' > common

rm -f lengthfile.bed
for com in `cat common`
do
echo "$com"
gene=`echo $com|cut -f2 -d"_"`
clade=`echo $com|cut -f1 -d"_"`
len=`awk -v g=$gene -v cl=$clade '$1==cl && $2==g {print $9}' coordinatecompiledwithalnlength|head -n1`
echo -e "$com\t$len" >> lengthfile.bed
done

#######now doing bedtools fisher####
rm -f fisher_result
while read j
do
grep "$j\b" lengthfile.bed > g.bed
grep "$j\b" mergedfile > b.bed
grep "$j\b" sitemodelcoordinates > a.bed
res=`bedtools fisher -a a.bed -b b.bed -g g.bed|tail -n1`
echo -e "$j\t$res" >> fisher_result
done < common
## Second p-value is about the favor of occuring together.
awk '$3<0.05{print $1,$3}' fisher_result OFS="\t"|sed 's/_/\t/g'|sed 's/ /\t/g' |sed -rn 's/./\u&/1p' > positive_sites_favoured_repeats
##The below command capitalizes the first alphabet
sed -rn 's/./\u&/1p' common > for_gedit

## First p-value tells about avoiding the same position
awk '$2<0.05{print $1,$2}' fisher_result OFS="\t"|sed 's/_/\t/g'|sed 's/ /\t/g' |sed -rn 's/./\u&/1p' > positive_sites_avoiding_repeats

###########The script ends here. Have a good day!##################
