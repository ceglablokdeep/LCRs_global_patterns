###This script is to test the overlap of PSSs and LCRs in a clade-wise manner using bedtools fisher.
##Correcting the names of the clades in all the initial files
sed -i 's/afrotheria/Afrotheria/g' sitemodel.outresults_together.txt
sed -i 's/amphibia/Amphibia/g' sitemodel.outresults_together.txt
sed -i 's/artiodactyla/Artiodactyla/g' sitemodel.outresults_together.txt
sed -i 's/aves/Aves/g' sitemodel.outresults_together.txt
sed -i 's/carnivore/Carnivora/g' sitemodel.outresults_together.txt
sed -i 's/chiroptera/Chiroptera/g' sitemodel.outresults_together.txt
sed -i 's/lagomorpha/Lagomorpha/g' sitemodel.outresults_together.txt
sed -i 's/marsupials/Marsupialia/g' sitemodel.outresults_together.txt
sed -i 's/perissodactyla/Perissodactyla/g' sitemodel.outresults_together.txt
sed -i 's/primates/Primates/g' sitemodel.outresults_together.txt
sed -i 's/rodents/Rodentia/g' sitemodel.outresults_together.txt
sed -i 's/squamata/Squamata/g' sitemodel.outresults_together.txt
sed -i 's/testudines/Testudines/g' sitemodel.outresults_together.txt

sed -i 's/afrotheria/Afrotheria/g' coordinatecompiledwithalnlength
sed -i 's/amphibia/Amphibia/g' coordinatecompiledwithalnlength
sed -i 's/artiodactyla/Artiodactyla/g' coordinatecompiledwithalnlength
sed -i 's/aves/Aves/g' coordinatecompiledwithalnlength
sed -i 's/carnivora/Carnivora/g' coordinatecompiledwithalnlength
sed -i 's/chiroptera/Chiroptera/g' coordinatecompiledwithalnlength
sed -i 's/marsupials/Marsupialia/g' coordinatecompiledwithalnlength
sed -i 's/perissodactyla/Perissodactyla/g' coordinatecompiledwithalnlength
sed -i 's/primates/Primates/g' coordinatecompiledwithalnlength
sed -i 's/rodentia/Rodentia/g' coordinatecompiledwithalnlength
sed -i 's/squamata/Squamata/g' coordinatecompiledwithalnlength
sed -i 's/testudines/Testudines/g' coordinatecompiledwithalnlength

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

sed -i 's/afrotheria/Afrotheria/g' mergedfile
sed -i 's/amphibia/Amphibia/g' mergedfile
sed -i 's/artiodactyla/Artiodactyla/g' mergedfile
sed -i 's/aves/Aves/g' mergedfile
sed -i 's/carnivore/Carnivora/g' mergedfile
sed -i 's/chiroptera/Chiroptera/g' mergedfile
sed -i 's/lagomorpha/Lagomorpha/g' mergedfile
sed -i 's/marsupials/Marsupialia/g' mergedfile
sed -i 's/perissodactyla/Perissodactyla/g' mergedfile
sed -i 's/primates/Primates/g' mergedfile
sed -i 's/rodents/Rodentia/g' mergedfile
sed -i 's/squamata/Squamata/g' mergedfile
sed -i 's/testudines/Testudines/g' mergedfile
sed -i 's/carnivora/Carnivora/g' mergedfile
sed -i 's/rodentia/Rodentia/g' mergedfile

####the output of PAML sitemodel.outresults_together.txt has amino acid coordinates, we need to convert it to codon aligned coordinates
awk '{print $1"_"$2,$3*3,($3*3)+2}' OFS='\t' sitemodel.outresults_together.txt > sitemodelcoordinates
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

for clade in Afrotheria Amphibia Artiodactyla Aves Carnivora Chiroptera Marsupialia Perissodactyla Primates Rodentia Squamata Testudines
do
echo $clade
grep "$clade" common > "$clade"_common
while read j
do
awk -v gn=$j '$1==gn{print $0}' OFS="\t" mergedfile >> "$clade"_rpt_merged_common
awk -v gn=$j '$1==gn{print $0}' OFS="\t" sitemodelcoordinates >> "$clade"_pss_common
awk -v gn=$j '$1==gn{print $0}' OFS="\t" lengthfile.bed >> "$clade"_length_common
done < "$clade"_common
res=`bedtools fisher -a "$clade"_rpt_merged_common -b "$clade"_pss_common -g "$clade"_length_common|tail -n1`
echo -e "$clade\t$res" >> whole_clade_fisher_result
done


