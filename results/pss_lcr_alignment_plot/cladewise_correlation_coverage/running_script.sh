######################To check for correlation between positively selected sites, repeat locations and coverage in cladewise manner#########

while read j
do
len=`echo $j|awk '{print $2}'`
nam=`echo $j|awk '{print $1}'|sed 's/_similar_length_orf.aln_oneliner//g'`
echo $nam
echo $j|awk '{$1=$2=$3="";print $0}' OFS="\n"|sed '/^$/d'|awk -v l=$len -v n=$nam '{print n,$1,(NR/l)*100}' >> alndensity_withnormpositions_cladewise.txt
done < alndensity_normalized.txt

mkdir making_averagefile
cp alndensity_withnormpositions_cladewise.txt making_averagefile
cd making_averagefile
for clades in afrotheria amphibia artiodactyla aves carnivore chiroptera marsupials perissodactyla primates rodents squamata testudines
do
grep "$clades" alndensity_withnormpositions_cladewise.txt > "$clades"_aln.txt
for i in `seq 1 100`
do
(j=`echo $i|awk '{print $1-1}'`
echo $clades $i $j
awk -v llim=$j -v ulim=$i '$3>=llim && $3<ulim{print $2,$3}' OFS="\t" "$clades"_aln.txt > temp_"$i".txt
awk -v inter=$i 'END{y=NR} {x+=$1} END{print inter-0.5,x/y}' OFS="\t" temp_"$i".txt >> "$clades"_withinterval_"$i".txt)&
done
wait
done
rm temp_*.txt
for clades in afrotheria amphibia artiodactyla aves carnivore chiroptera marsupials perissodactyla primates rodents squamata testudines
do
touch "$clades"_final_alndensity.txt
for i in `seq 1 100`
do
cat "$clades"_final_alndensity.txt "$clades"_withinterval_"$i".txt >> "$clades"_final_alndensity.txt
done
done

mkdir interval_files
mv *_withinterval_*.txt interval_files

for clades in afrotheria amphibia artiodactyla aves carnivore chiroptera marsupials perissodactyla primates rodents squamata testudines
do
echo $clades
grep "$clades" sitenorm.txt > cladesite.txt
grep "$clades" rpt_norm.txt > claderpt.txt
Rscript alndensity.r cladesite.txt claderpt.txt "$clades"_final_alndensity.txt "$clades"
done
###############alndensity cladewise ends here##########################
