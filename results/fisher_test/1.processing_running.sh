awk '{print $1"_"$2}' coordinatecompiledwithalnlength|sort -u > flps_output_unique
awk '{print $1"_"$2}' sitemodel.outresults_together.txt|sort -u > pss_result_unique
sed -i 's/afrotheria/Afrotheria/g' pss_result_unique
sed -i 's/amphibia/Amphibia/g' pss_result_unique
sed -i 's/artiodactyla/Artiodactyla/g' pss_result_unique
sed -i 's/aves/Aves/g' pss_result_unique
sed -i 's/carnivore/Carnivora/g' pss_result_unique
sed -i 's/chiroptera/Chiroptera/g' pss_result_unique
sed -i 's/lagomorpha/Lagomorpha/g' pss_result_unique
sed -i 's/marsupials/Marsupialia/g' pss_result_unique
sed -i 's/perissodactyla/Perissodactyla/g' pss_result_unique
sed -i 's/primates/Primates/g' pss_result_unique
sed -i 's/rodents/Rodentia/g' pss_result_unique
sed -i 's/squamata/Squamata/g' pss_result_unique
sed -i 's/testudines/Testudines/g' pss_result_unique
###as a precaution, we are running a command to uppercase the first character of each file
for file in *_result_unique
do
sed -i "s/\b\(.\)/\u\1/1" $file
done



sed -i 's/afrotheria/Afrotheria/g' flps_output_unique
sed -i 's/amphibia/Amphibia/g' flps_output_unique
sed -i 's/artiodactyla/Artiodactyla/g' flps_output_unique
sed -i 's/aves/Aves/g' flps_output_unique
sed -i 's/carnivore/Carnivora/g' flps_output_unique
sed -i 's/chiroptera/Chiroptera/g' flps_output_unique
sed -i 's/lagomorpha/Lagomorpha/g' flps_output_unique
sed -i 's/marsupials/Marsupialia/g' flps_output_unique
sed -i 's/perissodactyla/Perissodactyla/g' flps_output_unique
sed -i 's/primates/Primates/g' flps_output_unique
sed -i 's/rodents/Rodentia/g' flps_output_unique
sed -i 's/squamata/Squamata/g' flps_output_unique
sed -i 's/testudines/Testudines/g' flps_output_unique
for file in *_output_unique
do
sed -i "s/\b\(.\)/\u\1/1" $file
done

##The file "all_sp_gene_clade_aligned" is already provided in the folder "results/omega_comparison", so we are skipping it here.
awk '{print $1}' all_sp_gene_clade_aligned|sort -u|sed 's/_/\t/g'|awk '{print $2"_"$1}' > all_aln_files

sed -i 's/afrotheria/Afrotheria/g' all_aln_files
sed -i 's/amphibia/Amphibia/g' all_aln_files
sed -i 's/artiodactyla/Artiodactyla/g' all_aln_files
sed -i 's/aves/Aves/g' all_aln_files
sed -i 's/carnivore/Carnivora/g' all_aln_files
sed -i 's/chiroptera/Chiroptera/g' all_aln_files
sed -i 's/lagomorpha/Lagomorpha/g' all_aln_files
sed -i 's/marsupials/Marsupialia/g' all_aln_files
sed -i 's/perissodactyla/Perissodactyla/g' all_aln_files
sed -i 's/primates/Primates/g' all_aln_files
sed -i 's/rodents/Rodentia/g' all_aln_files
sed -i 's/squamata/Squamata/g' all_aln_files
sed -i 's/testudines/Testudines/g' all_aln_files

for clade in Afrotheria Amphibia Artiodactyla Aves Carnivora Chiroptera Marsupialia Perissodactyla Primates Rodentia Squamata Testudines
do
echo $clade
echo -e "#pss_genes\t#no_pss" > "$clade"_contingency
cmn=`paste <(grep $clade pss_result_unique|sort -u) <(grep $clade flps_output_unique|sort -u) -d "\n"|sort|uniq -c|awk '$1==2{print $0}'|wc -l`
flps_num=`grep -c "$clade" flps_output_unique`
rpt_only=`echo $flps_num $cmn|awk '{print $1-$2}'`
echo -e "$cmn\t$rpt_only" >> "$clade"_contingency
pss_only=`grep -c "$clade" pss_result_unique|awk -v c=$cmn '{print $1-c}'`
neither=`grep -c "$clade" all_aln_files|awk -v c=$cmn -v r=$rpt_only -v p=$pss_only '{print $1-(c+r+p)}'`
echo -e "$pss_only\t$neither" >> "$clade"_contingency
done

##Making the header of the result file
echo -e "Clade\tfisher_odds_ratio\tpval(fisher_greater)\tpval(fisher_lesser)\tpval(fisher_twoside)" > fisher_contingency_out.txt

Rscript fisher_contingency.r
