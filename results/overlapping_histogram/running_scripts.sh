Rscript overlap_histo_pss_rpt.r

chmod a+x overlap_histo.r

cp sitenorm.txt sitenorm_clades.txt
sed -i 's/_/\t/1' sitenorm_clades.txt
sed -i 's/ /\t/g' sitenorm_clades.txt
for clades in afrotheria amphibia artiodactyla aves carnivore chiroptera lagomorpha marsupials perissodactyla primates rodents squamata testudines
do
Rscript overlap_histo.r "$clades"_all_flps_combined "$clades"
done

rm -rf *_flps_morethanseventy_repeats
for files in `ls *_all_flps_combined`
do
echo $files
clade=`echo $files|sed 's/_all_flps_combined//g'`
awk '(($4-$3)+1)==$5{print $0}' "$files" > "$clade"_flps_pure_repeats
awk '($5/(($4-$3)+1))>0.7{print $0}' "$files" > "$clade"_flps_morethanseventy_repeats
done


chmod a+x overlap_histomorethanseventy.r
for clades in afrotheria amphibia artiodactyla aves carnivore chiroptera lagomorpha marsupials perissodactyla primates rodents squamata testudines
do
Rscript overlap_histomorethanseventy.r "$clades"_flps_morethanseventy_repeats "$clades"
done


chmod a+x overlap_histopure.r
for clades in afrotheria amphibia artiodactyla aves carnivore chiroptera lagomorpha marsupials perissodactyla primates rodents squamata testudines
do
Rscript overlap_histopure.r "$clades"_flps_pure_repeats "$clades"
done


chmod a+x overlap_histosimpleseqs.r
for clades in afrotheria amphibia artiodactyla aves carnivore chiroptera lagomorpha marsupials perissodactyla primates rodents squamata testudines
do
Rscript overlap_histosimpleseqs.r "$clades"
done


chmod a+x overlap_histo_all_clades_rpt_grad.r
Rscript overlap_histo_all_clades_rpt_grad.r


chmod a+x overlapping_histo_pss_rpt_all_genes.r
Rscript overlapping_histo_pss_rpt_all_genes.r
