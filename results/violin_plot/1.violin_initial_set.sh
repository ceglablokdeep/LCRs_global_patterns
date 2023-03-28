###NOTE: We have provided the file "all_flps_output_with_gene_clade_name" in "initial_dataset/4.flps_related" folder

rm -rf flps_all_output_with_xm_gene_sp_clade
while read j
do
echo $j
clade=`echo $j|awk '{print $1}'`
gene=`echo $j|awk '{print $2}'`
xp=`echo $j|awk '{print $3}'`
sp=`echo $j|awk '{print $4}'`
xm=`grep "$xp" all_datasets/"$gene"_dataset/data/data_table.tsv|cut -f11`
echo $j|awk -v xm=$xm '{print $0,xm}' OFS="\t" >> flps_all_output_with_xm_gene_sp_clade
done < all_flps_output_with_gene_clade_name

###arranging the file phylogenetically###
sed -i 's/afrotheria/Afrotheria/g' flps_all_output_with_xm_gene_sp_clade
sed -i 's/amphibia/Amphibia/g' flps_all_output_with_xm_gene_sp_clade
sed -i 's/artiodactyla/Artiodactyla/g' flps_all_output_with_xm_gene_sp_clade
sed -i 's/aves/Aves/g' flps_all_output_with_xm_gene_sp_clade
sed -i 's/carnivore/Carnivora/g' flps_all_output_with_xm_gene_sp_clade
sed -i 's/chiroptera/Chiroptera/g' flps_all_output_with_xm_gene_sp_clade
sed -i 's/lagomorpha/Lagomorpha/g' flps_all_output_with_xm_gene_sp_clade
sed -i 's/marsupials/Marsupialia/g' flps_all_output_with_xm_gene_sp_clade
sed -i 's/perissodactyla/Perissodactyla/g' flps_all_output_with_xm_gene_sp_clade
sed -i 's/primates/Primates/g' flps_all_output_with_xm_gene_sp_clade
sed -i 's/rodents/Rodentia/g' flps_all_output_with_xm_gene_sp_clade
sed -i 's/squamata/Squamata/g' flps_all_output_with_xm_gene_sp_clade
sed -i 's/testudines/Testudines/g' flps_all_output_with_xm_gene_sp_clade
####
for clade in Afrotheria Primates Rodentia Lagomorpha Chiroptera Artiodactyla Perissodactyla Carnivora Marsupialia Testudines Aves Squamata Amphibia
do
awk -v c=$clade '$1==c{print $0}' flps_all_output_with_xm_gene_sp_clade OFS="\t" >> flps_all_output_with_xm_gene_sp_clade_phylo_arranged
done
