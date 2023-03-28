sed -i 's/afrotheria/Afrotheria/g' all_omega_compiled
sed -i 's/amphibia/Amphibia/g' all_omega_compiled
sed -i 's/artiodactyla/Artiodactyla/g' all_omega_compiled
sed -i 's/aves/Aves/g' all_omega_compiled
sed -i 's/carnivore/Carnivora/g' all_omega_compiled
sed -i 's/chiroptera/Chiroptera/g' all_omega_compiled
sed -i 's/lagomorpha/Lagomorpha/g' all_omega_compiled
sed -i 's/marsupials/Marsupialia/g' all_omega_compiled
sed -i 's/perissodactyla/Perissodactyla/g' all_omega_compiled
sed -i 's/primates/Primates/g' all_omega_compiled
sed -i 's/rodents/Rodentia/g' all_omega_compiled
sed -i 's/squamata/Squamata/g' all_omega_compiled
sed -i 's/testudines/Testudines/g' all_omega_compiled

sed -i 's/afrotheria/Afrotheria/g' sorted_all_gene_pro_ac
sed -i 's/amphibia/Amphibia/g' sorted_all_gene_pro_ac
sed -i 's/artiodactyla/Artiodactyla/g' sorted_all_gene_pro_ac
sed -i 's/aves/Aves/g' sorted_all_gene_pro_ac
sed -i 's/carnivore/Carnivora/g' sorted_all_gene_pro_ac
sed -i 's/chiroptera/Chiroptera/g' sorted_all_gene_pro_ac
sed -i 's/lagomorpha/Lagomorpha/g' sorted_all_gene_pro_ac
sed -i 's/marsupials/Marsupialia/g' sorted_all_gene_pro_ac
sed -i 's/perissodactyla/Perissodactyla/g' sorted_all_gene_pro_ac
sed -i 's/primates/Primates/g' sorted_all_gene_pro_ac
sed -i 's/rodents/Rodentia/g' sorted_all_gene_pro_ac
sed -i 's/squamata/Squamata/g' sorted_all_gene_pro_ac
sed -i 's/testudines/Testudines/g' sorted_all_gene_pro_ac
###The above file contains entries from only those clade and genes which are aligned######################
##2904605 sorted_all_gene_pro_ac


sed -i 's/_/\t/1' all_sp_gene_clade_aligned

sed -i 's/afrotheria/Afrotheria/g' all_sp_gene_clade_aligned
sed -i 's/amphibia/Amphibia/g' all_sp_gene_clade_aligned
sed -i 's/artiodactyla/Artiodactyla/g' all_sp_gene_clade_aligned
sed -i 's/aves/Aves/g' all_sp_gene_clade_aligned
sed -i 's/carnivore/Carnivora/g' all_sp_gene_clade_aligned
sed -i 's/chiroptera/Chiroptera/g' all_sp_gene_clade_aligned
sed -i 's/lagomorpha/Lagomorpha/g' all_sp_gene_clade_aligned
sed -i 's/marsupials/Marsupialia/g' all_sp_gene_clade_aligned
sed -i 's/perissodactyla/Perissodactyla/g' all_sp_gene_clade_aligned
sed -i 's/primates/Primates/g' all_sp_gene_clade_aligned
sed -i 's/rodents/Rodentia/g' all_sp_gene_clade_aligned
sed -i 's/squamata/Squamata/g' all_sp_gene_clade_aligned
sed -i 's/testudines/Testudines/g' all_sp_gene_clade_aligned
sed -i 's/ /\t/g' all_sp_gene_clade_aligned




awk '{print $2"_"$1"_"$3,$4,$5}' OFS="\t" sorted_all_gene_pro_ac > sorted_all_gene_pro_ac_columns_combined
awk '{print $1"_"$2"_"$3}' all_sp_gene_clade_aligned > all_sp_gene_clade_aligned_combined

awk -F"\t" 'NR==FNR {a[$1]; next} $1 in a' all_sp_gene_clade_aligned_combined sorted_all_gene_pro_ac_columns_combined > sorted_all_gene_pro_ac_columns_combined_aligned
###selecting only those XP IDs which are aligned
awk '{print $8"_"$1"_"$2,$3,$4,$5,$6,$7}' OFS="\t" all_clades_flps_out_combined > all_clades_flps_out_columns_combined
awk '{print $1"_"$4"_"$3}' OFS="\t" sorted_all_gene_pro_ac > sorted_all_gene_pro_ac_columns_combined
awk -F"\t" 'NR==FNR {a[$1]; next} $1 in a' sorted_all_gene_pro_ac_columns_combined all_clades_flps_out_columns_combined > all_clades_flps_out_columns_combined_aligned
##18704333 all_clades_flps_out_columns_combined_aligned
awk '{print $1"_"$4"_"$3,$2}' OFS="\t" sorted_all_gene_pro_ac > sorted_all_gene_pro_ac_columns_combined_gene_name

awk -F"\t" 'NR == FNR { x[$1]=$0; next}($1 in x) { print x[$1],$0}' OFS="\t" sorted_all_gene_pro_ac_columns_combined_gene_name all_clades_flps_out_columns_combined_aligned|awk '{print $1,$2,$4,$5,$6,$7,$8}' OFS="\t" > flps_output_aligned_with_gene_name
sed -i 's/_/\t/1' flps_output_aligned_with_gene_name
sed -i 's/_/\t/2' flps_output_aligned_with_gene_name
##18704333 flps_output_aligned_with_gene_name

##Afrotheria	XP_006877394.1	Chrysochloris_asiatica	A1BG	519	395	423	9	P
##subsetting flps output according to purity
for pur in `seq 0 0.1 1`
do
awk -v p=$pur '($8/(($7-$6)+1))>=p{print $0}' OFS="\t" flps_output_aligned_with_gene_name > flps_output_morethan_"$pur"_aligned_with_gene_name
awk -v p=$pur '($8/(($7-$6)+1))<p{print $0}' OFS="\t" flps_output_aligned_with_gene_name > flps_output_lessthan_"$pur"_aligned_with_gene_name
done

###
awk '{print $1"_"$2"_"$3,$4}' OFS="\t" all_omega_compiled > omega_columns_combined

for file in `ls *than_*_aligned_with_gene_name|grep -v "_0.0_"`
do
echo $file
awk '{print $1"_"$4"_"$3,$2,$5,$6,$7,$8,$9}' OFS="\t" "$file" > "$file"_column_combined
awk -F"\t" 'NR == FNR { x[$1]=$0; next}($1 in x) { print x[$1],$0}' OFS="\t" omega_columns_combined "$file"_column_combined|awk '{print $1,$2,$4,$5,$6,$7,$8,$9}' OFS="\t" > "$file"_omega_values
done
########################VERY IMPORTANT####################
########741206 flps_output_lessthan_0.1_aligned_with_gene_name_column_combined
########733786 flps_output_lessthan_0.1_aligned_with_gene_name_omega_values
########Difference of sequences = 3901; Those are the sequences with NNNs in between sequence and are excluded from PAML analyses##########

####Since same omega value can repeat multiple times depending on the number of LPSs in the gene for the species, we need to select the unique value to not bias the distribution of Omega values
for file in `ls *_aligned_with_gene_name_omega_values`
do
echo $file
awk '{print $1,$2}' OFS="\t" $file|sort -u|sed 's/_/\t/1' > "$file"_unique_values
done
