##NOTE: I have already provided the processed files from the raw output
##I am providing the code to provide an insight of the steps included i nthe processing
#Hope you have a nice day ahead!! :)


##Making plots of cladewise repeat density on normalized gene position
mkdir cladewise_rpt_pos_density
cd cladewise_rpt_pos_density
cp ../*_all_flps_combined .
mv afrotheria_all_flps_combined Afrotheria_all_flps_combined
mv amphibia_all_flps_combined Amphibia_all_flps_combined
mv artiodactyla_all_flps_combined Artiodactyla_all_flps_combined
mv aves_all_flps_combined Aves_all_flps_combined
mv carnivore_all_flps_combined Carnivora_all_flps_combined
mv chiroptera_all_flps_combined Chiroptera_all_flps_combined
mv lagomorpha_all_flps_combined Lagomorpha_all_flps_combined
mv marsupials_all_flps_combined Marsupialia_all_flps_combined
mv perissodactyla_all_flps_combined Perissodactyla_all_flps_combined
mv primates_all_flps_combined Primates_all_flps_combined
mv rodents_all_flps_combined Rodentia_all_flps_combined
mv squamata_all_flps_combined Squamata_all_flps_combined
mv testudines_all_flps_combined Testudines_all_flps_combined

sed -i 's/{//g' Afrotheria_all_flps_combined
sed -i 's/{//g' Amphibia_all_flps_combined
sed -i 's/{//g' Artiodactyla_all_flps_combined
sed -i 's/{//g' Aves_all_flps_combined
sed -i 's/{//g' Carnivora_all_flps_combined
sed -i 's/{//g' Chiroptera_all_flps_combined
sed -i 's/{//g' Lagomorpha_all_flps_combined
sed -i 's/{//g' Marsupialia_all_flps_combined
sed -i 's/{//g' Perissodactyla_all_flps_combined
sed -i 's/{//g' Primates_all_flps_combined
sed -i 's/{//g' Rodentia_all_flps_combined
sed -i 's/{//g' Squamata_all_flps_combined
sed -i 's/{//g' Testudines_all_flps_combined
sed -i 's/}//g' Afrotheria_all_flps_combined
sed -i 's/}//g' Amphibia_all_flps_combined
sed -i 's/}//g' Artiodactyla_all_flps_combined
sed -i 's/}//g' Aves_all_flps_combined
sed -i 's/}//g' Carnivora_all_flps_combined
sed -i 's/}//g' Chiroptera_all_flps_combined
sed -i 's/}//g' Lagomorpha_all_flps_combined
sed -i 's/}//g' Marsupialia_all_flps_combined
sed -i 's/}//g' Perissodactyla_all_flps_combined
sed -i 's/}//g' Primates_all_flps_combined
sed -i 's/}//g' Rodentia_all_flps_combined
sed -i 's/}//g' Squamata_all_flps_combined
sed -i 's/}//g' Testudines_all_flps_combined
###Replace LPS of NA with AN otherwise R will read them as "NA" (Not-available)
sed -i 's/NA/AN/g' Afrotheria_all_flps_combined
sed -i 's/NA/AN/g' Amphibia_all_flps_combined
sed -i 's/NA/AN/g' Artiodactyla_all_flps_combined
sed -i 's/NA/AN/g' Aves_all_flps_combined
sed -i 's/NA/AN/g' Carnivora_all_flps_combined
sed -i 's/NA/AN/g' Chiroptera_all_flps_combined
sed -i 's/NA/AN/g' Lagomorpha_all_flps_combined
sed -i 's/NA/AN/g' Marsupialia_all_flps_combined
sed -i 's/NA/AN/g' Perissodactyla_all_flps_combined
sed -i 's/NA/AN/g' Primates_all_flps_combined
sed -i 's/NA/AN/g' Rodentia_all_flps_combined
sed -i 's/NA/AN/g' Squamata_all_flps_combined
sed -i 's/NA/AN/g' Testudines_all_flps_combined

for clade in Afrotheria Amphibia Artiodactyla Aves Carnivora Chiroptera Lagomorpha Marsupialia Perissodactyla Primates Rodentia Squamata Testudines
do
echo $clade
Rscript repeatdensity.r $clade
done


###############To make the plots of AA repeats types across normalized gene position across all clades

for clade in Afrotheria Amphibia Artiodactyla Aves Carnivora Chiroptera Lagomorpha Marsupialia Perissodactyla Primates Rodentia Squamata Testudines
do
echo $clade
awk -v c=$clade '{print $0,c}' OFS="\t" "$clade"_all_flps_combined > "$clade"_flps_out_with_clade
done

cat *_flps_out_with_clade > all_clades_flps_out_combined
##################
##To plot the images run the below command. We have provided the R script named aarptdensity.r to generate the necessary plots.
Rscript aarptdensity.r
