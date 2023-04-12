##NOTE: I have already provided the processed files from the raw output
##I am providing the code to provide an insight of the steps included in the processing
#Hope you have a nice day ahead!! :)
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
