###The output of fLPS2.0 are having extension "flpsout_withcomposition", we have kept all the output files in a directory
##Below, we are listing all the output of fLPS2.0 to a file "all_flps_out_files"
ls|grep "flpsout_withcomposition" > all_flps_out_files
for clades in afrotheria amphibia artiodactyla aves carnivore chiroptera lagomorpha marsupials perissodactyla primates rodents squamata testudines
do
grep "$clades" all_flps_out_files > "$clades"_flps_out
done

##XP_003128180.1_Sus_scrofa       500     SINGLE  4       7       61      16      8.822e-04       {L}     ----    2.234
for clades in afrotheria amphibia artiodactyla aves carnivore chiroptera lagomorpha marsupials perissodactyla primates rodents squamata testudines
do
echo $clades
rm -rf "$clades"_all_flps_combined
while read j
do
awk '{print $1,$2,$5,$6,$7,$9}' "$j" >> "$clades"_all_flps_combined
done < "$clades"_flps_out
done
cat *_all_flps_combined > all_flps_output_combined
