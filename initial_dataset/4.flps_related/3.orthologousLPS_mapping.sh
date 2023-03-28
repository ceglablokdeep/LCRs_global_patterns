## making the coordinate files for the respective clade having more than one species for a repeat and saving the rest in "Discarded.txt" file.
rm -f *_aa_coordinates.txt
for clade in `cut -f1 morethanseventy_combined|sort -u`
do
for gene in `cut -f2 "$clade"_morethanseventy_combined|sort -u`
do
for repeat in `awk -v g=$gene '$2==g{print $0}' "$clade"_morethanseventy_combined|awk '{print $11}'|sort -u`
do
awk -v g=$gene '$2==g{print $0}' "$clade"_morethanseventy_combined|awk -v r=$repeat '$11=r{print $3}' > species.list
numsp=`wc -l species.list|awk '{print $1}'`
if [ $numsp -gt 1 ]
then
awk -v g=$gene -v r=$repeat '$2==g && $11==r{print $0}' "$clade"_morethanseventy_combined|awk '{print $1,$2,$3,$7,$8,$11}' OFS="\t" >> "$gene"_"$clade"_"$repeat"_aa_coordinates.txt
echo "$gene $clade $numsp $repeat"
else
awk -v g=$gene -v r=$repeat '$2==g && $11==r{print $0}' "$clade"_morethanseventy_combined|awk '{print $1,$2,$3,$7,$8,$11}' OFS="\t" >> Discarded.txt
fi
done
done
done
ls|grep "_aa_coordinates.txt" > allcoordinatefiles
mkdir ../allcoordinateoutputfiles
while read j
do
echo $j
mv "$j" ../allcoordinateoutputfiles
done < allcoordinatefiles
cd ../allcoordinateoutputfiles
find . -name '*_aa_coordinates.txt' -type f | xargs wc -l|awk '$1>3'|cut -f2 -d"/"|grep "_aa_coordinates.txt" > allfileswithmorethanthree
mkdir morethanthreespecies
cut -f1,2 -d"_" allfileswithmorethanthree|sort -u > foralnfiles
###moving aln files to same folder that have more than three species
while read j
do
file=`echo "$j"_similar_length_orf.aln_oneliner`
cp /path/to/alnfiles_oneliner/"$file" morethanthreespecies
done < foralnfiles
###moving coordinate files with more than three lines in them to morethanthreespecies folder
while read j
do
cp "$j" morethanthreespecies
done < allfileswithmorethanthree
cd morethanthreespecies
####since all repeat files with more than three species does not have alignment (probably because some sequences had premature stop codon in between and after excluding them, the number of species in MSA files were less than four), now we are keeping aln files and repeats files in same folder where both are available
mkdir alnandrepeatcommon
ls|grep aa_coordinates.txt|cut -f1,2 -d_|sort -u > coordfiles
ls|grep aln_oneliner|cut -f1,2 -d_|sort -u > alnfiles
cat coordfiles alnfiles|sort|uniq -c|awk '$1==2{print $0}'
cat coordfiles alnfiles|sort|uniq -c|awk '$1==2{print $2}' > bothfilespresent
while read j
do
echo $j
cp "$j"* alnandrepeatcommon
done < bothfilespresent
cd alnandrepeatcommon
ls|grep "_aa_coordinates.txt" > allcoordinatefiles
split -d --additional-suffix=_coordinate_splits -l 8958 allcoordinatefiles
#####we have moved coordinate_mapper.sh here in tis folder. Now we gonna run it.
cp allcoordinatefiles allcoordinatefiles_backup
while [ -s allcoordinatefiles ]
do
j=`ps -ef|grep -c coordinate_mapper.sh`
if [ $j -lt 12 ]
then
i=`head -n1 allcoordinatefiles`
( bash coordinate_mapper.sh "$i" ) &
sed -i '1d' allcoordinatefiles
fi
done

#####concat all the coordinate_mapper.sh output together and save them in coordcompiled file
cat *coordinates_repeats.txt > coordoutcompiled
awk '$5!=0 && $6!=0 && $7!=$8{print $0}' coordoutcompiled > coordinateoutfiltered

while read j
do
echo $j
clade=`echo $j|awk '{print $1}'`
gene=`echo $j|awk '{print $2}'`
species=`echo $j|awk '{print $3}'`
alnlength=`head -n2 alnfiles_oneliner/"$gene"_"$clade"_similar_length_orf.aln_oneliner|tail -n1|awk '{print length}'`
echo "$j $alnlength" >> coordinatecompiledwithalnlength
done < coordinateoutfiltered
