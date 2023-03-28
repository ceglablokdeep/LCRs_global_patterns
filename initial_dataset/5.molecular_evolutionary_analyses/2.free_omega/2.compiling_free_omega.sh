###compiling the results of PAML free ratio on khorana
ls|grep "out" > alloutfiles
split -d --additional-suffix=_split.txt -l 35 alloutfiles

for files in `ls x*_split.txt`
do
for j in `cat $files`
do
(echo $j
name=`echo $j|sed 's/.out//g'`
grep -A1 "w ratios as labels for TreeView" "$j"|tail -n1|sed -e 's/(/\n/g' -e 's/)/\n/g' -e 's/,/\n/g' -e 's/;/\n/g' -e 's/ //g'|grep "^[A-Z]"|sed 's/#/\t/g' > "$name"_paml_perspecies ) &
done
wait
done
rm x*split*
ls|grep "paml_perspecies" > allpamlomegafiles
mkdir ../pamlomegafiles
split -d --additional-suffix=_split.txt -l 35 allpamlomegafiles
for files in `ls x*_split.txt`
do
for j in `cat $files`
do
(echo $j
cp "$j" ../pamlomegafiles ) &
done
wait
done
#to split the omega files according to "with repeats" and "without repeats", we will make use of coordinatecompiledwithalnlength file

ls|grep "_paml_perspecies" > allpamlfiles
split -d --additional-suffix=_split.txt -l 60 allpamlfiles

for files in `ls x*_split.txt`
do
for j in `cat $files`
do
(echo $j
col1=`echo $j|sed 's/_paml_perspecies//g'`
clade=`echo $col1|cut -f2 -d_`
gene=`echo $col1|cut -f1 -d_`
presabs=`awk -v c=$clade -v g=$gene '$1==c && $2==g{print $0}' coordinatecompiledwithalnlength|wc -l`
if [ $presabs -eq 0 ]
then
cat "$j" >> "$clade"_pamlomega_without_repeats
else
cat "$j" >> "$clade"_pamlomega_with_repeats
fi ) &
done
wait
done
#################compiling the data of paml omega
for clade in afrotheria amphibia artiodactyla aves carnivore chiroptera marsupials perissodactyla primates rodents squamata testudines
do
echo $clade
cl=`echo "${clade^}"`
awk -v c=$cl '{print c,$2,"With_repeats"}' OFS="\t" "$clade"_pamlomega_with_repeats >> omegafree_compiled_allclades
awk -v c=$cl '{print c,$2,"Without_repeats"}' OFS="\t" "$clade"_pamlomega_without_repeats >> omegafree_compiled_allclades
done
###correct the names of the clades
sed -i 's/Carnivore/Carnivora/g' omegafree_compiled_allclades
sed -i 's/Rodents/Rodentia/g' omegafree_compiled_allclades
