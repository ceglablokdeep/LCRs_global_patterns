cd $curr_dir
mkdir flps_output
cd all_datasets
for dir in `ls -d *_dataset`
do
cp "$dir"/data/*_flpsout_withcomposition ../flps_output
done
cd "$curr_dir"/flps_output
## remove {X} or any repeat with X, select more than 70% and stretch with more than 3aa
############################################################################################################################################################################################
## combine all the composition file to make a single file
rm -v !(*_flpsout_withcomposition)
for files in `ls|grep "flpsout_withcomposition"`
do
gene=`echo $files|cut -f1 -d"_"`
clade=`echo $files|cut -f2 -d"_"`
awk -v g=$gene -v c=$clade '$7>3 && $9 !~ /X/ && ($7/(($6-$5)+1)) > 0.7 {print c,g,$0}' OFS='\t' $files >> morethanseventy_combined
done
#########################
##subsetting the fLPS output cladewise to reduce runtime for each clade
for clades in afrotheria amphibia artiodactyla aves carnivore chiroptera lagomorpha marsupials perissodactyla primates rodents squamata testudines
do
grep "$clades" morethanseventy_combined > "$clades"_morethanseventy_combined
done
#############
####Removing all the LPS with more than 4 different amino acids in the composition
for clades in afrotheria amphibia artiodactyla aves carnivore chiroptera lagomorpha marsupials perissodactyla primates rodents squamata testudines
do
awk 'length($11) < 5 {print $0}' "$clades"_morethanseventy_combined > temp.txt
mv temp.txt "$clades"_morethanseventy_combined
done

###############
#making a file with all unique LPSs
for clades in afrotheria amphibia artiodactyla aves carnivore chiroptera lagomorpha marsupials perissodactyla primates rodents squamata testudines
do
cut -f11 "$clades"_morethanseventy_combined|sort -u > "$clades"_allrepeats
done

######Running the below script will generate permute.wk which identifies LPS of same composition that are written in different combinations, i.e. AGC=GAC=ACG=CAG
##save file in permute.awk
echo 'function permute(s, st,     i, j, n, tmp) {' > permute.awk
echo '    n = split(s, item,//)' >> permute.awk
echo '    if (st > n) {  print s; return }' >> permute.awk
echo '    for (i=st; i<=n; i++) {' >> permute.awk
echo '        if (i != st) {' >> permute.awk
echo '         tmp = item[st]; item[st] = item[i]; item[i] = tmp' >> permute.awk
echo '         nextstr = item[1]' >> permute.awk
echo '         for (j=2; j<=n; j++) nextstr = nextstr delim item[j]' >> permute.awk
echo '        }else {' >> permute.awk
echo '          nextstr = s' >> permute.awk
echo '        }' >> permute.awk
echo '       permute(nextstr, st+1)' >> permute.awk
echo '       n = split(s, item, //)' >> permute.awk
echo '   }' >> permute.awk
echo '}' >> permute.awk
echo '{ permute($0,1) }' >> permute.awk
###########################################
for clades in afrotheria amphibia artiodactyla aves carnivore chiroptera lagomorpha marsupials perissodactyla primates rodents squamata testudines
do
files="$clades"_morethanseventy_combined
while [ -s "$clades"_allrepeats ]
do
a=`head -n1 "$clades"_allrepeats`
sed -i "/\b$a\b/d" "$clades"_allrepeats
echo "$a" > "$clades"_permute_file
for newrep in `awk -f permute.awk "$clades"_permute_file`
do
  echo -e "$a\t$newrep"
  awk -v ol=$a -v k=$newrep '$11==k{$11=ol}1' OFS='\t' $files > "$clades"_temp && mv "$clades"_temp $files
  sed -i "/\b$newrep\b/d" "$clades"_allrepeats
done
done
rm "$clades"_allrepeats "$clades"_permute_file
done

####################################


