################coordinate_mapper.sh################
for file in `ls $1`
do
echo $file
gene=`echo $file|cut -d '_' -f1`
clade=`echo $file|cut -d '_' -f2`
repeat=`echo $file|cut -d '_' -f3`
alnfile=`echo "$gene"_"$clade"_similar_length_orf.aln_oneliner`
echo $gene $clade $repeat $alnfile
if [ -f $alnfile ]
then
while read oneline
do
species=`echo $oneline|awk '{print $3}'|cut -f3- -d "_"|sort -u`
(sequ=`sed -e '/>/s/$/#/g' -e 's/>/@>/g' "$gene"_"$clade"_similar_length_orf.aln_oneliner| sed -z 's/\n//g'|sed -e 's/@/\n/g' -e 's/#/\n/g' |grep -A1 "$species"|grep -v "^>" |sed -z 's/\n//g'`
j1=`echo $oneline|awk '{print $4}'`
j2=`echo $oneline|awk '{print $5}'`
len=`echo $sequ|awk '{print length}'`
init=`echo $j1|awk '{print ($1*3) - 2}'`
last=`echo $j2|awk '{print $1*3}'`
a1=0
a2=0
b1=0
b2=0
for char in `echo $sequ|sed 's/.\{1\}/&\t/g'`
do
if [ $a1 -eq $init ]
then
break
elif [ $char != '-' ]
then
a1=`echo $a1|awk '{print $1 + 1}'`
fi
a2=`echo $a2|awk '{print $1+1}'`
done
for char in `echo $sequ|sed 's/.\{1\}/&\t/g'`
do
if [ $b1 -eq $last ]
then
break
elif [ $char != '-' ]
then
b1=`echo $b1|awk '{print $1 + 1}'`
fi
b2=`echo $b2|awk '{print $1+1}'`
done
echo -e "$clade\t$gene\t$species\t$repeat\t$a1\t$b1\t$a2\t$b2"
echo -e "$clade\t$gene\t$species\t$repeat\t$a1\t$b1\t$a2\t$b2" >> "$clade"_"$gene"_"$species"_"$repeat"_coordinates_repeats.txt )&
done < $file
else
echo "$file" >> aln_file_not_present.txt
fi
done
