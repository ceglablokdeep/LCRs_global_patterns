##I will use "inputfilenames" file to make individual control file for each aln file
while read t
do
name=`echo $t|sed 's/.nwk//g'`
cp demo.ctl "$name".ctl
alnfile=`ls "$name"_similar_length_orf.aln_oneliner`
outfile=`echo "$name".out`
sed -i "s/ssssss/$alnfile/g" "$name".ctl
sed -i "s/tttttt/$t/g" "$name".ctl
sed -i "s/oooooo/$outfile/g" "$name".ctl
done < inputfilenames

ls|grep ".nwk"|sed 's/.nwk//g' > allcontrolfiles
ls|grep ".nwk"|sed 's/.nwk//g' > allcontrolfiles_backup

for i in `cat allcontrolfiles`
do
codeml $i
done
