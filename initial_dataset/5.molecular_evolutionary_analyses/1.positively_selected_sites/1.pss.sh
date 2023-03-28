#####IMPORTANT: Please make sure that you have PAML installed in your system and it is accessible

####Keep all the CDS multifasta aligned files (with extension aln) and their associated unrooted species trees in the same folder

###List all the alignment files first
ls|grep aln > allalnfiles


cm=F3x4
models=M7vsM8
cfreq=2
ns=`echo 7 8`

#####Making the control file (a config file for PAML) for each alignment file
for file in `cat allalnfiles`
do
echo $file
j=`echo $file|sed 's/_similar_length_orf.aln//g'`
tree=`echo "$j".nwk`
cp demo.ctl "$j"_F3x4_M7vsM8.ctl
sed -i "s/ssssss/$file/g" "$j"_F3x4_M7vsM8.ctl
sed -i "s/tttttt/$tree/g" "$j"_F3x4_M7vsM8.ctl
sed -i "s/CF/$cfreq/g" "$j"_F3x4_M7vsM8.ctl
sed -i "s/oooooo/$file.F3x4.M7vsM8.out/g" "$j"_F3x4_M7vsM8.ctl
sed -i "s/nnnnnn/$ns/g" "$j"_F3x4_M7vsM8.ctl
done

########The above script will make the control files for all the alignment files present in the directory

####to run codeml####
ls|grep "M7vsM8.ctl" > allctlfiles

for file in `cat allctlfiles`
do
codeml $file
done

###############Now compiling all the outputs of PAML for positively-selected sites##############
ls|grep "aln.F3x4.M7vsM8.out" > ../allcodemloutfiles
cd ..
while read j
do
gene=`echo $j|cut -f1 -d"_"`
clade=`echo $j|cut -f2 -d"_"`
outfile=`echo all_codeml_output/$j`
echo $outfile
j1=`grep -n "^The grid" $outfile |cut -f1 -d ':'` 
k1=`expr $j1 - 3`
m1=`grep -n "BEB" $outfile |cut -f1 -d ':'`
n1=`expr $m1 + 5`
o1=`expr $k1 - $n1`
if [ $o1 -gt 1 ]
then
head -n $k1 $outfile |tail -n $o1 |awk '{print $1"\t"$2"\t"$3}' |sed 's/*//g' |awk '$3>0.95{print c,g,$1,$2,$3}' g=$gene c=$clade OFS='\t' >> sitemodel.outresults_together.txt
fi
done < allcodemloutfiles 2> sitemodellogfile.Log

###########################################################


