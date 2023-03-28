####CompositionMaker is an accessory program alongwith fLPS2.0

##We run the CompositionMaker and fLPS2.0 for each amino acid fasta file in their respective datasets folder
ls all_datasets > all_genes
for dir in `cat all_genes`
do
gene=`echo $dir|sed 's/_dataset//g'`
cd "$dir"/data
for orfs in `ls -1 *_list_orf.fa`
do
prefix=`echo $orfs|cut -f1,2 -d"_"`
CompositionMaker "$orfs"
composfile=`echo "$orfs".COMPOSITION`
fLPS2 -c "$composfile" "$orfs" > "$prefix"_flpsout_withcomposition
done
cd ../..
done
