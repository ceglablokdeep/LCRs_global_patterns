####This script will go through the aligned CDS files and write their gene and protein accession numbers in a file
##working directory: cd /path/to/alnfiles_trees
for clades in afrotheria amphibia artiodactyla aves carnivore chiroptera marsupials perissodactyla primates rodents squamata testudines
do
ls|grep "$clades"_"similar_length_orf.aln_oneliner" > ../"$clades"_aligned_files
done

##example of a file name: PHF21B_afrotheria_similar_length_orf.aln_oneliner
cd /path/to/alnfiles_trees
for clade in afrotheria amphibia artiodactyla aves carnivore chiroptera marsupials perissodactyla primates rodents squamata testudines
do
while read j
do
echo "$j"
gene=`echo $j|cut -f1 -d_`
grep ">" "$j"|sed 's/>//g' > "$clade"_sp_list
while read k
do
pro_ac=`grep "$k" /path/to/all_datasets/"$gene"_dataset/data/"$gene"_"$clade"_list_orf.fa|cut -f1,2 -d_|sed 's/>//g'`
gene_ac=`grep "$pro_ac" /path/to/all_datasets/"$gene"_dataset/data/data_table.tsv|cut -f11`
echo -e "$clade\t$gene\t$k\t$pro_ac\t$gene_ac" >> ../"$clade"_proac_gene_ac
done < "$clade"_sp_list
done < ../"$clade"_aligned_files
done
cd ..
awk 'NF==5{print $0}' *_proac_gene_ac > all_genes_proac_geneac_combined
sort -b -k2,2 -k1,1 all_genes_proac_geneac_combined > sorted_all_gene_pro_ac
##2904605
###The spreadsheet is unable to load the full file as it exceeds the number of rows, so we provide it as a separate text file.

