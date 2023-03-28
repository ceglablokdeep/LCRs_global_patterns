#######Making XP and XM of all the sequences used in fLPS2.0
##cd /path/to/flps_output
##all_flps_out_files contains names of all the output files of fLPS: 214904

####Example of how fLPS2.0 output looks like
##XP_033781150.1	Geotrypetes_seraphini	368	SINGLE	2	304	348	12	7.949e-05	{L}	----	3.505
while read j
do
echo $j
gene=`echo $j|cut -f1 -d_`
clade=`echo $j|cut -f2 -d_`
sed 's/_/\t/2' "$j"|awk -v g=$gene -v c=$clade '{print c,g,$1,$2,$3,$6,$7,$8,$10}' OFS="\t" >> ../all_flps_output_with_gene_clade_name
done < all_flps_out_files
cd ..
###20384289 all_flps_output_with_gene_clade_name
###As the file has a large number of rows, we will uniquely sort the output
awk '{print $1,$2,$3,$4}' OFS="\t" all_flps_output_with_gene_clade_name|sort -u > all_flps_output_sp_xp_clade_gene
##2739712 number of rows
##afrotheria	A1BG	XP_004387000.1	Trichechus_manatus_latirostris
rm -rf all_flps_output_with_xp_xm_gene_sp_clade
while read j
do
echo $j
clade=`echo $j|awk '{print $1}'`
gene=`echo $j|awk '{print $2}'`
xp=`echo $j|awk '{print $3}'`
sp=`echo $j|awk '{print $4}'`
xm=`grep "$xp" all_datasets/"$gene"_dataset/data/data_table.tsv|cut -f11|head -n1`
echo $j|awk -v xm=$xm '{print $0,xm}' OFS="\t" >> all_flps_output_with_xp_xm_gene_sp_clade
done < all_flps_output_sp_xp_clade_gene
#####
paste <(awk '{print $1,$2,$3,$4}' OFS="\t" all_flps_output_sp_xp_clade_gene) <(awk '{print $1,$2,$3,$4}' OFS="\t" all_flps_output_with_xp_xm_gene_sp_clade) -d "\n"|sort|uniq -c|awk '$1!=2{print $2,$3,$4,$5}' OFS="\t" > notrun
sed -i '/^$/d' notrun
while read j
do
echo $j
clade=`echo $j|awk '{print $1}'`
gene=`echo $j|awk '{print $2}'`
xp=`echo $j|awk '{print $3}'`
sp=`echo $j|awk '{print $4}'`
xm=`grep "$xp" all_datasets/"$gene"_dataset/data/data_table.tsv|cut -f11|head -n1`
echo $j|awk -v xm=$xm '{print $0,xm}' OFS="\t" >> all_flps_output_with_xp_xm_gene_sp_clade_rerun
done < notrun

cat all_flps_output_with_xp_xm_gene_sp_clade all_flps_output_with_xp_xm_gene_sp_clade_rerun >> all_flps_output_with_xp_xm_gene_sp_clade
##2739712 all_flps_output_with_xp_xm_gene_sp_clade
