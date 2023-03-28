for pur in `seq 0.1 0.1 1`
do
echo $pur
lpsfile=`ls flps_output_morethan_"$pur"_aligned_with_gene_name_omega_values_unique_values`
nolpsfile=`ls flps_output_lessthan_"$pur"_aligned_with_gene_name_omega_values_unique_values`
Rscript omegascript.r "$lpsfile" "$nolpsfile" "$pur"
done
