sed -i 's/Marsupials/Marsupialia/g' gcpercompiled_allclades
awk '{print $1,$3,$2}' OFS="\t" gcpercompiled_allclades > gcpercompiled_arranged

chmod a+x gc_comparison.r
Rscript gc_comparison.r
