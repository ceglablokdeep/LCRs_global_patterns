while read j
do
len=`echo $j|awk '{print $2}'`
echo $j|awk '{$1=$2=$3="";print $0}' OFS="\n"|sed '/^$/d'|awk -v l=$len '{print $1,(NR/l)*100}' >> alndensity_withnormpositions.txt
done < alndensity_normalized.txt
#########
for i in `seq 1 100`
do
(j=`echo $i|awk '{print $1-1}'`
echo $i $j
awk -v llim=$j -v ulim=$i '$2>=llim && $2<ulim{print $0}' alndensity_withnormpositions.txt > temp_"$i".txt
awk -v inter=$i 'END{y=NR} {x+=$1} END{print inter-0.5,x/y}' OFS="\t" temp_"$i".txt >> withinterval_"$i".txt)&
done

###Making alignment coverage boxplot wise
rm -f aln_coverage_binwise.txt
for i in `seq 1 100`
do
echo $i
awk -v ran=$i '$2<ran && $2>=(ran-1){print (ran+ran-1)/2,$1}' alndensity_withnormpositions.txt >> aln_coverage_binwise.txt
done

rm -rf final_alndensity.txt
for i in `seq 0.5 99.5`
do
awk -v n=$i '$1==n{print $2}' aln_coverage_binwise.txt > temp.txt
awk -v n=$i '{ sum += $1 } END { if (NR > 0) print n,sum / NR }' temp.txt OFS="\t" >> final_alndensity.txt
done

##########

Rscript plotting.r
