fdr=0.01
#tag=DE_upper_Wosczyna_2015_PullDown
#tag1=DE_upper_Wosczyna_2015_PullDown
tag=DE_upper_Wosczyna_2015_PullDown
tag1=DE_upper_Wosczyna_2015_RNAseq
files=($(ls $tag/*edgeR.txt.gz))
files1=($(ls $tag1/*edgeR.txt.gz))
ofile=FF_final/$tag\_$tag1.all
echo -e "Id ${files[@]}" | sed 's/ /\t/g' > $ofile
for f in ${files[@]}; do
	echo -en "$f"
	less $f | grep -v "^#" | awk -v fdr=$fdr '$2>0 && $6<fdr{ print $1}' > a
	for f1 in ${files1[@]}; do
		less $f1 | grep -v "^#" | awk -v fdr=$fdr '$2>0 && $6<fdr{ print $1}' > b
		cat a b | sort | uniq -D | uniq > c
		echo `less c | wc -l` `less a | wc -l` `less b | wc -l` | awk '{ printf("\t%.0f|%.0f|%.0f|%.2f", $1,$2,$3,$1/$2) }'
	done #| sort -k 1,1gr
	echo
done >> $ofile
rm a b c
less $ofile | sed 's/OUT./DEG:/g' | sed "s/\/$tag1/|n:/g" | sed "s/\/$tag/|n:/g" | sed 's/_DE_edgeR.txt.gz//g' > tmp
echo "# files:${files[@]}" >> tmp
echo "# files1:${files1[@]}" >> tmp
mv tmp $ofile

exit

dir1=RT
dir2=OUT.TMM
less $dir1/Wosczyna_2015_PullDown_wospikeIn_UQ_RUVs_DE_edgeR.txt.gz  | awk '$2>0 && $6<0.01{ print $1}'  > a
less $dir2/Wosczyna_2015_PullDown_wospikeIn_UQ_RUVs_DE_edgeR.txt.gz  | awk '$2>0 && $6<0.01{ print $1}'  > b
cat a b | sort | uniq -D | uniq > c
wc -l a b c

exit

dir=RT
less $dir/Wosczyna_2015_PullDown_wospikeIn_UQ_RUVs_DE_edgeR.txt.gz  | awk '$2>0 && $6<0.01{ print $1}'  > a
less $dir/Wosczyna_2015_RNAseq_wospikeIn_UQ_RUVs_DE_edgeR.txt.gz  | awk '$2>0 && $6<0.01{ print $1}'  > b
cat a b | sort | uniq -D | uniq > c
wc -l a b c
