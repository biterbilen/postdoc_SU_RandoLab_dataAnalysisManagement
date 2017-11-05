d=final_$(date -I)
mkdir -p $d; 
cp README.md $d/.

for tag in RNAseq PullDown; do # for Mike
#for tag in RNAseq; do
	files=($(ls *$tag*/summary *$tag*/*pdf GS*$tag*/*gz ./fastTrack4LM/*pdf))
	for f in ${files[@]}; do
		ff=$(echo $f | sed "s/\//_/g");  
		cp $f $d/$ff; 
	done
done

tar -czvf $d.tar.gz $d
#for f in *RNAseq*/summary *PullDown*/summary; do ff=$(echo $f | sed "s/\//_/g");  cp $f $d/$ff; done 
#for f in *RNAseq*/*pdf *PullDown*/*pdf; do ff=$(echo $f | sed "s/\//_/g");  cp $f $d/$ff; done 
