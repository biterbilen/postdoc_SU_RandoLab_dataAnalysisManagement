less ~/Projects/biter_lingliu_DNAdamage/data/Liu_2016/ChIPseq/rawData/metadata | head -n 1 > a
less metadata | cut -f 3- | awk -F "/" '{ print $NF}' | \
	awk -F "\t" 'BEGIN{OFS="\t"}{ print $1,$2"_"$3"_"$8,$4$5,0,"unstranded",$6,$3,"NA","NA","NA","gDNA",$10,"male",$12,$8,$11,"NA", "NA", "NA", $9, "NA", $13, "HighSeq2000" }' | \
	sort -k 7,7 -k 3,3g -k 15,15 | \
	sed 's/-//g'  | sed 's/ //g' | sed 's/,//g' | sed 's/.fastq.gz//g' | sed 's/human//g' >> a

