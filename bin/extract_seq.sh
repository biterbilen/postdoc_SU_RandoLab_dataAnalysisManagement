# AP-1 components and Fosl1 targets 20170419
# Rbpj Heyl Junb have Fosl1 TFBS (ChIP-seq in C2C12) Rbpj and Heyl chnage expression in aging QSC (Ikeda Brett RNAseq)
# Ccnd1 chmaged protein expression in Fosl1 cKO QSCs (Brett Western (x1))
genes="Fosl1|Fosl2|Fosb|Fos|Junb|Ccnd1|Rbpj|Heyl"; ofle=AP1_Fosl1targets.CDS.fa
less mm10_refGene.seltab.gz | grep -wP "$genes" | awk -F "\t" '$4 != "n/a" { OFS="\t"; print $6,$2,$1 }' | \
	while read i; do
		id=$(echo $i | cut -d " " -f 1) 
		symbol=$(echo $i | cut -d " " -f 2)
		refid=$(echo $i | cut -d " " -f 3) 
		bioawk -c fastx -v id=$id -v symbol=$symbol -v refid=$refid '$name"," ~ id { gsub(/([[:lower:]]+)/, "", $seq); print ">"$name" "symbol" "refid"\n"$seq } ' mm10_ensGene.mRNA.fasta.gz
	done > $ofle 
# did not work  could be gawk awk incomplete implementation
#		bioawk -c fastx -v id=$id -v symbol=$symbol ' $name ~ id && match($seq,/([[:lower:]]+)([[:upper:]]+)([[:lower:]]+)/,k) > 0 { print $name, symbol, k[2], k[1], k[3]; } ' mm10_ensGene.mRNA.fasta.gz


less /home/biter/Projects/biter_jbrett_DNAmethylation/results/2017-02-16/DEOvYquiescenceWM_k1_scupper_Brett_2015_RNAseq/norm_DEG_scupper_RSRUVs_contrast000-100010_DE_edgeR.txt.gz | perl ~/Projects/biter_biter_EWSR1/bin/scripts/_PAPD5/innerJoin_stdin.pl /home/biter/Projects/biter_jbrett_epigenetics/results/2014-11-28/back.bs10/STAR_Q4_Q18_DE.txt 1 1 '' | cut -f 1,3,6-7,24- | awk '$NF<0.001 && $3<0.001{ print $4}' | sort | sort > filter2

b=~/Projects/biter_biter_stemCellFateRegulation/results/2016-06-19_TFBS_prediction_SCSB/Wold_2011/ChIPseq/FOSL1-mouse.bed.gz
over=/home/biter/PI_HOME/Data/casco/UCSC_tracks/mm10/mm9ToMm10.over.chain.gz
ob=mm10_${b##*/}
ob=${ob%.*}
less $b | cut -f 1-6 > tmp_$ob
liftOver tmp_$ob $over $ob $ob.unmapped;
less $ob | bedtools sort -i stdin | gzip -c > $ob.gz
less ~/PI_HOME/Data/casco/UCSC_tracks/mm10/mm10_refGene_canonical.v2.gtf.gz | grep start_codon | bedtools sort -i stdin | bedtools closest -D a -id -a stdin -b $ob.gz | awk ' $13 != "." && $NF > -10000 && $NF < 1 { print $10}' |  sed 's/";*//g' | sort | uniq > filter1

cat filter* | sort | uniq -D | uniq > Fosl1_targets.id

rm mm10_${ob} tmp_${ob}

