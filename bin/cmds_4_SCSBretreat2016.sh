#!/bin/bash - 
#===============================================================================
#
#          FILE: cmds_4_SCSBretreat2016.sh
# 
#         USAGE: ./cmds_4_SCSBretreat2016.sh 
# 
#   DESCRIPTION: 
# 
#       OPTIONS: ---
#  REQUIREMENTS: ---
#          BUGS: ---
#         NOTES: ---
#        AUTHOR: Dr. Biter Bilen (bb), biterbilen@yahoo.com
#  ORGANIZATION: Stanford, Palo Alto, USA
#       CREATED: 06/23/2016 14:20
#      REVISION:  ---
#===============================================================================

set -o nounset                              # Treat unset variables as an error
source ~/Projects/biter_biter_shared/bin/auxslurm.sh

declare -A o=( [jobnames]=liftOver2mm10 [queue]=normal [t]=1 [nc]=1 [selchrom]=chr19 [outdir]=datachr19 [over]=~/PI_HOME/Data/casco/UCSC_tracks/mm10/mm9ToMm10.over.chain.gz)
declare -A o=( [jobnames]=fimo [queue]=normal [t]=1 [nc]=1 [selchrom]=chr19 [outdir]=fimochr19 [msc]=1000000 [thres]=1-e5 [org]=mm10)
declare -A o=( [jobnames]=chipOverlap [queue]=normal [t]=1 [nc]=1 [selchrom]=chr19 [outdir]=overlap)
declare -A o=( [jobnames]=footprint [queue]=normal [t]=1 [nc]=1 [selchrom]=chr19 [outdir]=footprint [minreadcut]=20 )
declare -A o=( [jobnames]=prepPWMs [queue]=normal [t]=1 [nc]=1 [outdir]=cisdb [orgLatin]=Mus_musculus)

declare -A o=( [jobnames]=prepPWMs [queue]=normal [t]=1 [nc]=1 [outdir]=cisdb [orgLatin]=Homo_sapiens)
declare -A o=( [jobnames]=fimo [queue]=normal [t]=4 [nc]=1 [selchrom]=chr19 [outdir]=fimochr19 [msc]=1000000 [thres]=1-e5 [org]=hg38)


jobname=prepPWMs
if [[ ,${o[jobnames]:-X}, =~ ,$jobname, ]]; then 
	mkdir -p ${o[outdir]}; pushd ${o[outdir]};
	# Prep PWMs -mouse:641 Matrix_IDs ; human:748
	# grep -Piw "CEBPB|CTCF|FOSL1|MAX|MYOD1|MYOG|USF1"
	ln -sf ~/PI_HOME/Data/casco/CIS-DB/${o[orgLatin]}/TF_Information.txt .
	# Select Motif ids w/ logos
	ls ~/PI_HOME/Data/casco/CIS-DB/${o[orgLatin]}/logos_all_motifs/ | grep _fwd.png | sed 's/_fwd.png//g' | sort | uniq > Motif_ID.wlogos 
	# Select for each TF_Name one representative (alphabetically the last) Motif_ID and TF_Name among the PWMs w/ logos
	grep -w -f Motif_ID.wlogos TF_Information.txt | awk '$4 != "." && NR>1 { print }' | cut -f 4,7 | grep ^M | sort -k 2,2 -k 1,1r | uniq -f 1 | \
		sort -k 1,2 | bedtools groupby -g 1 -c 2 -o collapse -delim "|" > Motif_ID.TF_Names 
	allmotifs=~/PI_HOME/Data/casco/CIS-DB/${o[orgLatin]}/all_pwms.txt
	Rscript -e "source('~/Projects/biter_biter_shared/bin/auxTFBSTools.R')" \
		-e "getMatrixList(src='$allmotifs', gnm='hg38', matrixtype='PWM', ii=1, ni=1, byrow=T, select='Motif_ID.TF_Names', ofle='CIS-DB.selected')"
	popd
fi

# Check rate of ChIP-seq peaks in PWM mapping sites
#homer2 gave segmentation fault!!! using fimo instead
jobname=fimo
if [[ ,${o[jobnames]:-X}, =~ ,$jobname, ]]; then 
	mkdir -p ${o[outdir]}; pushd ${o[outdir]};
	less ~/PI_HOME/Data/casco/UCSC_tracks/${o[org]}/assembly/${o[selchrom]}.fa* > ${o[selchrom]}.fa
	grep ^MOTIF ../cisdb/CIS-DB.selected | \
		while read i; do
			m=$(echo $i | awk '{ print $2}')
			o=$(echo $i | awk '{ print $3}' | sed 's/|/__/g')
			echo $o
			sbatch -J $jobname -p ${o[queue]} -t ${o[t]}:0:0 -c ${o[nc]} -e slurm_$o.out -o slurm_$o.out <<- END
			#!/bin/bash -l
			fimo --thresh ${o[thres]} --max-stored-scores ${o[msc]} --motif $m --oc $o ../cisdb/CIS-DB.selected ${o[selchrom]}.fa
			less $o/fimo.txt | awk 'BEGIN{OFS="\t"} NR>1 { print \$2,\$3-1,\$4,\$1,\$6,\$5}' | \
				sort -k 1,1 -k 2,2g | gzip -c > $o.bed.gz
			#rm -rf $o
			echo DONE 
			END
	done
	popd
fi
# TODO run after the jobs finish
# ls */*txt | xargs wc -l > PWM_match_counts

# liftover mm9 data to mm10
jobname=liftOver2mm10
if [[ ,${o[jobnames]:-X}, =~ ,$jobname, ]]; then 
	mkdir -p ${o[outdir]}; pushd ${o[outdir]};
	# TODO set
	cdir=../../../data/Wold_2011/ChIPseq;                 cs=($(ls $cdir/*-mouse.bed.gz));     ln -sf $cdir .;declare -A tags=([CEBPB]=Cebpb_MB [CTCF]=Ctcf_MB [FOSL1]=Fosl1_MB [MAX]=Max_MB [MYOD1]=Myod1_MB [MYOG]=Myog_MB [USF1]=Usf1_MB)
	cdir=../../../data/DellOrso_2016_CellReports; cs=($(ls $cdir/*PBX1*peaks.bed.gz)); ln -sf $cdir .;declare -A tags=([PBX1_MB]=Pbx1_MB [PBX1_MT]=Pbx1_MT [MB_ATACSeq_rep1]=ATAC_MB_1 [MB_ATACSeq_rep2]=ATAC_MB_2 [MT_ATACSeq_rep1]=ATAC_MT_1 [MT_ATACSeq_rep2]=ATAC_MT_2)
	for tag in ${!tags[@]}; do
		otag=${tags[$tag]}
		if [[ $otag =~ ATAC ]]; then
			f=$(ls $cdir/*$tag*peaks.txt.gz)
			zless $f | grep -vP "^#|^$|^chr\t" | awk 'BEGIN{ OFS="\t";} { print $1,$2,$3,$9,$8,"." }'
		else
			f=$(ls $cdir/*$tag*bed.gz)
			zless $f | awk 'BEGIN{ OFS="\t";} { print $1,$2,$3,$4,$7,$6 }'
		fi | grep -w ${o[selchrom]} | sort -k 1,1 -k 2,2g > $otag\_mm9.bed
		liftOver $otag\_mm9.bed $(auxDirNameTune ${o[over]}) $otag.bed $otag.unmapped
		gzip -f $otag.bed
		rm $otag\_mm9.bed
		if [ ! -s $otag.unmapped ]; then rm $otag.unmapped; fi 
	done
	popd
fi

# 40M single locus mapping reads does not seem to be sufficient for the modeling
jobname=footprint
if [[ ,${o[jobnames]:-X}, =~ ,$jobname, ]]; then 
	mkdir -p ${o[outdir]}; pushd ${o[outdir]};
	adir=~/Projects/biter_biter_stemCellFateRegulation/data/DellOrso_2016_CellReports/ATACseq/RM/; as=($(ls $adir/*sorted.bam)); ln -sf $adir .; 
	fdir=../fimochr19; declare -A fs=([CEBPB]=Cebpb [CTCF]=Ctcf [FOSL1]=Fosl1 [MAX]=Max [MYOD1]=Myod1 [MYOG]=Myog [USF1]=Usf1); ln -sf $fdir .; 
	fii=MYOD1
	fii=MAX
	ai=2



	echo ${as[$ai]}  ../fimochr19/${fs[$fii]}.bed.gz

	# TODO Get 5' end??? TODO
	# Gets 1st pair w -f flag
	samtools view -f 0x0040  -b ${as[$ai]} chr19 | bedtools bamtobed | \
		bedtools flank -s -l 4 -r 0 -g ~/aux/genomes/mouse.mm10.genome | bedtools flank -s -l 1 -r 0 -g ~/aux/genomes/mouse.mm10.genome > tmp.bed
	span=50
	l=$(less ../fimochr19/${fs[$fii]}.bed.gz | head -n 1 | awk -v span=$span '{ print $3-$2+span+span}')
	bedtools slop -l $span -r $span -i ../fimochr19/${fs[$fii]}.bed.gz -g ~/aux/genomes/mouse.mm10.genome | \
		awk -v l=$l '$3-$2 == l{ print }' | \
		bedtools coverage -counts -a stdin -b tmp.bed > processme
	less processme | \
		awk -v m=${o[minreadcut]} 'BEGIN{OFS="\t"}$7>m{ print $1,$2,$3,$4,$5,$6 }' | \
		bedtools coverage -d -a stdin -b tmp.bed | \
		cut -f 6- | \
		awk -v l=$l -v ps=1 'BEGIN{OFS="\t";}{if ($1=="-"){ $2=l-$2+1; }; print $2,log($3+ps)/log(2); }' | \
		sort -k 1,1g -k 2,2g > processme2

	ops=min,max,sum,count,mode,median,mean,sstdev
	echo pos,$ops | sed 's/,/\t/g' > footprint.${fs[$fii]}
	less processme2 | bedtools groupby -g 1 -c 2 -o min,max,sum,count,mode,median,mean,sstdev >> footprint.${fs[$fii]}

	# TODO do it with bedmap	
	# http://bedops.readthedocs.io/en/latest/content/reference/statistics/bedmap.html
	#bedmap --chrom ${o[selchrom]}   ${as[0]} $fdir/${fs[0]}.bed.gz 
	popd
fi


# Check rate of ChIP-seq peaks in ATAC-seq peaks
jobname=chipOverlap
if [[ ,${o[jobnames]:-X}, =~ ,$jobname, ]]; then 
	mkdir -p ${o[outdir]}; pushd ${o[outdir]};
	tag=_MT
	tag=_MB
	cs=($(ls ../datachr19/*bed.gz | grep -v ATAC | grep $tag))
	as=($(ls ../datachr19/*bed.gz | grep ATAC | grep $tag))
	for c in ${cs[@]}; do
		p=$(echo $c | sed "s/$tag//g" | sed "s/datachr19/fimochr19/g")
		ccount=$(less $c | wc -l)
		pcount=$(bedtools intersect -u -a $c -b $p | wc -l | awk -v t=$ccount '{ printf("%.2f", $1/t) }' )
#		pcount=$(bedtools intersect -u -a $c -b $p | wc -l )
		for a in ${as[@]}; do
			acount=$(bedtools intersect -u -a $c -b $a | wc -l | awk -v t=$ccount '{ printf("%.2f", $1/t) }' )
			k=$(less $a | wc -l | awk '{ printf("%.0f", $1/2)}')
			ahcount=$(less $a | sort -k 5,5gr | awk -v k=$k 'NR<k{ print }' | bedtools intersect -u -a $c -b stdin | wc -l | awk -v t=$ccount '{ printf("%.2f", $1/t) }' )
#			acount=$(bedtools intersect -u -a $c -b $a | wc -l )
			#bedtools fisher -a $a -b $c -g mm10.genome.sorted 
			echo $(basename $c .bed.gz) $(basename $p .bed.gz) $(basename $a .bed.gz) $ccount $pcount $acount $ahcount
		done
	done | sed 's/ /\t/g' > log.intersect$tag
	#less log.fisher$tag | grep -P "Doing|in -a" | awk '{if (NR%3){ ORS="\t"; } else {ORS="\n";} print }' | \
	#	awk 'BEGIN{OFS="\t"}{ t=$8+$17; print $2,$3, t, $8/t }' | \
	#	sort -k 1,2 | bedtools groupby -g 1 -c 4 -o mean | sort -k 2,2gr| \
	#	awk '{ printf("%s\t%.2f\n", $1,$2) }' > rateOfChIPinATAC$tag.txt
	#popd
fi


