#!/bin/bash -l 
#===============================================================================
#
#          FILE: sbatch_processRNAseqwbaySeq.sh
# 
#         USAGE: ./sbatch_processRNAseqwbaySeq.sh 
# 
#   DESCRIPTION: 
# 
#       OPTIONS: ---
#  REQUIREMENTS: ---
#          BUGS: ---
#         NOTES: ---
#        AUTHOR: Dr. Biter Bilen (bb), biterbilen@yahoo.com
#  ORGANIZATION: Stanford, Palo Alto, USA
#       CREATED: 12/04/2014 20:44
#      REVISION:  ---
#===============================================================================

set -o nounset                              # Treat unset variables as an error

projectDir=`pwd | awk -F "results" '{ print $1}'`
gtf=$projectDir/data/UCSC_tracks/mm10/mm10_refGene_canonical.gtf.gz

tag=TopHat_Q4_Q18; p1=Q18; p2=Q4; FWER=0.00001; metafile=../../data/Ikeda_2014/metadataMika1; s=0
tag=STAR_Q4_Q18; p1=Q18; p2=Q4; FWER=0.00001; metafile=../../data/Ikeda_2014/metadataMika1; s=0
tag=STAR_Q18_X18; p1=Q18; p2=X18; FWER=0.01; metafile=../../data/Ikeda_2014/metadataMika1; s=0
tag=STAR_Q18_X18; p1=Q18; p2=X18; FWER=0.01; metafile=../../data/Ikeda_2014/metadataMika1; s=0
tag=STAR_Q4_X4; p1=Q4; p2=X4; FWER=0.01; metafile=../../data/Ikeda_2014/metadataMika1; s=0
tag=STAR_Liu_Q24_Q3; p1=Q24; p2=Q3; FWER=0.0001; metafile=../../data/Liu_2014/RNAseq/metadata; s=0
tag=STAR_Liu_Q12_Q3; p1=Q12; p2=Q3; FWER=0.01; metafile=../../data/Liu_2014/RNAseq/metadata; s=0
tag=STAR_Liu_Q2_Q3; p1=Q2; p2=Q3; FWER=0.01; metafile=../../data/Liu_2014/RNAseq/metadata; s=0
tag=STAR_Liu_Q12_Q24; p1=Q12; p2=Q24; FWER=0.01; metafile=../../data/Liu_2014/RNAseq/metadata; s=0
tag=STAR_Cheung_A3_Q3; p1=A3; p2=Q3; FWER=0.00001; metafile=../../data/Cheung_2014/RNAseq/metadata; s=0
tag=STAR_Leeman_PG34_PG1922; p1=PG34; p2=PG1922; FWER=0.001; metafile=../../data/Leeman_2015/RNAseq/metadata; s=0
tag=STAR_Leeman_PGE34_PGE1922; p1=PGE34; p2=PGE1922; FWER=0.001; metafile=../../data/Leeman_2015/RNAseq/metadata; s=0
tag=STAR_Ryall_SIRT1mKOCul_SIRT1mKOFI; p1=SIRT1mKOCul; p2=SIRT1mKOFI; FWER=0.001; metafile=../../data/Ryall_2015_Cell/RNAseq/metadata; s=0
tag=STAR_Ryall_WTCul_WTFI; p1=WTCul; p2=WTFI; FWER=0.001; metafile=../../data/Ryall_2015_Cell/RNAseq/metadata; s=0
tag=STAR_Sun_mo4_mo24; p1=m04; p2=m24; FWER=0.001; metafile=../../data/Sun_2014_CellStemCell/RNAseq/metadata; s=0

#TODO correct above for indata
tag=HISAT_Wosczyna_KO_WT; p1=KO; p2=WT; FWER=0.05; indata=$projectDir/data/Wosczyna_2015/RNAseq; metafile=$indata/metadata; s=1
tag=HISAT_Wosczyna_KO_WT; p1=KO; p2=WT; FWER=0.05; indata=$projectDir/data/Wosczyna_2015/RNAseq; metafile=$indata/metadata; s=0; outdir=$tag\s$s
tag=HISAT_Wosczyna_KO_WT; p1=KO; p2=WT; FWER=0.05; indata=$projectDir/data/Wosczyna_2015/RNAseq; metafile=$indata/metadata; s=2; outdir=$tag\s$s
tag=HISAT_Liu_KO_WT; p1=KO; p2=WT; FWER=0.05; indata=$projectDir/data/Liu_2015/RNAseq; metafile=$indata/metadata; s=0; outdir=$tag

#tag=STAR_Ikeda_Liu_ALL; FWER=0.01;

nc=16; bs=40; t=14
nc=8; bs=20; t=14
nc=8; bs=10; t=14
nc=16; bs=20; t=16
nc=8; bs=20; t=16
nc=8; bs=40; t=16
nc=16; bs=10; t=48
nc=8; bs=20; t=16
nc=8; bs=20; t=48
norm=edgeR

#sbatch -J featureCountsBayseq -p normal -t $t:0:0 --cpus-per-task=$nc <<SBATCH
##!/bin/bash -l
#
#	paste STAR_Q4_Q18.out STAR_Liu_Q24_Q3.out STAR_Liu_Q12_Q3.out STAR_Liu_Q2_Q3.out | \
#		cut -f 1-19,26-29,36-38,47-48 | \
#		sed 's/STAR_Q4_Q18/STAR_Ikeda_Liu_ALL/g' | sed 's/STAR_Liu_Q24_X3/STAR_Ikeda_Liu_ALL/g' | sed 's/STAR_Liu_Q12_Q3/STAR_Ikeda_Liu_ALL/g' | sed 's/STAR_Liu_Q2_Q3/STAR_Ikeda_Liu_ALL/g' > $tag.out
#
#	echo "R --no-save --args $tag.out $nc $bs $FWER $norm < ./baySeq_generic.R"
#	R --no-save --args $tag.out $nc $bs $FWER $norm < ./baySeq_generic.R  
#
#	echo DONE
#SBATCH
#exit

mkdir -p $outdir; pushd $outdir

zless $gtf > $tag.gtf
ln -sf $metafile $tag.meta

less $tag.meta | awk '$4==0 && ($3==p1||$3==p2){ print $1}' p1=$p1 p2=$p2 | \
	while read sample; do 
		d=`awk 'c==$1{ print $2}'	c=$sample $tag.meta`; 
		if [[ $tag =~ STAR_Ryall ]]; then
			ln -sf ../../data/Ryall_2015_Cell/RNAseq/mapSTAR/${s}Aligned.out.bam $tag.$d.bam
		elif [[ $tag =~ STAR_Leeman ]]; then
			ln -sf ../../data/Leeman_2015/RNAseq/mapSTAR/${s}Aligned.out.bam $tag.$d.bam
		elif [[ $tag =~ STAR_Sun ]]; then
			ln -sf ../../data/Sun_2014_CellStemCell/RNAseq/mapSTAR/${s}Aligned.out.bam $tag.$d.bam
		elif [[ $tag =~ STAR_Liu ]]; then
			ln -sf ../../data/Liu_2014/RNAseq/starMapOut/${s}Aligned.out.bam $tag.$d.bam
		elif [[ $tag =~ STAR_Cheung ]]; then
			ln -sf ../../data/Cheung_2014/RNAseq/mapSTAR/${s}Aligned.out.bam $tag.$d.bam
		elif [[ $tag =~ STAR ]]; then
			ln -sf ../../data/Ikeda_2014/starMapOut/${s}Aligned.out.bam $tag.$d.bam
		elif [[ $tag =~ TopHat ]]; then
			ln -sf ../../data/Ikeda_2014/mapOut/$sample/accepted_hits.bam $tag.$d.bam
		elif [[ $tag =~ HISAT ]]; then
			ln -sf $indata/HISAT_mapOut/$d.sorted.bam $tag.$d.bam
		else	
			echo set files for $tag
		fi
	done

ls $tag.*.bam
bamFiles="`ls $tag*bam`"

sbatch -J featureCountsBayseq -p normal -t $t:0:0 --cpus-per-task=$nc <<SBATCH
#!/bin/bash -l

	cmd="/scratch/users/jbrett/RNASeq2_2/featureCounts/subread-1.4.6-Linux-x86_64/bin/featureCounts \
		-s $s -T $nc -p -B -C -t exon -g transcript_id -a $tag.gtf -o $tag.out $bamFiles" 
	echo \$cmd

	if [ ! -e $tag.out ]; then
		\$cmd
	fi

	echo "R --no-save --args $tag.out $nc $bs $FWER $norm < ~/Projects/biter_wosczyna_miRNAsInMuscleCellFate/bin/baySeq_generic.R"
	R --no-save --args $tag.out $nc $bs $FWER $norm < ~/Projects/biter_wosczyna_miRNAsInMuscleCellFate/bin/baySeq_generic.R

	rm $tag.gtf

	echo DONE
SBATCH

popd


