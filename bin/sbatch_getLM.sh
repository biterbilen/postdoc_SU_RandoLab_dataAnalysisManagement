#!/bin/bash - 
#===============================================================================
#
#          FILE: sbatch_getLM.sh
# 
#         USAGE: ./sbatch_getLM.sh 
# 
#   DESCRIPTION: 
# 
#       OPTIONS: ---
#  REQUIREMENTS: ---
#          BUGS: ---
#         NOTES: ---
#        AUTHOR: Dr. Biter Bilen (bb), biterbilen@yahoo.com
#  ORGANIZATION: Stanford, Palo Alto, USA
#       CREATED: 02/18/2015 16:15
#      REVISION:  ---
#===============================================================================

set -o nounset                              # Treat unset variables as an error

config() {

	debug=T
	gnm=mm10
	minScore=95%
	nc=8; t=16; queue=normal
	txdb=~/Projects/biter_lingliu_DNAdamage/data/UCSC_tracks/mm10/mm10_refGene_canonical.gtf.gz
	mirSeqFile=~/Projects/biter_jbrett_epigenetics/data/miRBase/mature.fa.gz
	designFile=design
	prj=""
	DEG=NULL
	DEG=DESeq

	# depracated
	tag=Cheung
	regionsStr=promoters,threeUTRsByTranscript
	mirExpFile=~/Projects/biter_jbrett_epigenetics/data/Cheung_2014_Nature/TaqManRodentmiRNA_QSC_ASC.txt.gz
	expFilter=T
	expressedTxFile=~/Projects/biter_jbrett_epigenetics/results/2014-11-28/back.bs20_edgeRnorm_Cheung_A3_Q3/STAR_Cheung_A3_Q3_DE.txt
	featureCountFile=~/Projects/biter_jbrett_epigenetics/results/2014-11-28/STAR_Cheung_A3_Q3.out
	metadata=~/Projects/biter_jbrett_epigenetics/data/Cheung_2014/RNAseq/metadata
	controlGroup=Q3
	padj=0.0000001

	tag=Cheung
	featureCountFile=~/Projects/biter_jbrett_epigenetics/results/2014-11-28/STAR_Cheung_A3_Q3.out
	metadata=~/Projects/biter_jbrett_epigenetics/data/Cheung_2014/RNAseq/metadata
	prj=SCinjury3mo
	controlGroup=QSC3
	padj=0.01
#	regionsStr=promoters,threeUTRsByTranscript
#	mirExpFile=~/Projects/biter_jbrett_epigenetics/data/Cheung_2014_Nature/TaqManRodentmiRNA_QSC_ASC.txt.gz

	tag=Liu
	featureCountFile=~/Projects/biter_lingliu_DNAdamage/results/2015-06-15/HISAT_Liu_KO_WT/HISAT_Liu_KO_WT.out
	metadata=~/Projects/biter_lingliu_DNAdamage/data/Liu_2015/RNAseq/metadata
	prj=Hr_cKO
	controlGroup=WT
	padj=0.01
#	regionsStr=promoters,threeUTRsByTranscript
#	mirExpFile=~/Projects/biter_jbrett_epigenetics/data/Cheung_2014_Nature/TaqManRodentmiRNA_QSC_ASC.txt.gz
#	expFilter=F
#	expressedTxFile=~/Projects/biter_lingliu_DNAdamage/results/2015-06-15/LM_wDESeq_Liu_KO_WT/DEseq_WT_KO_DE.txt

	tag=Ikeda
	bamdir=TODO_map_w_HISAT
	featureCountFile=~/Projects/biter_jbrett_epigenetics/results/2014-11-28/STAR_Q4_Q18.out
	metadata=~/Projects/biter_jbrett_epigenetics/data/Ikeda_2014/RNAseq/metadata
	prj=Aging
	controlGroup=QSC4
	padj=0.01

	tag=Wosczyna
	regionsStr=promoters,threeUTRsByTranscript
	regionsStr=threeUTRsByTranscript
	mirExpFile="mmu-miR-206-3p"
	mirExpFile=~/Projects/biter_wosczyna_miRNAsInMuscleCellFate/data/Wosczyna_2013/selected
	mirExpFile=~/Projects/biter_wosczyna_miRNAsInMuscleCellFate/data/Wosczyna_2013/selected2
	mirExpFile=~/Projects/biter_wosczyna_miRNAsInMuscleCellFate/data/Wosczyna_2013/selected3
	mirExpFile=~/Projects/biter_wosczyna_miRNAsInMuscleCellFate/data/Wosczyna_2013/selected4
	expFilter=F
	expressedTxFile=~/Projects/biter_wosczyna_miRNAsInMuscleCellFate/results/2015-05-07/LM_wDESeq_Wosczyna_KO_WT/head_DEseq_WT_KO_DE.txt
	expressedTxFile=~/Projects/biter_wosczyna_miRNAsInMuscleCellFate/results/2015-05-07/LM_wDESeq_Wosczyna_KO_WT/DEseq_WT_KO_DE.txt
	featureCountFile=~/Projects/biter_wosczyna_miRNAsInMuscleCellFate/results/2015-05-07/HISAT_Wosczyna_KO_WTs0/HISAT_Wosczyna_KO_WT.out
	metadata=~/Projects/biter_wosczyna_miRNAsInMuscleCellFate/data/Wosczyna_2015/RNAseq/metadata
	controlGroup=WT
	padj=0.01



	metadata=$bamdir/metadata
	odir=LM_w$DEG\_$tag\_$prj
}

cmd() {
	ln -sf $featureCountFile $metadata .
	txdbUncomp=`basename $txdb .gz`
	less $txdb > $txdbUncomp
	#ln -sf $mirSeqFile $mirExpFile
	#ln -sf $expressedTxFile expressedTxFile
	less $metadata | awk -v prj="$prj" 'NR==1 || (($8 == prj || prj == "NA") && ($4==0) && ($7=="NA" || $7 == 1)) { OFS="\t"; print $2,$3}' > $designFile
	strandSpecific=`less $metadata | awk -v prj="$prj" 'NR>1 || (($8 == prj || prj == "NA") && ($4==0) && ($7=="NA" || $7 == 1)) { print $6}' | sort | uniq`

	if [ $strandSpecific == unstranded ]; then 
		strandSpecific=0
	elif [ $strandSpecific == RF ]; then
		strandSpecific=2
	else 
		echo "ERROR strandSpecific=$strandSpecific"
	fi
	featureCountFile=`basename $featureCountFile`
#	    -e 'source("~/Projects/biter_biter_shared/bin/auxRsubread.R")' \
#			-e 'fc <- caller.featureCounts(bamdir="$bamdir", param.featureCounts=$featureCounts, ofleTag="", mc.cores=$nc)'
#			-e 'caller.getFeatures(expressedTx.file=expressedTx, features=$featureCounts, mc.cores=$nc, debug=$debug)'

	# R params
	featureCounts="list(annot.inbuilt=\"$gnm\", annot.ext=\"$txdbUncomp\",  isGTFAnnotationFile=T, GTF.featureType=\"exon\", GTF.attrType=\"transcript_id\", isPairedEnd=T, requireBothEndsMapped=T, checkFragLength=T, minFragLength=50, maxFragLength=50000, strandSpecific=$strandSpecific, read2pos=NULL)"
	DEG="list(method=\"$DEG\", padj=$padj, control.group=\"$controlGroup\", mfle=\"$designFile\", gid=\"$txdb\")"
	features="list(regions.str=\"$regionsStr\", gnm=\"$gnm\", exp.filter=T, mirExp.file=\"$mirExpFile\", mirSeq.file=\"$mirSeqFile\", min.score=\"$minScore\")"
	Rcmd=$(cat <<END
Rscript \
	    -e 'source("~/Projects/biter_biter_shared/bin/auxML.R")' \
			-e 'expressedTx <- caller.getExpressed(fle="$featureCountFile", dta=fc, param.DEG=$DEG, mc.cores=$nc, debug=$debug)'
END
)
sbatch -J $odir -p $queue -t $t:0:0 --cpus-per-task=$nc <<SBATCH
#!/bin/bash -l 
echo $Rcmd
$Rcmd
fcho DONE
SBATCH
}

# Main
config
mkdir -p $odir; pushd $odir
cmd
popd



