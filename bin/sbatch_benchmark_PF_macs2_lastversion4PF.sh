#!/bin/bash - 
#===============================================================================
#
#          FILE: sbatchLike_deepTools.sh
# 
#         USAGE: ./sbatchLike_deepTools.sh 
# 
#   DESCRIPTION: 
# 
#       OPTIONS: ---
#  REQUIREMENTS: ---
#          BUGS: ---
#         NOTES: ---
#        AUTHOR: Dr. Biter Bilen (bb), biterbilen@yahoo.com
#  ORGANIZATION: Stanford, Palo Alto, USA
#       CREATED: 05/05/2016 16:47
#      REVISION:  ---
#===============================================================================

set -o nounset                              # Treat unset variables as an error
source ~/Projects/biter_biter_shared/bin/auxslurm.sh

auxCalcNucSpacing() {
	fle=$1
	name=$2
	ofle=$3
	nc=${4:-1}

	samtools view -@ $nc $fle $selchr | \
		bioawk -c sam -v name=$name '$tlen>0 { h[$tlen]++; } END { OFS="\t"; for (e in h) { print e,h[e],name } }' > $ofle
}

# http://stackoverflow.com/questions/6660010/bash-how-to-assign-an-associative-array-to-another-variable-name-e-g-rename-t
PF() {
	echo FUNCTION $FUNCNAME $* $(date):

	local analysisdirbase=$1
	local indir=${o[indir]:-RM}
	local selchrom=${o[selchrom]:-chr19}
	local names=($(bioawk -tc hdr 'NR>1 && $Pair != 2 { print $Name; }' metadata))
	local cnames=($(bioawk -tc hdr 'NR>1 && $Pair != 2 { print $Cname; }' metadata))
#	local name=$(bioawk -tc hdr -v i=$((${o[stiv]}+1)) 'NR>1 && $Pair != 2 { k+=1; if (k==i) print $Name; }' metadata)
#	local pair=$(bioawk -tc hdr -v i=$((${o[stiv]}+1)) 'NR>1 && $Pair != 2 { k+=1; if (k==i) print $Pair; }' metadata)
	ln -sf ${o[datadir]}/$indir $indir

	if [ $((${o[stiv]}+1)) -eq 1 ]; then
		# calculate nuc spacing
		fles=(${names[@]/%/.sorted.bam})
		fles=(${fles[@]/#/$indir\/})
		local ofleTag=nucSpacing
		local r=$RANDOM
		export -f auxCalcNucSpacing
		parallel --delay 1 --plus --tmpdir . -j0 --header : --tag --xapply \
			"auxCalcNucSpacing {fle} {name} $r.{name} ${o[nc]}" \
			::: fle ${fles[@]} ::: name ${names[@]} 
		cat $r.* > $ofleTag.txt 
		less $ofleTag.txt | bedtools groupby -g 3 -c 2 -o sum > $ofleTag.sum.txt
		Rscript -e "source('$PI_SCRATCH/Projects/biter_biter_shared/bin/auxggplot2.R')" \
						-e "ggplot <- set.ggplot()" \
						-e "a <- read.table('"$ofleTag.txt"')" \
						-e "b <- read.table('"$ofleTag.sum.txt"')" \
						-e "dta <- merge(a,b,by.x='V3',by.y='V1')" \
						-e 'dta$group = sub("_(\\d+)$", "", dta$V3)' \
						-e 'ggplot(subset(dta,V1<1000), aes(x=V1, y=V2.x/V2.y, color=V3, shape=group)) + geom_line() + labs(y="Normalized Read Count", x="Fragment length(bp)")' \
					  -e "ggsave('"$ofleTag.pdf"')"
		#rm $r.*

		# calculate similarity
		# non-overlapping windows of 1K with >20 reads
		local ofleTag=similarityW${mwW}S${mwS}MRC${minReadCutoff}
		local genomef=$HOME/aux/genomes/${genomel:-mouse.mm10}.genome
		grep -w $selchrom $genomef > genome.$selchrom
		bedtools makewindows -g genome.$selchrom -w ${o[makewindowsW]:-1000} -s ${o[makewindowsS]:-1000} > genome.$selchrom.bed
		fles=(${names[@]/%/.sorted.bam})
		echo Chr Start End ${fles[@]} | sed 's/ /\t/g' > genome.$selchrom.bedcov
		fles=(${fles[@]/#/$indir\/})
		samtools bedcov genome.$selchrom.bed ${fles[@]} >> genome.$selchrom.bedcov
		gzip -f genome.$selchrom.bedcov
		Rscript -e "source('$PI_SCRATCH/Projects/biter_biter_shared/bin/auxML.R')" \
			-e "param.preprocess <- list(col.sel.pat='.sorted.bam$', colname.split.rank=2, colname.split.pat='.')" \
			-e "spike <- list(spikeId='ERCC', rm=T)" \
			-e "param.process <- list(method='filterRowMeanMin', cutoff=${o[minReadCutoff]:-20})" \
			-e "dta.exp <- caller.getExpressed('genome.$selchrom.bedcov.gz', param.preprocess=param.preprocess, ofle.tag='$ofleTag', nc=${o[nc]}, spike=spike, density.norm.colname=NULL, param.process=param.process)"

	fi
#	source ~/PI_HOME/PythonSandbox/MACS2/bin/activate
#	param="-g ${o[gsize]:-mm} -f AUTO"
#	# pooled index 0
#	if [ ${o[poolReps]:-T} == T ]; then
#		ofle=pooled
#	fi
 #	param+=" -n $ofle"

#	macs2 callpeak -t XX -c XX
#	macs2 --version
#	deactivate
	
	bioawk --version
	samtools --version
}

# 4 Ling - Hairless paper submission
slurmMain PF --datadir=$PI_SCRATCH/Projects/biter_lingliu_DNAdamage/data/Liu_2017/ATACseq --t1=HrCKO_R_SC_1,HrCKO_R_SC_2,HrCKO_R_SC_3,HrCKO_S_SC_1 --t2=HrWT_R_SC_1,HrWT_R_SC_2,HrWT_R_SC_3,HrWT_S_SC_1 
slurmMain PF --datadir=$PI_SCRATCH/Projects/biter_lingliu_DNAdamage/data/Liu_2016/ATACseq
#slurmMain PF --datadir=$PI_SCRATCH/Projects/biter_lingliu_DNAdamage/data/Liu_2016/ATACseq --indir=RM_TR
#slurmMain PF --jobnames=calcsim --datadir=~/Projects/biter_lingliu_DNAdamage/data/Liu_2016/ATACseq --indir=RM_TR --t1=HrCKO_SC_1,HrCKO_SC_2 --t2=HrWT_SC_1,HrWT_SC_2 
#slurmMain PF --jobnames=calcsim,predictd,filterdup,pileup,locallambda,callpeak --datadir=~/Projects/biter_demorree_PAS/data/Soleimani_2012_DevCell/ChIPseq --t1=C2C12_Pax3TAP_3xFLAG_1,C2C12_Pax7TAP_3xFLAG_1

# 4 Antoine
#slurmMain PF --jobnames=calcsim,predictd,filterdup,pileup,locallambda,callpeak --datadir=~/Projects/biter_demorree_PAS/data/Soleimani_2012_DevCell/ChIPseq --nc=4 --otag=poolreps --t1=C2C12_Pax3TAP_3xFLAG,C2C12_Pax7TAP_3xFLAG
#slurmMain PF --jobnames=calcsim,predictd,filterdup,pileup,locallambda,callpeak --datadir=~/Projects/biter_demorree_PAS/data/Soleimani_2012_DevCell/ChIPseq --t1=C2C12_Pax3TAP_3xFLAG_1,C2C12_Pax7TAP_3xFLAG_1
	   
# 4 Jay
#slurmMain PF --jobnames=calcsim,predictd,filterdup,pileup,locallambda,callpeak --datadir=~/Projects/biter_jsalvi_spuriousReplication/data/Dellino_2013_GenomeResearch/ChIPseq --gsize=hs --genomel=human.hg38 \
#	--t1=HeLa_ORC1_LDF_1,HeLa_NA_HDF_1
#slurmMain PF --jobnames=calcsim,predictd,filterdup,pileup,locallambda,bdgdiff --datadir=~/Projects/biter_jsalvi_spuriousReplication/data/Salvi_2016/ChIPseq \
#	--t1=NA_SC_ATR_1,05hdPI_SC_ATR_1,NA_SC_ORC1_1,05hdPI_SC_ORC1_1,NA_SC_ATR_1,05hdPI_SC_ATR_1,NA_SC_ATR_2,05hdPI_SC_ATR_2,NA_SC_ATR_1,NA_SC_ORC1_1,NA_SC_ATR_2,NA_SC_ORC1_2 \
#	--t2=NA_SC_ATR_2,05hdPI_SC_ATR_2,NA_SC_ORC1_2,05hdPI_SC_ORC1_2,NA_SC_ORC1_1,05hdPI_SC_ORC1_1,NA_SC_ORC1_2,05hdPI_SC_ORC1_2,05hdPI_SC_ATR_1,05hdPI_SC_ORC1_1,05hdPI_SC_ATR_2,05hdPI_SC_ORC1_2
#slurmMain PF --jobnames=callpeak --datadir=~/Projects/biter_jsalvi_spuriousReplication/data/Salvi_2016/ChIPseq \
#	--t1=NA_SC_ATR_1,05hdPI_SC_ATR_1,NA_SC_ORC1_1,05hdPI_SC_ORC1_1,NA_SC_ATR_2,05hdPI_SC_ATR_2,NA_SC_ORC1_2,05hdPI_SC_ORC1_2

# 4 Ling
#slurmMain PF --jobnames=calcsim,predictd,filterdup,pileup,locallambda,bdgdiff --datadir=~/Projects/biter_lingliu_DNAdamage/data/Liu_2016/ChIPseqV2 \
#	  --t1=HrCKO_young_TMX_SC_NA_H3K9me2_1,HrCKO_young_TMX_SC_NA_H3K9me3_1,HrWT_young_TMX_SC_NA_H3K9me2_1,HrWT_young_TMX_SC_NA_H3K9me3_1 \
#	  --t2=HrCKO_young_TMX_SC_NA_H3K9me2_2,HrCKO_young_TMX_SC_NA_H3K9me3_2,HrWT_young_TMX_SC_NA_H3K9me2_2,HrWT_young_TMX_SC_NA_H3K9me3_2

#slurmMain PF --jobnames=callpeak --datadir=~/Projects/biter_lingliu_DNAdamage/data/Liu_2016/ChIPseqV2 \
#	  --t1=HrCKO_young_TMX_SC_NA_H3K9me2_1,HrCKO_young_TMX_SC_NA_H3K9me3_1,HrWT_young_TMX_SC_NA_H3K9me2_1,HrWT_young_TMX_SC_NA_H3K9me3_1,HrCKO_young_TMX_SC_NA_H3K9me2_2,HrCKO_young_TMX_SC_NA_H3K9me3_2,HrWT_young_TMX_SC_NA_H3K9me2_2,HrWT_young_TMX_SC_NA_H3K9me3_2

#slurmMain PF --jobnames=filterdup,pileup,locallambda,bdgdiff --datadir=~/Projects/biter_lingliu_DNAdamage/data/Liu_2016/ChIPseqV2 --otag=poolreps \
#	  --t1=HrCKO_young_TMX_SC_NA_H3K9me2,HrCKO_young_TMX_SC_NA_H3K9me3,HrCKO_young_TMX_SC_NA_H3K9me2,HrWT_young_TMX_SC_NA_H3K9me2 \
#	  --t2=HrWT_young_TMX_SC_NA_H3K9me2,HrWT_young_TMX_SC_NA_H3K9me3,HrCKO_young_TMX_SC_NA_H3K9me3,HrWT_young_TMX_SC_NA_H3K9me3

#slurmMain PF --jobnames=callpeak --datadir=~/Projects/biter_lingliu_DNAdamage/data/Liu_2016/ChIPseqV2 --otag=poolreps \
#	  --t1=HrCKO_young_TMX_SC_NA_H3K9me2,HrCKO_young_TMX_SC_NA_H3K9me3,HrCKO_young_TMX_SC_NA_H3K9me2,HrWT_young_TMX_SC_NA_H3K9me2,HrWT_young_TMX_SC_NA_H3K9me2,HrWT_young_TMX_SC_NA_H3K9me3,HrCKO_young_TMX_SC_NA_H3K9me3,HrWT_young_TMX_SC_NA_H3K9me3 


