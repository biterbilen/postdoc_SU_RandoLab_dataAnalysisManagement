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

PF2() {
	files=($(ls $( auxDirNameTune ${o[datadir]}/${o[indir]:-RM}, | sed 's/,/\/*.sorted.bam /g' )))
	labels=($( echo ${files[@]##*/} | awk -F / '{ print $NF}' | sed 's/.sorted.bam//g' ))
	# 1.
	jobname=multiBamSummary; qcmd=$(cat <<- END
	multiBamSummary bins -p ${o[nc]} -bs 1000 -b ${files[@]} -l ${labels[@]} \
		-r ${o[selchr]} --outFileName ${o[selchr]}.npz --outRawCounts ${o[selchr]}.raw
	plotPCA -in ${o[selchr]}.npz -o pca.pdf
	plotCorrelation -in ${o[selchr]}.npz -p heatmap -o heatmap.pdf -c spearman
	plotCorrelation -in ${o[selchr]}.npz -p scatterplot -o scatter.pdf -c spearman
	END
	)
	if [[ ,${o[jobnames]}, =~ ,$jobname, ]]; then qcmds+=( [$jobname]="$qcmd" ); fi

	# 2.
	jobname=computeGCBias 
	if [[ ${o[gsize]:-mm} == mm ]]; then
		gsize=1870000000
	elif [[ ${o[gsize]:-mm} == hs ]]; then
		gsize=2700000000
	else
		echo "ERROR: gsize ${o[gsize]:-mm} is not defined"
	fi
	for i in ${!files[@]}; do
		FN=${files[$i]}
		FT=${labels[$i]}
		qcmd=$(cat <<- END
		computeGCBias -p 1 -b ${FN} -r ${o[selchr]} --effectiveGenomeSize $gsize -g \
			~/PI_HOME/Data/casco/UCSC_tracks/${o[gnm]:-mm10}/assembly/all_2bit/chr.2bit \
			-l ${o[flen]:-300} -freq ${FT}_GCbias.txt --biasPlot ${FT}_GCbias.pdf --regionSize ${o[rsize]:-300}
		correctGCBias -p 1 -b ${FN} --effectiveGenomeSize $gsize -g ~/PI_HOME/Data/casco/UCSC_tracks/${o[gnm]:-mm10}/assembly/all_2bit/chr.2bit \
			--GCbiasFrequenciesFile ${FT}_GCbias.txt -o ${FT}_GCbiasCorrected.sorted.bam
		END
		)
		if [[ ,${o[jobnames]}, =~ ,$jobname, ]]; then qcmds+=( [$jobname$FT]="$qcmd" ); fi
	done

	# 3.
	jobname=plotFingerprint; qcmd=$(cat <<- END
	plotFingerprint -p ${o[nc]} -bs ${o[bs]:-1000} --skipZeros \
		-r ${o[selchr]} --plotFile fingerprints.pdf --outRawCounts fingerprints.txt \
		-b ${files[@]} -l ${labels[@]}
	END
	)
	if [[ ,${o[jobnames]}, =~ ,$jobname, ]]; then qcmds+=( [$jobname]="$qcmd" ); fi
	# 4.
	jobname=plotCoverage; qcmd=$(cat <<- END
	plotCoverage -p ${o[nc]} -b ${files[@]} -l ${labels[@]} -r ${o[selchr]} --plotFile coverage.pdf --outRawCounts coverage.txt
	END
	)
	if [[ ,${o[jobnames]}, =~ ,$jobname, ]]; then qcmds+=( [$jobname]="$qcmd" ); fi

	# 5.
	jobname=plotProfile; qcmd=$(cat <<- END
	less ~/PI_HOME/Data/casco/UCSC_tracks/mm10/canonical/transcripts.bed.gz | grep "^${o[selchr]}" > genes_${o[selchr]}.bed
	less ~/PI_HOME/Data/casco/UCSC_tracks/mm10/canonical/transcripts.bed.gz | grep "^chrX" > genes_chrX.bed
	parallel -k --delay 1 --plus --tmpdir . -j0 --header : --tag --joblog $jobname.jlog --xapply \
		bamCoverage -p 1 -b {FN} -o {FT}.bw  ::: FT ${labels[@]} ::: FN ${files[@]}
	computeMatrix scale-regions -p ${o[nc]} -S *bw -R genes_${o[selchr]}.bed genes_chrX.bed -b 1000 -a 1000 -m 1000 -o matrix.mat.gz
	plotProfile -m matrix.mat.gz --numPlotsPerRow 2 -out profile.pdf
	END
	)
	if [[ ,${o[jobnames]}, =~ ,$jobname, ]]; then qcmds+=( [$jobname]="$qcmd" ); fi

}

slurmMain PF2 --skippostproc=F --datadir=~/Projects/biter_lingliu_DNAdamage/data/Liu_2015/ChIPseq2 --indir=RM_TR_ada \
	--selchr=chr10 --gsize=mm --t=10 \
	--jobnames=plotFingerprint
	#--jobnames=multiBamSummary,plotFingerprint,plotCoverage,plotProfile 

#declare -A o=( [jobnames]=computeGCBias [queue]=normal [t]=8 [nc]=8 [indir]=~/Projects/biter_lingliu_DNAdamage/data/Liu_2016/ChIPseq2/RM/ [selchr]=chr10 [gsize]=mm [outdir]=Liu_2016_ChIPseq2 )
#declare -A o=( [jobnames]=multiBamSummary,computeGCBias,plotFingerprint,plotCoverage,plotProfile [queue]=normal [t]=8 [nc]=8 [indir]=~/Projects/biter_lingliu_DNAdamage/data/Liu_2016/ChIPseq/RM/ [selchr]=chr10 [gsize]=mm [outdir]=Liu_2016_ChIPseq )
#declare -A o=( [jobnames]=multiBamSummary,plotFingerprint,plotCoverage,plotProfile [queue]=normal [t]=6 [nc]=8 [indir]=~/Projects/biter_biter_shared/results/2016-05-04/Liu_2016_ChIPseq [selchr]=chr10 [gsize]=mm [outdir]=Liu_2016_ChIPseq_wGCbiasCorrected )
#declare -A o=( [jobnames]=multiBamSummary,plotFingerprint,plotCoverage,plotProfile [queue]=normal [t]=6 [nc]=8 [indir]=~/Projects/biter_biter_shared/results/2016-05-04/Liu_2016_ChIPseq2 [selchr]=chr10 [gsize]=mm [outdir]=Liu_2016_ChIPseq2_wGCbiasCorrected )
#declare -A o=( [jobnames]=multiBamSummary,computeGCBias,plotFingerprint,plotCoverage,plotProfile [queue]=normal [t]=10 [nc]=8 [indir]=~/Projects/biter_lingliu_DNAdamage/data/Liu_2016/ChIPseq/RM/,~/Projects/biter_lingliu_DNAdamage/data/Liu_2016/ChIPseq2/RM/ [selchr]=chr10 [gsize]=mm [outdir]=Liu_2016_ChIPseq_ChIPseq2 )
#declare -A o=( [jobnames]=computeGCBias [queue]=normal [t]=10 [nc]=8 [indir]=~/Projects/biter_lingliu_DNAdamage/data/Liu_2016/ChIPseq/RM/,~/Projects/biter_lingliu_DNAdamage/data/Liu_2016/ChIPseq2/RM/ [selchr]=chr10 [gsize]=mm [outdir]=Liu_2016_ChIPseq_ChIPseq2 )
#declare -A o=( [jobnames]=multiBamSummary,plotFingerprint,plotCoverage,plotProfile [queue]=normal [t]=6 [nc]=8 [indir]=~/Projects/biter_biter_shared/results/2016-05-04/Liu_2016_ChIPseq_ChIPseq2 [selchr]=chr10 [gsize]=mm [outdir]=Liu_2016_ChIPseq_ChIPseq2_wGCbiasCorrected )
#declare -A o=( [jobnames]=multiBamSummary,computeGCBias,plotFingerprint,plotCoverage,plotProfile [queue]=normal [t]=10 [nc]=8 [indir]=~/Projects/biter_biter_EWSR1/data/Bilen_2015/ChIPseq/RM/ [selchr]=chr10 [gsize]=hs [gnm]=hg38 [outdir]=Bilen_2015_ChIPseq )


