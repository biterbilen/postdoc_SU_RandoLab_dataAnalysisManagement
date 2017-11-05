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

PF() {
	local analysisdirbase=$1
	local indir=${o[indir]:-RM}; mkdir -p $indir;
	ln -sf ${o[datadir]}/$indir/*.sorted.bam $indir/.

	# indexes of name and cname
	local namei=$(auxGetIndexFromColName metadata ${o[name]:-Name})
	local cnamei=$(auxGetIndexFromColName metadata ${o[name]:-Cname})
	# treatment and control data tags 
	local -a ttagscs=($(awk -v namei=$namei -v cnamei=$cnamei 'NR>1 && $cnamei!="NA"{ print $cnamei}' metadata))
	local -a ttags=($(awk -v namei=$namei -v cnamei=$cnamei 'NR>1 && $cnamei!="NA"{ print $namei}' metadata))
	local -a ctags=($(awk -v namei=$namei -v cnamei=$cnamei 'NR>1 && $cnamei=="NA"{ print $namei}' metadata))
	local -a tags=(${ttags[@]} ${ctags[@]})
	local -a ttagswor=($(echo ${ttags[@]} | sed 's/_[0-9]*$//g' | sort | uniq))
	local -a ctagswor=($(echo ${ctags[@]} | sed 's/_[0-9]*$//g' | sort | uniq))
	local -a tagswor=(${ttagswor[@]} ${ctagswor[@]})

	local qcmd
	local jobname

	jobname=bdgdiff
	for FT in ${stags[@]}; do
		qcmd=$(cat <<- END
			mkdir -p $jobname; pushd $jobname

			# set control samples of the treatment files
			o[c1]=\$(grep -w \$t1 ../metadata | cut -f $cnamei)
			o[c2]=\$(grep -w \$t2 ../metadata | cut -f $cnamei)

			t1=../pileup/${o[t1]:-X}.bdg
			t2=../pileup/${o[t2]:-X}.bdg
			lsp="tags after filtering in alignment file"
			d1=\$(grep "\$lsp" ../slurm*filterdup${o[t1]:-X}.out | awk '{ print \$NF}')
			d2=\$(grep "\$lsp" ../slurm*filterdup${o[t2]:-X}.out | awk '{ print \$NF}')
			k1=\$(grep "\$lsp" ../slurm*filterdup${o[c1]:-X}.out | awk -v d=\$d1 '{ print d/\$NF}')
			k2=\$(grep "\$lsp" ../slurm*filterdup${o[c2]:-X}.out | awk -v d=\$d2 '{ print d/\$NF}')
			clbr1=../pileup/${o[c1]}_local_bias_raw.bdg 
			clbr2=../pileup/${o[c2]}_local_bias_raw.bdg 
			c1=../pileup/${o[t1]:-X}_${o[c1]}_local_lambda.bdg 
			c2=../pileup/${o[t2]:-X}_${o[c2]}_local_lambda.bdg

			macs2 bdgopt -i \$clbr1 -m multiply -p \$k1 -o \$c1
			macs2 bdgopt -i \$clbr2 -m multiply -p \$k2 -o \$c2
			
			macs2 bdgdiff --cutoff ${o[cutoff]:-999} --t1 \$t1 --c1 \$c1 --t2 \$t2 --c2 \$c2 --d1 \$d1 --d2 \$d2 -g ${o[g]:-100} -l ${o[l]:-146} --o-prefix diff_${o[t1]}_vs_${o[t2]}
			gzip *.bed

			popd
		END
		)
		if [[ ,${o[jobnames]}, =~ ,$jobname, ]]; then qcmds+=( [$analysisdirbase$jobname$FT]="$qcmd" ); fi
	done

	jobname=GB
	qcmd=$(cat <<- END
		bb_scaleTruncBedScore() {
			# Scales scores to 1-1000 proportionally; like in Thales\' theorem
			local bedp=\$1
			local bedt=\$(basename \$bedp | awk -F . '{ print $1}')
			local min=\$(less \$bedp | grep ^chr | cut -f 5 | sort -g | head -n 1)
			local max=\$(less \$bedp | grep ^chr | cut -f 5 | sort -gr| head -n 1)
			less \$bedp | cut -f 1-6 | awk -v min=\$min -v max=\$max '\$1 ~ /chr/ { OFS="\t"; \$5=int((1000-1)*(\$5-min)/(max-min)+1); print \$0;}' > tmp\$bedt
			bedToBigBed tmp\$bedt ~/aux/genomes/mouse.${o[gnm]:-mm10}.genome \$bedt.bb; 
			rm tmp\$bedt
		}

		#export -f bb_scaleTruncBedScore
		#parallel --tmpdir . 'bb_scaleTruncBedScore {}' ::: ../bdgdiff/*bed.gz

		bb_track() {
			local d=\$1
			local ind=\$2

			local bigDataUrl="http://www-dev.stanford.edu/group/rando_lab/biter/\$d"
			mkdir -p \$d
			cat > \$d/hub.txt <<- END1
				hub \$d 
				shortLabel \$d 
				longLabel \$d 
				genomesFile genomes.txt
				email biter@stanford.edu
			END1
			cat > \$d/genomes.txt <<- END1
				trackDb tracks.txt
				genome mm10 
			END1

			echo browser position chr10:212000-250000 > \$d/tracks.txt
			for s in \$(find \$ind -name *bb 2> /dev/null); do
				f=\$(basename \$s)
				cat >> \$d/tracks.txt <<- END1 

					track \$f
					shortLabel \$f
					longLabel \$f
					type bigBed 3
					bigDataUrl \$bigDataUrl/\$f
					visibility full
					nextItemButton on
				END1
				mv \$s \$d/\$f
			done
			for s in \$(find \$ind -name *treat*bw 2> /dev/null); do
				f=\$(basename \$s)
				cat >> \$d/tracks.txt <<- END1 

					track \$f
					shortLabel \$f
					longLabel \$f
					type bigWig
					bigDataUrl \$bigDataUrl/\$f
					maxHeightPixels 50:30:11
					viewLimits 3:50
					visibility full
					color 0,0,0
				END1
				mv \$s \$d/\$f
			done
		}
		bb_track GB_$analysisdirbase .

			popd
	END
	)
	if [[ ,${o[jobnames]}, =~ ,$jobname, ]]; then qcmds+=( [$analysisdirbase$jobname$FT]="$qcmd" ); fi

	# TODO calculate enrichments
	jobname=annotateRegion
	qcmd=$(cat <<- END
		mkdir -p $jobname; pushd $jobname
		less ~/aux/genomes/mouse.${o[gnm]:-mm10}.genome | awk -v gnm=${o[gnm]:-mm10} '{ print gnm"\t"\$2}' | \
			bedtools groupby -g 1 -c 2 -o sum | awk '{ OFS="\t"; print \$1,int(\$2/1000000)"M"}' > \$r.stat
		annotations=(~/PI_HOME/Data/casco/UCSC_tracks/${o[gnm]:-mm10}/canonical/promoters.bed.gz 
								~/PI_HOME/Data/casco/UCSC_tracks/${o[gnm]:-mm10}/canonical/transcripts_wopromoters.bed.gz)
		for an in \${annotations[@]}; do
			less \$an | sort -k 1,1 -k 2,2g | bedtools merge -i stdin | awk -v tag=\$(basename \$an .bed.gz) '{ print tag"\t"\$3-\$2 }' | \
				bedtools groupby -g 1 -c 2 -o sum | awk '{ OFS="\t"; print \$1,int(\$2/1000000)"M"}'
		done >> $jobname.stat
		for f in ../bdgdiff/*common.bed.gz; do
			echo
			tag=\$(basename \$f _common.bed.gz)
			#parallel -k --delay 1 --plus --tmpdir -j0 --header : --tag --joblog \$r.jlog \
			#	::: TP cond1 cond2 common
			for TP in cond1 cond2 common; do
				less ../bdgdiff/\$tag\_\$TP.bed.gz | grep "^chr" | bedtools annotate -i stdin -files \${annotations[@]} | \
					awk -v tag1=\$tag -v tag2=\$TP '{ r=""; if (\$6>0) { r="promoters"; } if (\$7>0) { r=r"_genic"; } if (r=="") { r="intergenic"}  OFS="\t"; print r, tag1, tag2; }' | \
					sed 's/^_//g' | sort | uniq -c 
				echo
			done
			done >> \$r.stat
			popd
		fi

		# TODO continue here
		r=callpeak
		if [[ ${o[callpeak]:-F} == T ]]; then
			n=${o[stag]}
			macs2 callpeak --nomodel --shift ${o[shift]:-0} --extsize ${o[extsize]:-73} --trackline --SPMR -B -g ${o[gsize]:-mm} -n \$n $ct -t ${samples[@]}; Rscript \$n\_model.r
		fi &>> slurm_\$r.out
		echo \$r

		echo
		Rscript --version
		parallel --version
		macs2 --version

		END
		)
	if [[ ,${o[jobnames]}, =~ ,$jobname, ]]; then qcmds+=( [$analysisdirbase$jobname$FT]="$qcmd" ); fi
}

slurmMain PF --jobnames=calcsim,poolreps,predictd,filterdup,pileup,bdgdiff,GB,annotateRegion,callpeak --datadir=~/Projects/biter_lingliu_DNAdamage/data/Liu_2016/ChIPseq2 --nc=16 --stag=WT \
	  --t1=WT_old_NA_SC_H3K9me2 --t2=WT_old_NA_SC_H3K9me3 --c1=WT_old_NA_SC_NA --c2=WT_old_NA_SC_NA    --GB=$gb --annotateRegion=$ar     --callpeak=$cp

#slurmMain PF --datadir=~/Projects/biter_lingliu_DNAdamage/data/Liu_2016/ChIPseq2 --nc=16 --stag=WT --calcsim=$cs --poolreps=$pr --predictd=$pd --filterdup=$fd --pileup=$pu \
#	  --bdgdiff=$bd --t1=WT_old_NA_SC_H3K9me2 --t2=WT_old_NA_SC_H3K9me3 --c1=WT_old_NA_SC_NA --c2=WT_old_NA_SC_NA    --GB=$gb --annotateRegion=$ar     --callpeak=$cp


exit
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


