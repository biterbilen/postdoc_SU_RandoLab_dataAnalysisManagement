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
	name=$1
	selchrom=$2
	ofle=$3

	samtools view -@ ${o[nc]} $indir/$name.sorted.bam $selchr | \
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
		# calculate similarity
		# non-overlapping windows of 1K with >20 reads
		local ofleTag=similarityW${mwW}S${mwS}MRC${minReadCutoff}
		local genomef=$HOME/aux/genomes/${genomel:-mouse.mm10}.genome
		grep -w $selchrom $genomef > genome.$selchrom
		bedtools makewindows -g genome.$selchrom -w ${o[makewindowsW]:-1000} -s ${o[makewindowsS]:-1000} > genome.$selchrom.bed
		echo Chr Start End ${names[@]}.sorted.bam | sed 's/ /\t/g' > genome.$selchrom.bedcov
		samtools bedcov genome.$selchrom.bed $indir/${names[@]}.sorted.bam >> genome.$selchrom.bedcov
		gzip -f genome.$selchrom.bedcov
		Rscript -e "source('$PI_SCRATCH/Projects/biter_biter_shared/bin/auxML.R')" \
			-e "param.preprocess <- list(col.sel.pat='.sorted.bam$', colname.split.rank=2, colname.split.pat='.')" \
			-e "spike <- list(spikeId='ERCC', rm=T)" \
			-e "param.process <- list(method='filterRowMeanMin', cutoff=${o[minReadCutoff]:-20})" \
			-e "dta.exp <- caller.getExpressed('genome.$selchrom.bedcov.gz', param.preprocess=param.preprocess, ofle.tag='$ofleTag', nc=${o[nc]}, spike=spike, density.norm.colname=NULL, param.process=param.process)"

		# calculate nuc spacing
		local ofleTag=nucSpacing
		local r=$RANDOM
		parallel --delay 1 --plus --tmpdir . -j0 --header : --tag --xapply \
			"auxCalcNucSpacing {name} $selchrom > $r.{name}" \
			::: name ${names[@]##*/}
		cat $r.* > $ofleTag.txt 
		Rscript -e "source('$PI_SCRATCH/Projects/biter_biter_shared/bin/auxggplot2.R')" \
						-e "a <- read.table('"$ofleTag.txt"')" \
						-e "ggplot(subset(a,V1<2000), aes(x=V1, y=V2, color=V3)) + geom_line()"
					  -e "ggsave('"$ofleTag.pdf"')"
		rm $r.*
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
}
#	# predictd
#	macs2 predictd -g ${o[gsize]:-mm} -m 5 50 -i $indir/$name.sorted.bam --rfile $name\_predictd.R; 
#	Rscript $name\_predictd.R

	# filterdup
#	macs2 filterdup -g ${o[gsize]:-mm} --keep-dup=1 -i $indir/$name$ft.sorted.bam -o $name.bed; 
#	gzip -f $name.bed


PFdump() {
	local analysisdirbase=$1
	local indir=${o[indir]:-RM}

	# treatment and control data tags 
	# form associative array from string
	# http://stackoverflow.com/questions/6660010/bash-how-to-assign-an-associative-array-to-another-variable-name-e-g-rename-t
	local ttagsworcsS="$(bioawk -c hdr 'NR>1 && $Cname != "NA" { ORS=" "; gsub("_[0-9]*$","",$Name); gsub("_[0-9]*$","",$Cname); print "["$Name"]="$Cname; }' metadata)"
	local ttagscsS="$(bioawk -c hdr 'NR>1 && $Cname != "NA"{ ORS=" "; print "["$Name"]="$Cname; }' metadata)"
	local ttags=($(bioawk -c hdr 'NR>1 && $Pair != 2 && $Cname != "NA"{ print $Name}' metadata))
	local ctags=($(bioawk -c hdr 'NR>1 && $Pair != 2 && $Cname == "NA"{ print $Name}' metadata))
	local tags=(${ttags[@]} ${ctags[@]})
	# replace pattern of replicate information
	local ttagswor=($(echo ${ttags[@]/%_[0-9]/} | sed 's/ /\n/g' | sort | uniq))
	local ctagswor=($(echo ${ctags[@]/%_[0-9]/} | sed 's/ /\n/g' | sort | uniq))
	local tagswor=(${ttagswor[@]} ${ctagswor[@]})

	ln -sf ${o[datadir]}/$indir $indir

	local qcmd
	qcmd=$(cat <<- END
		jobnames=${o[jobnames]:-NA}
		selchrom=${o[selchrom]:-chr19}
		ttagsworcsS="$ttagsworcsS"
		ttagscsS="$ttagscsS"
		ttags=(${ttags[@]})
		ctags=(${ctags[@]})
		tags=(${tags[@]})
		ttagswor=(${ttagswor[@]})
		ctagswor=(${ctagswor[@]})
		tagswor=(${tagswor[@]})

		indir=$indir
		selchrom=$selchrom
		genomef=~/aux/genomes/${o[genomel]:-mouse.mm10}.genome
		mwW=${o[makewindowsW]:-1000}
		mwS=${o[makewindowsS]:-1000}
		minReadCutoff=${o[minReadCutoff]:-20}
		ofleTag=similarityW\${mwW}S\${mwS}MRC\${minReadCutoff}
		nc=${o[nc]}
		if [[ ,${o[jobnames]}, =~ ,calcsim, ]]; then 
			$(declare -f calcsim)
			calcsim $indir $mwW $mwS $minReadCutoff $selchrom $genomef $ofleTag $nc
		fi

		END
		)

  jobnames+=($analysisdirbase)
  qcmds[$analysisdirbase]="$qcmd"
  #qarrs[$analysisdirbase]=$(( ${#names[@]} - 1 ))
  qarrs[$analysisdirbase]=0
  qwpts[$analysisdirbase]=F

	#if [[ ,${o[jobnames]}, =~ ,$jobname, ]]; then qcmds+=( [$analysisdirbase$jobname]="$qcmd" ); fi

	jobname=predictd; 
	for FT in ${tags[@]}; do
		qcmd=$(cat <<- END
			mkdir -p $jobname; pushd $jobname
			source ~/PI_HOME/PythonSandbox/MACS2/bin/activate
			macs2 predictd -g ${o[gsize]:-mm} -m 5 50 -i ../$indir/$FT.sorted.bam --rfile $FT\_predictd.R; 
			Rscript $FT\_predictd.R
			deactivate
			popd
			END
			)
		if [[ ,${o[jobnames]}, =~ ,$jobname, ]]; then qcmds+=( [$analysisdirbase$jobname$FT]="$qcmd" ); fi
	done

	jobname=filterdup
	if [ "${o[otag]:-}" == poolreps ]; then
		A=${tagswor[@]}
		ft="_*"
	else
		A=${tags[@]}
		ft=""
	fi
	for FT in ${A[@]}; do
		qcmd=$(cat <<- END
			mkdir -p $jobname; pushd $jobname
			macs2 filterdup -g ${o[gsize]:-mm} --keep-dup=1 -i ../$indir/$FT$ft.sorted.bam -o $FT.bed; 
			gzip -f $FT.bed
		popd
		END
		)
		if [[ ,${o[jobnames]}, =~ ,$jobname, ]]; then qcmds+=( [$analysisdirbase$jobname$FT]="$qcmd" ); fi
	done

	jobname=pileup
	if [ "${o[otag]:-}" == poolreps ]; then
		A=${ttagswor[@]}
		B=${ctagswor[@]}
	else
		A=${ttags[@]}
		B=${ctags[@]}
	fi
	for FT in ${A[@]}; do
		qcmd=$(cat <<- END
			mkdir -p $jobname; pushd $jobname
			macs2 pileup -i ../filterdup/$FT.bed.gz --extsize ${o[extsize]:-73} -o tmp$FT.bdg
			bedClip tmp$FT.bdg ~/aux/genomes/${o[genomel]:-mouse.mm10}.genome $FT.bdg
			rm tmp$FT.bdg
		popd
		END
		)
		if [[ ,${o[jobnames]}, =~ ,$jobname, ]]; then dqcmds+=( [$analysisdirbase$jobname$FT]="$qcmd" ); fi
	done
	for FT in ${B[@]}; do
		qcmd=$(cat <<- END
			mkdir -p $jobname; pushd $jobname

			if [[ ${o[gsize]:-mm} == mm ]]; then
				gsize=1870000000
			elif [[ ${o[gsize]:-mm} == hs ]]; then
				gsize=2700000000
			else
				echo "ERROR: gsize ${o[gsize]:-mm} is not defined"
			fi
			bedf=../filterdup/$FT.bed.gz

			# pileup
			macs2 pileup -B -i \$bedf --extsize ${o[extsize]:-73} -o $FT.bdg
			macs2 pileup -B -i \$bedf --extsize ${o[slocal]:-500} -o ${FT}_1k_bg.bdg
			macs2 pileup -B -i \$bedf --extsize ${o[llocal]:-5000} -o ${FT}_10k_bg.bdg

			# local bias raw
			mslocal=\$(echo ${o[extsize]:-73} / ${o[slocal]:-500} | bc -l)
			macs2 bdgopt -i ${FT}_1k_bg.bdg -m multiply -p \$mslocal -o ${FT}_1k_bg_norm.bdg

			mllocal=\$(echo ${o[extsize]:-73} / ${o[llocal]:-5000} | bc -l)
			macs2 bdgopt -i ${FT}_10k_bg.bdg -m multiply -p \$mllocal -o ${FT}_10k_bg_norm.bdg

			macs2 bdgcmp -m max -t ${FT}_1k_bg_norm.bdg -c ${FT}_10k_bg_norm.bdg -o ${FT}_1k_10k_bg_norm.bdg
			macs2 bdgcmp -m max -t ${FT}_1k_10k_bg_norm.bdg -c ${FT}.bdg -o ${FT}_d_1k_10k_bg_norm.bdg

			gb=\$(grep "tags after filtering in alignment file" ../slurm*filterdup${FT}.out | \
				awk -v gsize=\$gsize -v extsize=${o[extsize]:-73} '{ print \$NF*2*extsize/gsize}')
			macs2 bdgopt -i ${FT}_d_1k_10k_bg_norm.bdg -m max -p \$gb -o tmp${FT}_local_bias_raw.bdg
			less tmp${FT}_local_bias_raw.bdg | grep -v "^track" > tmp2${FT}_local_bias_raw.bdg
			bedClip tmp2${FT}_local_bias_raw.bdg ~/aux/genomes/${o[genomel]:-mouse.mm10}.genome ${FT}_local_bias_raw.bdg
			rm tmp${FT}_local_bias_raw.bdg tmp2${FT}_local_bias_raw.bdg 
		popd
		END
		)
		if [[ ,${o[jobnames]}, =~ ,$jobname, ]]; then dqcmds+=( [$analysisdirbase$jobname$FT]="$qcmd" ); fi
	done

	jobname=locallambda
	if [ "${o[otag]:-}" == poolreps ]; then
		A=${ttagswor[@]}
		eval "local -A B=("${ttagsworcsS}")"
	else
		A=${ttags[@]}
		eval "local -A B=("${ttagscsS}")"
	fi
	for FT in ${A[@]}; do
		qcmd=$(cat <<- END
			mkdir -p $jobname; pushd $jobname
			c=${B[$FT]}
			t1=../pileup/$FT.bdg
			lsp="tags after filtering in alignment file"
			d=\$(grep "\$lsp" ../slurm*filterdup$FT.out | awk '{ print \$NF}')
			k=\$(grep "\$lsp" ../slurm*filterdup\$c.out | awk -v d=\$d '{ print d/\$NF}')

			macs2 bdgopt -m multiply -p \$k -i ../pileup/\$c\_local_bias_raw.bdg -o $FT\_\$c\_local_lambda.bdg
		END
		)
		if [[ ,${o[jobnames]}, =~ ,$jobname, ]]; then dqcmds2+=( [$analysisdirbase$jobname$FT]="$qcmd" ); fi
	done

	# replace all occurances of , w/ " " using //
	local t1s=(${o[t1]//,/ })
	if [ ${o[t2]:-X} == X ]; then
		jobname=callpeak
		if [ "${o[otag]:-}" == poolreps ]; then
			eval "local -A A=("${ttagsworcsS}")"
		else
			eval "local -A A=("${ttagscsS}")"
		fi
		for FT in ${t1s[@]}; do
			qcmd=$(cat <<- END
				mkdir -p $jobname; pushd $jobname
				macs2 bdgcmp -t ../pileup/$FT.bdg -c ../locallambda/$FT*.bdg -m qpois -o $FT\_qvalue.bdg
				macs2 bdgpeakcall -i $FT\_pvalue.bdg -c ${o[pvaluecutof]:-3} -g ${o[g]:-100} -l ${o[l]:-146} -o $FT\_qvalue.bed
				gzip -f $FT\_qvalue.bed
				macs2 bdgcmp -t ../pileup/$FT.bdg -c ../locallambda/$FT*.bdg -m ppois -o $FT\_pvalue.bdg
				macs2 bdgpeakcall -i $FT\_pvalue.bdg -c ${o[pvaluecutof]:-3} -g ${o[g]:-100} -l ${o[l]:-146} -o $FT\_pvalue.bed
				gzip -f $FT\_pvalue.bed
				macs2 bdgcmp -t ../pileup/$FT.bdg -c ../locallambda/$FT*.bdg -m FE -o $FT\_FE.bdg
				macs2 bdgpeakcall -i $FT\_FE.bdg -c ${o[FEcutoff]:-3} -g ${o[g]:-100} -l ${o[l]:-146} -o $FT\_FE.bed
				gzip -f $FT\_FE.bed

				n=nomodel$FT
				macs2 callpeak --nomodel --shift ${o[shift]:-0} --extsize ${o[extsize]:-73} --trackline --SPMR -B -g ${o[gsize]:-mm} -n \$n \
					-c ../filterdup/${A[$FT]}.bed.gz -t ../filterdup/$FT.bed.gz
				n=model$FT
				macs2 callpeak --shift ${o[shift]:-0} --extsize ${o[extsize]:-73} --trackline --SPMR -B -g ${o[gsize]:-mm} -n \$n \
					-c ../filterdup/${A[$FT]}.bed.gz -t ../filterdup/$FT.bed.gz
				Rscript \$n\_model.r
			END
			)
			if [[ ,${o[jobnames]}, =~ ,$jobname, ]]; then dqcmds3+=( [$analysisdirbase$jobname$FT]="$qcmd" ); fi
		done
	else
		jobname=bdgdiff
		# TODO do custm bdgdff
		# replace all occurances of , w/ " " using //
		local t2s=(${o[t2]//,/ })
		for FTi in ${!t1s[@]}; do
			FT1=${t1s[$FTi]}
			FT2=${t2s[$FTi]}
			qcmd=$(cat <<- END
				mkdir -p $jobname; pushd $jobname
				FT1=$FT1
				FT2=$FT2
				lsp="tags after filtering in alignment file"
				d1=\$(grep "\$lsp" ../slurm*filterdup\$FT1.out | awk '{ print \$NF}')
				d2=\$(grep "\$lsp" ../slurm*filterdup\$FT2.out | awk '{ print \$NF}')
				c1=../locallambda/\$FT1*.bdg
				c2=../locallambda/\$FT2*.bdg
				t1=../pileup/\$FT1.bdg
				t2=../pileup/\$FT2.bdg
				macs2 bdgdiff --t1 \$t1 --c1 \$c1 --t2 \$t2 --c2 \$c2 --d1 \$d1 --d2 \$d2 --cutoff ${o[cutoff]:-3} -g ${o[g]:-100} -l ${o[l]:-146} --o-prefix diff_\$FT1\_vs_\$FT2
				gzip -f *.bed
				#TODO customize with bdgcmp???
			END
			)
			if [[ ,${o[jobnames]}, =~ ,$jobname, ]]; then dqcmds3+=( [$analysisdirbase$jobname$FT1\vs$FT2]="$qcmd" ); fi
		done
	fi

}

# 4 Ling - Hairless paper submission
slurmMain PF --datadir=$PI_SCRATCH/Projects/biter_lingliu_DNAdamage/data/Liu_2017/ATACseq --indir=RM --t1=HrCKO_R_SC_1,HrCKO_R_SC_2,HrCKO_R_SC_3,HrCKO_S_SC_1 --t2=HrWT_R_SC_1,HrWT_R_SC_2,HrWT_R_SC_3,HrWT_S_SC_1 
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


