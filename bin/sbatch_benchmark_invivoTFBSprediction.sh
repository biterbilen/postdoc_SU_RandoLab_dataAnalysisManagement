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

ivTFBS() {
	local analysisdirbase=$1
	local indir=${o[indir]:-indir}; mkdir -p $indir;
	ln -sf ${o[datadir]}/ChIPseq $indir/.
	ln -sf ${o[datadir]}/RNAseq $indir/.
	ln -sf ${o[datadir]}/DNAse $indir/.
	ln -sf ${o[datadir]}/labels $indir/.
	ln -sf ${o[datadir]}/misc $indir/.
	local tfnamei=$(auxGetIndexFromColName metadata ${o[TFName]:-TFName})
	local celltypei=$(auxGetIndexFromColName metadata ${o[CellType]:-CellType})
	local TFNames=($( grep ${o[TFSelPat]:-train } metadata | cut -f $tfnamei | sort | uniq ))
	local cellNames=($( grep ${o[TFSelPat]:-train } metadata | cut -f $celltypei | sort | uniq ))

	local qcmd
	local jobname

	jobname=PWMPrep
	qcmd=$(cat <<- END
		mkdir -p $jobname; pushd $jobname

		# grep -Piw "CEBPB|CTCF|FOSL1|MAX|MYOD1|MYOG|USF1"
		# TODO remove licenced PWMs???

		ln -sf ~/PI_HOME/Data/casco/CIS-DB/${o[orgLatin]:-mm10}/TF_Information.txt .
		ln -sf ~/PI_HOME/Data/casco/CIS-DB/${o[orgLatin]:-mm10}/all_pwms.txt .

		# PWMs with licence
		ls ~/PI_HOME/Data/casco/CIS-DB/${o[orgLatin]}/logos_all_motifs/ | grep _fwd.png | sed 's/_fwd.png//g' | sort | uniq > Motif_ID.wlogos
		grep -f Motif_ID.wlogos TF_Information.txt | \
			awk 'BEGIN{OFS="\t"; print "TFName", "PWMids", "PWMsources"; } \
			\$4 != "." && NR>1 { OFS="\t"; print \$7,\$4,\$15 }' | \
			sort -k 1,1 | awk '{ OFS="\t"; print \$2,\$1,\$3}' > TF.all
		bedtools groupby -i TF.all -g 2 -o collapse,collapse -c 1,3 > TF.grouped.all

		# selected PWMs currently interested in
		awk -v i=$tfnamei 'NR>1{ print \$i }' ../metadata | sort | uniq > tfpatternfile 
		grep -w -f tfpatternfile TF.all > TF.selected
		grep -w -f tfpatternfile TF.grouped.all > TF.grouped.selected

		# PWMs representing more than one TF
		echo "ERROR Start - Any IDs representing more than one TF type?"
		less TF.selected  | cut -f 1 | sort | uniq -c | sort -k 1,1gr | awk '\$1>1{ print }'
		echo "ERROR End - Any IDs representing more than one TF type?"

		# PWMs with matrix 
		Rscript -e "source('~/Projects/biter_biter_shared/bin/auxTFBSTools.R')" \
			-e "getMatrixList(src='all_pwms.txt', gnm='${o[org]:-mm10}', matrixtype='PWM', ii=1, ni=1, byrow=T, select='TF.selected', ofle='CIS-DB.selected')"

		# some stats
		wc -l TF_Information.txt Motif_ID.wlogos TF.all TF.grouped.all tfpatternfile TF.selected TF.grouped.selected

		popd
		END
		)
	if [[ ,${o[jobnames]}, =~ ,$jobname, ]]; then qcmds+=( [$analysisdirbase$jobname]="$qcmd" ); fi

	# TODO 
	# calculate shape for each 200-bp window
	# current time and memory limit might be insufficient
	# each sequence in a file seems to be evaluated independently
	# the software looks for shapes in pentamers; what about frame shifts?
	jobname=DNAshapes
	qcmd=$(cat <<- END
		mkdir -p $jobname; pushd $jobname
		less ~/PI_HOME/Data/casco/UCSC_tracks/${o[org]:-mm10}/assembly/${o[selchrom]:-chr19}.fa* > tmp.fa
		Rscript -e 'library(DNAshapeR);' \
			-e 'pred <- getShape("tmp.fa"); trackShape("tmp.fa", pred); feat <- c("1-mer", "1-shape");' \
			-e 'save(pred, file="pred.rda"); sessionInfo();'
		rm tmp.fa
		popd
		END
		)
	if [[ ,${o[jobnames]}, =~ ,$jobname, ]]; then qcmds+=( [$analysisdirbase$jobname]="$qcmd" ); fi


	# Check rate of ChIP-seq peaks in PWM mapping sites
	#homer2 gave segmentation fault!!! using fimo instead
	jobname=PWMSearch
	for TF in ${TFNames[@]}; do
		qcmd=$(cat <<- END
			mkdir -p $jobname; pushd $jobname
			less ~/PI_HOME/Data/casco/UCSC_tracks/${o[org]:-mm10}/assembly/${o[selchrom]:-chr19}.fa* > tmp$TF.fa
			less ../PWMPrep/CIS-DB.selected | grep -w $TF | awk '{ print \$2; }' | \
				while read m; do
					echo Doing $TF \$m
					fimo --thresh ${o[thres]:-1-e5} --max-stored-scores ${o[mss]:-3000000} --motif \$m --oc \$m ../PWMPrep/CIS-DB.selected tmp$TF.fa 
					# TODO fix FIMOs memory issue and remove break
					break
				done
			echo fimo DONE
			rm tmp$TF.fa
			fimo --version
			popd
		END
		)
		if [[ ,${o[jobnames]}, =~ ,$jobname, ]]; then dqcmds+=( [$analysisdirbase$jobname$TF]="$qcmd" ); fi
	done

	jobname=PWMChIPseqOverlap
	for TF in ${TFNames[@]}; do
		qcmd=$(cat <<- END
			mkdir -p $jobname; pushd $jobname
			less ../PWMPrep/CIS-DB.selected | grep -w $TF | awk '{ print \$2; }' | \
				while read m; do
					PWMf=\$m.bed.gz
					ChIPseqfiles=(\$(ls ../$indir/ChIPseq/*.$TF.conservative.train.narrowPeak.gz))
					less ../PWMSearch/\$m/fimo.txt | awk 'BEGIN{OFS="\t"} NR>1 { print \$2,\$3-1,\$4,\$1,\$6,\$5}' | \
						sort -k 1,1 -k 2,2g | gzip -c > \$PWMf
					s1=\$(less \$PWMf | wc -l)
					for c in \${ChIPseqfiles[@]}; do
						s2=\$(less \$c | grep -w ${o[selchrom]:-chr19} | wc -l)
						i=\$(less \$c | grep -w ${o[selchrom]:-chr19} | bedtools intersect -u -a stdin -b \$PWMf | wc -l)
						echo \$i \$s2 \$s1 | awk -v TF=$TF -v tag="$jobname\t\$(basename \$c)\t\$m" \
							'{ printf("%s\t%.2f\t%.f\t%.2f\t%.f\t%s\n", tag,\$1/\$2,\$2,\$1/\$3,\$3,TF) }'
					done
					# TODO fix FIMOs memory issue and remove break
					break
				done
			popd
		END
		)
		if [[ ,${o[jobnames]}, =~ ,$jobname, ]]; then dqcmds2+=( [$analysisdirbase$jobname$TF]="$qcmd" ); fi
	done
		# TODO afterstats
#		tag=PWMChIPseqOverlap
#echo "tag fileName PWMid peakOverlapRate peakCount PWMOverlapRate PWMcount TF" | sed 's/ /\t/g' > $tag.stat 
#grep "^$tag" slurm_ivTFBS_randomrReprPWM_hg19_Synapse_syn6131484$tag* | awk -F ":" '{print $2}' >> $tag.stat 

	# TODO
	jobname=ChIPseqDNAAccOverlap
	for cell in ${cellNames[@]}; do
		qcmd=$(cat <<- END
			mkdir -p $jobname; pushd $jobname
			DNAsefiles=(\$(ls ../$indir/DNAse/*.$cell.conservative.narrowPeak.gz))
			ChIPseqfiles=(\$(ls ../$indir/ChIPseq/*.$cell.*.conservative.train.narrowPeak.gz))
			for d in \${DNAsefiles[@]}; do
				s1=\$(less \$d | grep -w ${o[selchrom]:-chr19} | wc -l)
				for c in \${ChIPseqfiles[@]}; do
					s2=\$(less \$c | grep -w ${o[selchrom]:-chr19} | wc -l)
					i=\$(less \$c | grep -w ${o[selchrom]:-chr19} | bedtools intersect -u -a stdin -b \$d | wc -l)
					echo \$i \$s2 \$s1 | awk -v cell=$cell -v tag="$jobname\t\$(basename \$c)\t\$(basename \$d)" \
						'{ printf("%s\t%.2f\t%.f\t%.2f\t%.f\t%s\n", tag,\$1/\$2,\$2,\$1/\$3,\$3,cell) }'
				done
			done
			popd
		END
		)
		if [[ ,${o[jobnames]}, =~ ,$jobname, ]]; then dqcmds2+=( [$analysisdirbase$jobname$cell]="$qcmd" ); fi
	done
		# TODO afterstats
#tag=ChIPseqDNAAccOverlap
#echo "tag ChIPseqFileName DNAseFileName ChIPseqPeakOverlapRate ChIPseqPeakCount DNAsePeakOverlapRate DNAsePeakCount cell" | sed 's/ /\t/g' > $tag.stat
#grep "^$tag" slurm_ivTFBS_randomrReprPWM_hg19_Synapse_syn6131484$tag* | awk -F ":" '{print $2}' >> $tag.stat 

	jobname=TFexpressionRank
	qcmd=$(cat <<- END
		mkdir -p $jobname; pushd $jobname
		ln -sf ~/PI_HOME/Data/casco/Gencode/gencode.v19.chr_patch_hapl_scaff.annotation.gtf.gz .
		less ../metadata | grep train | cut -f 2 | sort | uniq > a
		# Extract gene_id and gene_name
		less gencode.v19.chr_patch_hapl_scaff.annotation.gtf.gz | grep -v exon_number | \
			cut -f 9- | awk '{ OFS="\t"; print \$2,\$10}' | sed 's/";*//g' | sort | uniq > TFNameMap.all 
		grep -w -f a TFNameMap.all | grep -v "-" > TFNameMap.selected
		cut -f 1 TFNameMap.selected > id
		# Rank expression based of TPM
		for f in ../indir/RNAseq/*; do 
			b=\$(basename \$f .tsv)
			grep -v gene_id \$f | sort -k 6,6gr | awk 'BEGIN{i=1; OFS="\t"; }{ print \$1,i++; }' | > \$b.rank 
			grep -w -f id \$b.rank | awk -v tag="$jobname\t\$b" 'BEGIN{OFS="\t";} {print tag,\$0}'
		done
		grep ^$jobname *.rank | sort -k 3,3 > stat
		rm a id
	END
	)
	if [[ ,${o[jobnames]}, =~ ,$jobname, ]]; then dqcmds2+=( [$analysisdirbase$jobname]="$qcmd" ); fi
		# TODO afterstats
#tag=TFexpressionRank
#less slurm_ivTFBS_randomrReprPWM_hg19_Synapse_syn6131484$tag.out | grep ^$tag | sort -k 3,3 > tmp$tag.stat 
#echo "tag RNAseqFileName TFexpressionRank TF" | sed 's/ /\t/g' >> $tag.stat
#join tmp$tag.stat $tag/TFNameMap.selected -1 3 -2 1 -t $'\t' | cut -f 2- > $tag.stat
#rm tmp$tag.stat
	# TODO  isChIPRateInDNAseAnticorrelatedwTFExpression
#	r=$RANDOM
#	less ChIPseqDNAAccOverlap.stat | sed 's/ChIPseq\.//g' | sed 's/.conservative.train.narrowPeak.gz//g' | sort -k 2,2 > a.$r
#	less TFexpressionRank.stat | sed 's/gene_expression.//g' | sed 's/\./\t/' | awk '{OFS="\t"; print $1,$2"."$5,$4,$3}' | sort -k 2,2 > b.$r
#	echo "tag minExpressionRank maxExpressionRank medianExpressionRank" | sed 's/ /\t/g' > $jobname.stat2
#	join a.$r b.$r -1 2 -2 2 -t $'\t' | sort -k 4,4g | head -n 80 | bedtools groupby -g 2 -c 10,10,10 -o min,max,median | sed 's/Overlap/LowOverlap/g' >> $jobname.stat2
#	join a.$r b.$r -1 2 -2 2 -t $'\t' | sort -k 4,4g | head -n 160 | tail -n 80 | bedtools groupby -g 2 -c 10,10,10 -o min,max,median | sed 's/Overlap/IntermediateOverlap/g' >> $jobname.stat2
#	join a.$r b.$r -1 2 -2 2 -t $'\t' | sort -k 4,4g | tail -n 80 | bedtools groupby -g 2 -c 10,10,10 -o min,max,median | sed 's/Overlap/HighOverlap/g' >> $jobname.stat2
#	rm a.$r b.$r



	# TODO baseline with H2O
	# https://github.com/nboley/DREAM_invivo_tf_binding_prediction_challenge_baseline
	# start from plot_logistic_l1_l2_sparsity.py


}

#For 3M instances on chr10, nc=3 tmin=5 (time limit might be critical?)
#For 1M instances on chr10, nc=1 
#slurmMain ivTFBS --jobnames=PWMPrep,PWMSearch,PWMChIPseqOverlap,ChIPseqDNAAccOverlap,TFexpressionRank --datadir=~/PI_HOME/Data/casco/ENCODE/hg19/Synapse_syn6131484 --otag=randomrReprPWM --useTFName=T --nc=3 --t=0  --tmin=5 \
slurmMain ivTFBS --jobnames=DNAshapes --datadir=~/PI_HOME/Data/casco/ENCODE/hg19/Synapse_syn6131484 --otag=randomrReprPWM --useTFName=T \
	--orgLatin=Homo_sapiens --org=hg19 --selchrom=chr10 --nc=3 --t=1  --tmin=5

#declare -A o=( [jobnames]=liftOver2mm10 [queue]=normal [t]=1 [nc]=1 [selchrom]=chr19 [outdir]=datachr19 [over]=~/PI_HOME/Data/casco/UCSC_tracks/mm10/mm9ToMm10.over.chain.gz)
#declare -A o=( [jobnames]=fimo [queue]=normal [t]=1 [nc]=1 [selchrom]=chr19 [outdir]=fimochr19 [msc]=1000000 [thres]=1-e5 [org]=mm10)
#declare -A o=( [jobnames]=chipOverlap [queue]=normal [t]=1 [nc]=1 [selchrom]=chr19 [outdir]=overlap)
#declare -A o=( [jobnames]=footprint [queue]=normal [t]=1 [nc]=1 [selchrom]=chr19 [outdir]=footprint [minreadcut]=20 )
#declare -A o=( [jobnames]=PWMPrep [queue]=normal [t]=1 [nc]=1 [outdir]=cisdb [orgLatin]=Mus_musculus)

#declare -A o=( [jobnames]=PWMPrep [queue]=normal [t]=1 [nc]=1 [outdir]=cisdb [orgLatin]=Homo_sapiens)
#declare -A o=( [jobnames]=fimo [queue]=normal [t]=4 [nc]=1 [selchrom]=chr19 [outdir]=fimochr19 [msc]=1000000 [thres]=1-e5 [org]=hg19)
