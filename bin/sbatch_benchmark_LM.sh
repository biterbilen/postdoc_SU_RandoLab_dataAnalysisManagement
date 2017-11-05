#!/bin/bash - 
#===============================================================================
#
#          FILE: sbatch_map_wSTAR.sh
# 
#         USAGE: ./sbatch_map_wSTAR.sh 
# 
#   DESCRIPTION: 
# 
#       OPTIONS: ---
#  REQUIREMENTS: ---
#          BUGS: ---
#         NOTES: ---
#        AUTHOR: Dr. Biter Bilen (bb), biterbilen@yahoo.com
#  ORGANIZATION: Stanford, Palo Alto, USA
#       CREATED: 04/02/2015 17:53
#      REVISION:  ---
#===============================================================================

set -o nounset                              # Treat unset variables as an error

t=12
nc=16

RTAGs=(denemePromoters)

RTAGs=(allPromoters canPromoters)
MTAGs=(JASPAR2014 JASPARCOREVERT HOMER)


RTAGs=(Acc4yAyQoQwSPMRwSHIFT0 Acc4yAyQoQwSPMRBroadwSHIFT0); n=ATACseq
RTAGs=(Acc4yAyQoQwSPMRwSHIFT0); n=ATACseq_$RTAGs
#RTAGs=(Acc4yAyQoQwSPMRBroadwSHIFT0); n=ATACseq_$RTAGs
MTAGs=(JASPAR2014); n=ATACseq_$RTAGs
BTAGs=(Acc4yAyQoQ); paired=T

odir=LM_featureCount_${n}
mkdir -p $odir; pushd $odir;

#echo $cmd
#callpeak,bdgpeakcall,bdgbroadcall,bdgcmp,bdgopt,cmbreps,bdgdiff,filterdup,predictd,pileup,randsample,refinepeak
#sbatch -J $odir -p normal -t $t:0:0 -c $nc -e slurm_$n.out -o slurm_$n.out <<- END
# Too much load on one node; could be sent to separate nodes!!!
cat > cmds.sh <<- END
	#!/bin/bash -l
	bb_setBAMf() {
		local btag=\$1
		local bamf
		if [[ \$btag == Acc4yAyQoQ ]]; then
			bamf=(\$(ls ~/Projects/biter_lingliu_DNAdamage/data/Liu_2015/ATACseq/RM_TR_ada/*sorted.bam))
		else
			die ERROR: \$btag not known
		fi
		echo \${bamf[@]}
	}

	bb_setRegf() {
		local rtag=\$1
		local regf
		if [[ \$rtag == denemePromoters ]]; then
			regf=~/PI_HOME/Data/casco/UCSC_tracks/mm10/canonical/deneme.bed.gz
		elif [[ \$rtag ==	canPromoters ]]; then
			regf=~/PI_HOME/Data/casco/UCSC_tracks/mm10/canonical/promoters.bed.gz
		elif [[ \$rtag ==	allPromoters ]]; then
			regf=~/Projects/biter_jbrett_epigenetics/results/2015-01-28/promoters_A4Q4R5/allPromoters/TM_allPromoters.tab.gz
		elif [[ \$rtag ==	Acc4yAyQoQwSPMRwSHIFT0 ]]; then
			regf=~/Projects/biter_biter_shared/results/2015-05-19/PF_wSPMR_Acc4yAyQoQ/callpeak/TAG/wSPMR/SHIFT/0/callpeak_peaks.narrowPeak 
		elif [[ \$rtag ==	Acc4yAyQoQwSPMRBroadwSHIFT0 ]]; then
			regf=~/Projects/biter_biter_shared/results/2015-05-19/PF_wSPMRBroad_Acc4yAyQoQ/callpeak/TAG/wSPMRBroad/SHIFT/0/callpeak_peaks.broadPeak
		else
			die ERROR: \$rtag not known
		fi
		less \$regf | awk '\$1 ~ /^chr/ { if (\$6 == ".") { \$6="+"; } OFS="\t"; print \$1,\$2,\$3,\$4,\$5,\$6}' | sort -k 5,5gr | head -n 20000 > \$(basename \$regf).bed
		regf=\$(basename \$regf).bed
		echo \$regf
	}

	XXXbb_associateRegionWPromoters() {
		local rtag=\$1
		local minDirDist=\$2
		local odir=\$3
		local maxDirDist=-100000
		local regf=\$(bb_setRegf \$rtag)

		mkdir -p \$odir; pushd \$odir
		ln -sf \$regf .
		of=\$(basename \$regf .tab.gz)
		less \$regf | awk 'NR>1{ OFS="\t"; print \$1,\$2,\$3}' > tmp\$of.bed
		less ~/PI_HOME/Data/casco/UCSC_tracks/mm10/canonical/promoters.bed.gz | \
			bedtools closest -a stdin -b tmp\$of.bed -D a -t first | \
			awk '\$13 != "." && \$NF<\$minDirDist && \$NF>\$maxDirDist{ OFS="\t"; print \$13,\$14,\$15,\$12}' | sed 's/";*//g' > tmp2\$of.bed
		wc -l tmp2\$of.bed

		echo "R --no-save --args tmp2\$of.bed \$f \$of < a.R &> \$of.log"
		R --no-save --args tmp2\$of.bed \$f \$of < a.R &> \$of.log
		less \$of.log | grep -P "^Coef|Pearson's product-moment correlation" -A 9 -B 6 > \$of.stat
		rm tmp\$of.bed tmp2\$of.bed
		popd
	}

	bb_prepTFBStable() {
		local mtag=\$1
		local rtag=\$2
		local odir=\$3
		local nc=\$4
		local regions=self

		# prep TFBS table
		#motifSrc=JASPAR2014; motifSrcName=JASPAR2014; matrixtype="PWM"; all_versions=T;
		#motifSrc=~/PI_HOME/Data/casco/Homer/custom.motifs; motifSrcName=Homer; matrixtype="PWM"; ii=2; byrow=T;
		#motifSrc=~/PI_HOME/Data/casco/JASPAR/JASPAR_CORE_nonredundant_pfm_vertebrates.txt; motifSrcName=JASPARcoreVert; matrixtype="PFM"; ii=1; byrow=F;

		declare -A p=()
		if [[ \$mtag == JASPAR2014 ]]; then
			p=( [motifSrc]=JASPAR2014                                                               [matrixtype]=PWM [all_versions]=T  [ii]=NA [byrow]=NA)
		elif [[ \$mtag == JASPARCOREVERT ]]; then
			p=( [motifSrc]=~/PI_HOME/Data/casco/JASPAR/JASPAR_CORE_nonredundant_pfm_vertebrates.txt [matrixtype]=PFM [all_versions]=NA [ii]=1 [byrow]=F )
		elif [[ \$mtag == HOMER ]]; then
			p=( [motifSrc]=~/PI_HOME/Data/casco/Homer/custom.motifs                                 [matrixtype]=PWM [all_versions]=NA [ii]=2 [byrow]=T )
		else
			die ERROR: \$mtag not known
		fi

		mkdir -p \$odir; pushd \$odir
		local regf=\$(bb_setRegf \$rtag)
		Rscript -e "source('~/Projects/biter_biter_shared/bin/auxML.R')" \
			  -e "caller.getMotifCountsInRegion('\$regf', regions='\$regions',  src='\${p[motifSrc]}', matrixtype='\${p[matrixtype]}', all_versions=\${p[all_versions]}, ii=\${p[ii]}, byrow=\${p[byrow]}, ofleTag='tfbs', mc.cores=\$nc)"
		popd

		Rscript --version;
	}

	# start of first pair or single end read is counted with bedtools multicov
	bb_bam2start2bam() {
		local bamf=\$1
		local paired=\$2
		local obamf=\$3
		local nc=\$4

		if [[ \$paired == T ]]; then
			samtools view -b -@ \$nc -f 0x40 \$bamf | bedtools bamtobed -i stdin
		else 
			bedtools bamtobed -i \$bamf
		fi | awk 'BEGIN{OFS="\t";} { if (\$6=="+") { \$3=\$2+1; } else { \$2=\$3-1; } print \$0 }' | \
			bedtools bedtobam -i stdin -g ~/aux/genomes/mouse.mm10.genome | samtools sort -@ \$nc - > \$obamf
		samtools index \$obamf
	}

	bb_featureCounts() {
		local btag=\$1
		local rtag=\$2
		local odir=\$3
		local paired=\$4
		local nc=\$5

		export -f bb_bam2start2bam bb_setBAMf
		mkdir -p \$odir; pushd \$odir
		local bamf=(\$(bb_setBAMf \$btag))
		local regf=\$(bb_setRegf \$rtag)
		#parallel --delay 1 --plus --tmpdir . -j0 --tag bb_bam2start2bam {} \$paired {/} \$nc ::: \${bamf[@]} 
		#bamf=(\$(ls *sorted.bam))
		echo "\$(echo seqname start end name score strand \${bamf[@]##*/} | sed 's/.sorted.bam//g' | sed 's/ /\t/g')" > featureCounts
		bedtools multicov -bams \${bamf[@]} -bed \$regf >> featureCounts
		popd
	}

	export -f bb_prepTFBStable bb_setRegf bb_bam2start2bam bb_setBAMf bb_featureCounts
	
	mkdir -p tmp
	r=featureCount
	parallel --delay 1 --plus --tmpdir tmp -j0 --header : --tag --joblog joblog.\$r \
		"bb_featureCounts {BTAG} {RTAG} \$r/BTAG/{BTAG}/RTAG/{RTAG} $paired $nc" \
		::: BTAG ${BTAGs[@]} ::: RTAG ${RTAGs[@]}

#	r=TFBSinRegs
#	parallel --delay 1 --plus --tmpdir tmp -j0 --header : --tag --joblog joblog.\$r \
#		"bb_prepTFBStable {MTAG} {RTAG} \$r/MTAG/{MTAG}/RTAG/{RTAG} $nc" \
#		::: MTAG ${MTAGs[@]} ::: RTAG ${RTAGs[@]}

#	r=RegionPromoterAssoc
#	parallel --delay 1 --plus --tmpdir tmp -j0 --header : --tag --joblog joblog.\$r \
#		"bb_associateRegionWPromoters {RTAG} {MINDIRDIST} \$r/RTAG/{RTAG}/MINDIRDIST/{MINDIRDIST}" \
#		::: RTAG ${RTAGs[@]} ::: MINDIRDIST 1 0


	rm -rf tmp
END

sbatch -J $odir -p normal -t $t:0:0 -c $nc -e slurm_$n.out -o slurm_$n.out cmds.sh
popd

exit

