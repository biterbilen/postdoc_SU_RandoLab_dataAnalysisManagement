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

t=6
nc=16

# single-end 
# chromatin data --extsize could be set to 73 for histone marks; --shift could be set to -HALFofEXTSIZE 
indir=~/Projects/biter_lingliu_DNAdamage/data/Liu_2016/ChIPseq/RM/
stag=HrCKO_young_TMX_SC_NA_H3K9me2;ctag=HrCKO_young_TMX_SC_NA_NA;n=HrCKOK9me2

indir=~/Projects/biter_lingliu_DNAdamage/data/Wenemoser_2015/ChIPseq/RM_TR_ada/
stag=WT_11mo_NA_SC_NA_H3K27me3_1;ctag=WT_11mo_NA_SC_NA_NA_1;n=Wene11moK27

# ATAC-seq
indir=~/Projects/biter_lingliu_DNAdamage/data/Liu_2015/ATACseq/RM_TR_ada/
stag=WT;ctag=NULL;n=Acc4yAyQoQ

samples=($(ls $indir/*$stag*.sorted.bam))
controls=($(ls $indir/*$ctag*.sorted.bam 2> /dev/null))
ct=""
if [ ${#controls[@]} -gt 0 ]; then
	ct="-c ${controls[@]}"; echo Control is: $ct
	ln -sf ${controls[@]} .
fi

tool=macs2
odir=PF_$tool\_$n

mkdir -p $odir; pushd $odir;
ln -sf ${samples[@]} .

#echo $cmd
#callpeak,bdgpeakcall,bdgbroadcall,bdgcmp,bdgopt,cmbreps,bdgdiff,filterdup,predictd,pileup,randsample,refinepeak
# work on --rpl below
#cat > cmds.sh <<- END
sbatch -J $odir -p normal -t $t:0:0 -c $nc -e slurm_$n.out -o slurm_$n.out <<- END
	#!/bin/bash -l
	callpeak() {
		ct=$ct
		local params="macs2 callpeak -t ${samples[@]} \$ct --nomodel --shift \$2 --extsize 73 -g mm -n \$3 --trackline -B"
		declare -A cmds=()
		cmds[woCtlNoModelwSPMR]="--SPMR"
		cmds[woCtlNoModelwoSPMR]=""
		cmds[woCtlNoModelBroad]="--broad"
		\$params \${cmds[\$1]}

	}
	mkdir -p tmp
	r=filterdup
	parallel --delay 1 --plus --tmpdir tmp -j0 --results \$r --header : \
		"mkdir -p \$r/BAM; macs2 filterdup -g mm -i {BAM}.sorted.bam -o \$r/BAM/{BAM}/filterdup.bed" \
		::: BAM \$(ls *sorted.bam | sed 's/.sorted.bam//g')

	r=predictd
	parallel --delay 1 --plus --tmpdir tmp -j0 --results \$r --header : \
		"mkdir -p \$r/BAM; macs2 predictd -g mm -i {BAM}.sorted.bam --rfile \$r/BAM/{BAM}/predictd.R; Rscript \$r/BAM/{BAM}/predictd.R" \
		::: BAM \$(ls *sorted.bam | sed 's/.sorted.bam//g')

	export -f callpeak
	r=callpeak
	parallel --delay 1 --plus --tmpdir tmp -j0 --results $\r --header : \
		"mkdir -p \$r/TAG/{TAG}/SHIFT/{SHIFT}; eval callpeak {TAG} {SHIFT} \$r/TAG/{TAG}/SHIFT/{SHIFT}/callpeak" \
		::: TAG woCtlNoModelwSPMR woCtlNoModelwoSPMR woCtlNoModelBroad ::: SHIFT 0 -38

	r=bedGraphToBigWig
	parallel --delay 1 --plus --tmpdir tmp -j0 --results \$r --header : --xapply \
		"mkdir -p \$r/{DIR}; bedGraphToBigWig {BDG}.bdg ~/aux/genomes/mouse.mm10.genome.txt \$r/{DIR}/{BDG/}.bw" \
		::: BDG \$(find callpeak -name "*bdg" | sed 's/.bdg//g') ::: DIR \$(find callpeak -name "*bdg" | sed 's:[^/]*$::g' | sed 's:callpeak/::g' )
	rm -rf tmp
END
popd

exit



r=bedGraphToBigWig
parallel --delay 1 --plus --tmpdir tmp -j0 --results \$r --header : --xapply \
	"echo mkdir -p $r/{DIR}; echo bedGraphToBigWig {BDG}.bdg ~/aux/genomes/mouse.mm10.genome.txt $r/{DIR}/{BDG/}.bw" \
	::: BDG $(find callpeak -name "*bdg" | sed 's/.bdg//g') ::: DIR $(find callpeak -name "*bdg" | sed 's:[^/]*$::g' | sed 's:callpeak/::g' )



