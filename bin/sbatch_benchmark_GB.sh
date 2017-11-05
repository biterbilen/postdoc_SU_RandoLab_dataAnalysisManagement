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
#sbatch -J $odir -p normal -t $t:0:0 -c $nc -e slurm_$n.out -o slurm_$n.out <<- END
cat > cmds.sh <<- END
	#!/bin/bash -l
	bedGraphToBigWigTrack() {
		ct=$ct
		local params="macs2 callpeak -t ${samples[@]} \$ct --nomodel --shift \$2 --extsize 74 -g mm -n \$3 --trackline -B"
		declare -A cmds=()
		cmds[woCtlNoModelwSPMR]="--SPMR"
		cmds[woCtlNoModelwoSPMR]=""
		cmds[woCtlNoModelBroad]="--broad"
		\$params \${cmds[\$1]}

	}
	mkdir -p tmp

	export -f callpeak
	r=callpeak
	parallel --delay 1 --plus --tmpdir tmp -j0 --results \$r --header : \
		"mkdir -p \$r/TAG/{TAG}/SHIFT/{SHIFT}; callpeak {TAG} {SHIFT} \$r/TAG/{TAG}/SHIFT/{SHIFT}/callpeak" \
		::: TAG woCtlNoModelwSPMR woCtlNoModelwoSPMR woCtlNoModelBroad ::: SHIFT 0 -36

	r=bedGraphToBigWig
	parallel --delay 1 --plus --tmpdir tmp -j0 --results \$r --header : --xapply \
		"mkdir -p \$r/{DIR}; bedGraphToBigWig {BDG}.bdg ~/aux/genomes/mouse.mm10.genome.txt \$r/{DIR}/{BDG/}.bw" \
		::: BDG \$(find callpeak -name "*bdg" | sed 's/.bdg//g') ::: DIR \$(find callpeak -name "*bdg" | sed 's:[^/]*$::g' | sed 's:callpeak/::g' )

	rm -rf tmp
END

sbatch -J $odir -p normal -t $t:0:0 -c $nc -e slurm_$n.out -o slurm_$n.out cmds.sh
popd

exit


GB() {
	local analysisdirbase=$1
	#-----------------------
	# copy with links archive compress permissions 
	if [[ "${o[copy2server]:-F}" == T ]]; then
		rm -rf ${o[analysisdir]}/$analysisdirbase
		rsync -e 'ssh -p 22' -Lazp ${o[analysisdir]} \
			${o[webserverdir]:-"biter@cardinal.stanford.edu:/afs/ir/group/rando_lab/WWW/biter"}
		return
	fi

	local -a names=($(awk 'NR>1{print $NF}' metadata | sort | uniq))
	local bigDataUrl="${o[bigDataUrl]:-"http://www-dev.stanford.edu/group/rando_lab/biter"}/$analysisdirbase"
	local genomefafile=${o[genomefafile]:-"/home/biter/PI_HOME/Data/casco/HISAT_indexes/${o[genome]:-mm10}/all/chr.fa"}
	local selchr=${o[selchr]:-chr19}
	#-----------------------
	# hub.txt
	cat > hub.txt <<- END
		hub $analysisdirbase
		shortLabel $analysisdirbase
		longLabel $analysisdirbase
		genomesFile genomes.txt
		email ${o[email]:-"biter@stanford.edu"}
	END

	# genomes.txt
	cat > genomes.txt <<- END
		trackDb ${o[tracksfle]:=tracks.txt}
		genome ${o[genome]:-mm10}
	END

	# tracks.txt
	echo browser position ${o[position]:-"${selchr}:212000-250000"} > ${o[tracksfle]}

	local name
	# TODO bw currently treats paired-end data as single-end
	for name in ${names[@]}; do
		local qcmd=$(cat <<- END
			trackf() {
				local name=\$1
				local type=\$2
				local h="track name='$name' color=0,0,255"
				if [[ \$type == bigWig ]]; then
					echo type=bigWig bigDataUrl=$bigDataUrl/\$name.bw visibility=${o[visibility]:-squish} \$h
					bedtools genomecov -bg -ibam a.bam -g ~/aux/genomes/mouse.${o[genome]:-mm10}.genome > \$name.bg
					bedGraphToBigWig \$name.bg ~/aux/genomes/mouse.${o[genome]:-mm10}.genome \$name.bw
				elif [[ \$type == bam ]]; then
					echo type=bam bigDataUrl=$bigDataUrl/${selchr}_\$name.sorted.bam visibility=${o[visibility]:-hide} \$h
					samtools view -@ ${o[nc]} -b -o ${selchr}_\$name.sorted.bam ${o[datadir]}/${o[indir]:-RM}/\$name.sorted.bam $selchr
					samtools index ${selchr}_\$name.sorted.bam
					rm -rf ${selchr}_\$name.bam;
				fi >> ${o[tracksfle]}
				bedtools --version;
				samtools --version;
			}
			export -f trackf
			trackf $name ${o[
		END
		)
		qcmds+=( [$name]="$qcmd" )
	done

}

