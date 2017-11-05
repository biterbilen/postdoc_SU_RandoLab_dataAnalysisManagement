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

#Wenemoser_2015
indir=~/Projects/biter_lingliu_DNAdamage/data/Wenemoser_2015/ChIPseq/RM_TR_ada/
stag=WT_11mo_NA_SC_NA_H3K27me3_1;ctag=WT_11mo_NA_SC_NA_NA_1;n=Wene11moK27

# ATAC-seq
indir=~/Projects/biter_lingliu_DNAdamage/data/Liu_2015/ATACseq/RM_TR_ada/
stag=WT; ctag=NULL;n=Acc4yAyQoQ; TAGs=(wSPMR wSPMRBroad);  SHIFTs=(0)
#stag=WT;ctag=NULL;n=Acc4yAyQoQ; TAGS=woSPMR
#stag=WT;ctag=NULL;n=Acc4yAyQoQ; TAGS=woSPMRBroad
#stag=WT;ctag=NULL;n=Acc4yAyQoQ; TAGS=wSPMRBroad
#stag=WT_old_NA_SC_1;ctag=NULL;n=oQ1; TAGS=wSPMRBroad

#odir=PF_${TAGS}_$n
odir=PF_SHIFT${SHIFTs}_$n

mkdir -p $odir; pushd $odir;

#echo $cmd
#callpeak,bdgpeakcall,bdgbroadcall,bdgcmp,bdgopt,cmbreps,bdgdiff,filterdup,predictd,pileup,randsample,refinepeak
#sbatch -J $odir -p normal -t $t:0:0 -c $nc -e slurm_$n.out -o slurm_$n.out <<- END
# Too much load on one node; could be sent to separate nodes!!!
cat > cmds.sh <<- END
	#!/bin/bash -l
	bb_scaleTruncBedScore() {
		# Scales scores to 1-1000 proportionally; like in Thales' theorem
		local bedp=\$1
		local sbedp=\$2
		local min=\$(less \$bedp | grep ^chr | cut -f 5 | sort -g | head -n 1)
		local max=\$(less \$bedp | grep ^chr | cut -f 5 | sort -gr| head -n 1)
		less \$bedp | cut -f 1-6 | awk -v min=\$min -v max=\$max '\$1 ~ /chr/ { OFS="\t"; \$5=int((1000-1)*(\$5-min)/(max-min)+1); print \$0;}' > \$sbedp
	}

	bb_callpeak() {
		local pkey=\$1
		local shift=\$2
		local odir=\$3
		local indir=\$4
		local stag=\$5
		local ctag=\$6

		local samples=(\$(ls \$indir/*\$stag*.sorted.bam))
		local controls=(\$(ls \$indir/*\$ctag*.sorted.bam 2> /dev/null))
		local ct=""
		if [ \${#controls[@]} -gt 0 ]; then
			ct="-c \${controls[@]}"; echo Control is: \$ct
			ln -sf \${controls[@]} .
		fi
		ln -sf \${samples[@]} .


		local params="--nomodel --shift \$shift --extsize 74 -g mm -n callpeak --trackline -B -t \${samples[@]} \$ct"
		declare -A params2=()
		params2[wSPMR]="--SPMR"
		params2[woSPMR]=""
		params2[woSPMRBroad]="--broad"
		params2[wSPMRBroad]="--SPMR --broad"

		mkdir -p \$odir; pushd \$odir
		macs2 callpeak \$params \${params2[\$pkey]}
		parallel --tmpdir . 'bedGraphToBigWig {} ~/aux/genomes/mouse.mm10.genome.txt {.}.bw' ::: *bdg 
		export -f bb_scaleTruncBedScore
		parallel --tmpdir . 'bb_scaleTruncBedScore {} tmp{}; bedToBigBed tmp{} ~/aux/genomes/mouse.mm10.genome.txt {}.bb; rm tmp{}' ::: *bed
		parallel --tmpdir . 'bb_scaleTruncBedScore {} tmp{}; bedToBigBed tmp{} ~/aux/genomes/mouse.mm10.genome.txt {}.bb; rm tmp{}' ::: *narrowPeak
		parallel --tmpdir . 'bb_scaleTruncBedScore {} tmp{}; bedToBigBed tmp{} ~/aux/genomes/mouse.mm10.genome.txt {}.bb; rm tmp{}' ::: *broadPeak
		parallel --tmpdir . 'bb_scaleTruncBedScore {} tmp{}; bedToBigBed tmp{} ~/aux/genomes/mouse.mm10.genome.txt {}.bb; rm tmp{}' ::: *gappedPeak
		popd

	}

	bb_track() {
		local d=\$1
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
		for s in \$(find . -name *bb 2> /dev/null); do
			f=\$(echo \$s | sed 's:./callpeak/::' | sed 's:/:_:g')
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
		for s in \$(find . -name *treat*bw 2> /dev/null); do
			f=\$(echo \$s | sed 's:./callpeak/::' | sed 's:/:_:g')
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

	mkdir -p tmp

	export -f bb_callpeak bb_scaleTruncBedScore
	r=callpeak
	parallel --delay 1 --plus --tmpdir tmp -j0 --header : --tag --joblog joblog.\$r \
		"bb_callpeak {TAG} {SHIFT} \$r/TAG/{TAG}/SHIFT/{SHIFT} $indir $stag $ctag" \
		::: TAG ${TAGs[@]} ::: SHIFT ${SHIFTs[@]}

	bb_track GB_$odir

	rm -rf tmp
END

sbatch -J $odir -p normal -t $t:0:0 -c $nc -e slurm_$n.out -o slurm_$n.out cmds.sh
popd

exit

