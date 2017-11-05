#!/bin/bash - 
#===============================================================================
#
#          FILE: sbatchLike_w_segtools.sh
# 
#         USAGE: ./sbatchLike_w_segtools.sh 
# 
#   DESCRIPTION: 
# 
#       OPTIONS: ---
#  REQUIREMENTS: ---
#          BUGS: ---
#         NOTES: ---
#        AUTHOR: Dr. Biter Bilen (bb), biterbilen@yahoo.com
#  ORGANIZATION: Stanford, Palo Alto, USA
#       CREATED: 04/20/2016 10:21
#      REVISION:  ---
#===============================================================================

set -o nounset                              # Treat unset variables as an error
source ~/Projects/biter_biter_shared/bin/auxslurm.sh

SG() {
	local jobname
	local genomedir=${o[genomedir]:-X}
	# 0.
	local jobname=prepGenomedata; qcmd=$(cat <<- END
	genomefafile=~/PI_HOME/Data/casco/UCSC_tracks/${o[genome]:-mm10}/assembly/all/chr.fa
	genomedata-info --version
	genomedata-load -s \$genomefafile $genomedir
	cp -rf $genomedir $genomedir.back
	TNs=(DAyQSC DAoQSC DAyASC)
	TFs=(${o[datadir]}/${o[indir]:-PF}/pileup/WT_young_NA_SC_filterdup_pileup.bdg
		${o[datadir]}/${o[indir]:-PF}/pileup/WT_old_NA_SC_filterdup_pileup.bdg
		${o[datadir]}/${o[indir]:-PF}/pileup/WT_young_60hPI_SC_filterdup_pileup.bdg)
	genomedata-open-data $genomedir \${TNs[@]}
	for i in \$(seq 0 \$((\${#TNs[@]} - 1))); do
		tn=\${TNs[\$i]}
		tf=\${TFs[\$i]}
		less \$tf | genomedata-load-data $genomedir \$tn
	done
	genomedata-close-data $genomedir 
	genomedata-info tracknames_continuous $genomedir 
	END
	)
	if [[ ,${o[jobnames]:-X}, =~ ,$jobname, ]]; then qcmds+=( [$jobname]="$qcmd" ); fi

	# 2.
	# TODO add recover from a directory
	local traindir=${o[traindir]:-X}
	jobname=segwayTrain; qcmd=$(cat <<- END
	segway --version
	grep -w ${o[selchr]:-chr1} ~/aux/genomes/mouse.${o[genome]:-mm10}.genome | \
		awk '{ OFS="\t"; print \$1, 0, \$2 }' > includedInTraining.bed 
	SEGWAY_RAND_SEED=1498730685
	recover=${o[recover]:-F}
	if [ \$recover == F ]; then
		segway -c -v 3 \
			--track DAyASC --track DAyQSC --track DAoQSC --include-coords=includedInTraining.bed \
			--minibatch-fraction=0.01 --split-sequences=2000000 --resolution=200 \
			--distribution=asinh_norm --num-labels=8 --num-instances=${o[numInstances]:-1} \
			train $genomedir $traindir 
	else
		segway -r $traindir -v 3 \
			--track DAyASC --track DAyQSC --track DAoQSC --include-coords=includedInTraining.bed \
			--minibatch-fraction=0.01 --split-sequences=2000000 --resolution=200 \
			--distribution=asinh_norm --num-labels=8 --num-instances=${o[numInstances]:-1} \
			train $genomedir $traindir\2 
	fi
	END
	)
	if [ ${o[recover]:-F} == F ]; then
		if [[ ,${o[jobnames]:-X}, =~ ,$jobname, ]]; then dqcmds+=( [$jobname]="$qcmd" ); fi
	else
		if [[ ,${o[jobnames]:-X}, =~ ,$jobname, ]]; then dqcmds+=( [$jobname\2]="$qcmd" ); fi
	fi

	# 3.
	local identifydir=${o[identifydir]:-X}
	jobname=segwayIdentify; 
	for i in `seq 0 $((${o[numInstances]:-1}-1))`; do
		qcmd=$(cat <<- END
		mkdir -p $identifydir$i; pushd $identifydir$i
		segway --version
		SEGWAY_RAND_SEED=1498730685
		traindir=\$(ls -d ../$traindir* | tail -n 1)
		genomedir=../$genomedir
		ln -sf \$traindir . 
		ln -sf \$genomedir .
		traindir=\$(basename \$traindir)
		genomedir=$(basename \$genomedir)
		params=\$(ls \$traindir/params/params.$i.params.* | tail -n 1)
		segway -c -v 3 \
			--track DAyASC --track DAyQSC --track DAoQSC \
			--resolution=200 \
			--distribution=asinh_norm \
			-p \$params \
			identify \$genomedir \$traindir .
		popd
		END
		)
		if [[ ,${o[jobnames]:-X}, =~ ,$jobname, ]]; then dqcmds2+=( [$jobname$i]="$qcmd" ); fi
	done
}

#			-s \$traindir/segway.str \
#			-i \$traindir/params/input.$i.master \

slurmMain SG --skippostproc=F --datadir=~/Projects/biter_biter_stemCellFateRegulation/data/Liu_2015/ATACseq --indir=PF_RM_TR_ada \
	--jobnames=segwayTrain,segwayIdentify --genomedir=genomedata --traindir=segwayTrain --identifydir=segwayIdentify \
	--selchr=chr1 --numInstances=10 --nc=1 --t=6 --queue=bigmem

#slurmMain SG --skippostproc=F --datadir=~/Projects/biter_biter_stemCellFateRegulation/data/Liu_2015/ATACseq --indir=PF_RM_TR_ada \
#	--jobnames=prepGenomedata,segwayTrain,segwayIdentify --genomedir=genomedata --traindir=segwayTrain --identifydir=segwayIdentify \
#	--genomefafile=~/PI_HOME/Data/casco/UCSC_tracks/mm10/assembly/all/chr.fa --numInstances=2 --otag=numInstances2 --nc=16 --t=48
#	--jobnames=segwayTrain,segwayIdentify --genomedir=genomedata --traindir=segwayTrain --identifydir=segwayIdentify \
#	--genomefafile=~/PI_HOME/Data/casco/UCSC_tracks/mm10/assembly/all/chr.fa --numInstances=2 --nc=16 --t=48 --recover=T

# FIXME issues:
# The resolution feature is not implemented for posterior use in the Segway 1.4 release.
# they use 25 states (i.e. num-label)
#	experimental 
#	--num-sublabels=sublabels 
# --dry-run 
# --track dinucleotide
# -r traindir_res200.mm10
# TODO include and exclude coords files
#
#wget http://pmgenomics.ca/hoffmanlab/proj/segway/2011/test.genomedata
#genomedata-info tracknames_continuous test.genomedata
#genomedata-info contigs test.genomedata
#segway --num-labels=4 train test.genomedata traindir
#segway identify test.genomedata traindir identifydir


