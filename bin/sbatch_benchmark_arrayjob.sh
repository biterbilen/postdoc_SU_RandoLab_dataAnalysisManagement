#!/bin/bash - 
#===============================================================================
#
#          FILE: sbatch_qc.sh
# 
#         USAGE: ./sbatch_qc.sh
# 
#   DESCRIPTION: 
# 
#       OPTIONS: ---
#  REQUIREMENTS: ---
#          BUGS: ---
#         NOTES: ---
#        AUTHOR: Dr. Biter Bilen (bb), biterbilen@yahoo.com
#  ORGANIZATION: Stanford, Palo Alto, USA
#       CREATED: 07/21/2015 09:51
#      REVISION:  ---
#===============================================================================

source ~/Projects/biter_biter_shared/bin/auxslurm.sh

#slurmMain QA --skippostproc=T --verbose=T --nc=2 --t=2 --datadir=~/Projects/biter_demorree_PAS/data/deMorree_2012/directRNAseq --pattern=.fasta.gz 
#slurmMain TR --skippostproc=T --verbose=T --nc=8 --t=2 --datadir=~/Projects/biter_biter_shared/results/2016-10-12 --indir=QA_deMorree_2012_directRNAseq --dryrun=T
slurmMain RM --skippostproc=T --verbose=T --nc=2 --t=2 --datadir=~/Projects/biter_biter_shared/results/2016-10-12 --indir=QA_deMorree_2012_directRNAseq --dryrun=F --skipfastqc=T
