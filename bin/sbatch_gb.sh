#!/bin/bash - 
#===============================================================================
#
#          FILE: sbatch_gb.sh
# 
#         USAGE: ./sbatch_gb.sh 
# 
#   DESCRIPTION: 
# 
#       OPTIONS: ---
#  REQUIREMENTS: ---
#          BUGS: ---
#         NOTES: ---
#        AUTHOR: Dr. Biter Bilen (bb), biterbilen@yahoo.com
#  ORGANIZATION: Stanford, Palo Alto, USA
#       CREATED: 08/09/2015 16:56
#      REVISION:  ---
#===============================================================================

source ~/Projects/biter_biter_shared/bin/auxslurm.sh

copy2server=${1:-F} 
# 4 Antoine
# decide which RM to chose
#slurmMain GB --datadir=~/Projects/biter_demorree_PAS/data/deMorree_2012/directRNAseq --copy2server=$copy2server --selchr=chr10 --viewLimitsmin=1 --indir=RM_TR_trim5termT_L04minscoreNSC --nc=2 --t=1 

# 4 Mike
#slurmMain GB --datadir=~/Projects/biter_wosczyna_miRNAsInMuscleCellFate/data/Wosczyna_2015/RNAseq --copy2server=$copy2server --nc=2 --t=1
#slurmMain GB --datadir=~/Projects/biter_wosczyna_miRNAsInMuscleCellFate/data/Wosczyna_2015/PullDown --copy2server=$copy2server --nc=2 --t=1

# 4 Jamie
#slurmMain GB --datadir=~/Projects/biter_jbrett_DNAmethylation/data/Brett_2015/RNAseq --copy2server=$copy2server --nc=2 --t=1

# 4 Ling - Hairless paper submission
slurmMain GB --datadir=~/Projects/biter_lingliu_DNAdamage/data/Liu_2016/ATACseq --copy2server=$copy2server --nc=2 --t=1 --selchr=chr10 --indir=RM_TR

# 4 Ling
#slurmMain GB --datadir=~/Projects/biter_lingliu_DNAdamage/data/Liu_2016/ChIPseq --copy2server=$copy2server --nc=2 --t=1 --selchr=chr10
#slurmMain GB --datadir=~/Projects/biter_lingliu_DNAdamage/data/Liu_2015/ATACseq --copy2server=$copy2server --nc=2 --t=1 --indir=RM_TR_ada --selchr=chr10

#slurmMain GB --datadir=~/Projects/biter_lingliu_DNAdamage/data/Liu_2014/MNaseseq --copy2server=$copy2server --nc=2 --t=1 --indir=RM_TR_ada
#slurmMain GB --datadir=~/Projects/biter_lingliu_DNAdamage/data/Liu_2015/ChIPseq --copy2server=$copy2server --nc=2 --t=1 --indir=RM_TR_ada
#slurmMain GB --datadir=~/Projects/biter_lingliu_DNAdamage/data/Liu_2011/ChIPseq --copy2server=$copy2server --nc=2 --t=1
#slurmMain GB --datadir=~/Projects/biter_lingliu_DNAdamage/data/Wenemoser_2015/ChIPseq --copy2server=$copy2server --nc=2 --t=1 --indir=RM_TR_ada

