#!/bin/bash - 
#===============================================================================
#
#          FILE: sbatch_pf.sh
# 
#         USAGE: ./sbatch_pf.sh 
# 
#   DESCRIPTION: 
# 
#       OPTIONS: ---
#  REQUIREMENTS: ---
#          BUGS: ---
#         NOTES: ---
#        AUTHOR: Dr. Biter Bilen (bb), biterbilen@yahoo.com
#  ORGANIZATION: Stanford, Palo Alto, USA
#       CREATED: 08/17/2015 12:03
#      REVISION:  ---
#===============================================================================

set -o nounset                              # Treat unset variables as an error
source ~/Projects/biter_biter_shared/bin/auxslurm.sh
po=${1:-F}

# 4 Fosl1 project
#slurmMain PF --datadir=$PI_SCRATCH/Projects/biter_demorree_Fosl1/data/Wold_2011/ChIPseq --nc=2 --t=2 --parseOutput=$po

# 4 Ling 2017
slurmMain PF --datadir=$PI_SCRATCH/Projects/biter_lingliu_DNAdamage/data/Liu_2016/ChIPseqV2 --nc=2 --t=2 --parseOutput=$po #--broad=True

exit

cs=${1:-F}
pr=${2:-F}
pd=${3:-F}
fd=${4:-F}
pu=${5:-F}
bd=${6:-F}
gb=${7:-F}
ar=${8:-F}

cp=${9:-F}


# 4 Biter
#slurmMain PF --datadir=~/Projects/biter_biter_stemCellFateRegulation/data/Liu_2015/ATACseq --cutoff=10 --nc=16 --stag=WT --indir=RM_TR_ada --calcsim=$cs --poolreps=$pr --predictd=$pd --filterdup=$fd --pileup=$pu \
#	--bdgdiff=$bd --t1=WT_old_NA_SC      --t2=WT_young_NA_SC     --c1=NULL --c2=NULL                     --GB=$gb --annotateRegion=$ar     --callpeak=$cp
#	--bdgdiff=$bd --t1=WT_old_NA_SC      --t2=WT_young_60hPI_SC  --c1=WT_young_NA_SC --c2=WT_young_NA_SC --GB=$gb --annotateRegion=$ar     --callpeak=$cp

# 4 Ling
#slurmMain PF --datadir=~/Projects/biter_lingliu_DNAdamage/data/Liu_2016/ChIPseq --nc=16 --stag=Hr --calcsim=$cs --poolreps=$pr --predictd=$pd --filterdup=$fd --pileup=$pu \
slurmMain PF --datadir=~/Projects/biter_lingliu_DNAdamage/data/Liu_2016/ChIPseqV2 --nc=16 --stag=Hr --calcsim=$cs --poolreps=$pr --predictd=$pd --filterdup=$fd --pileup=$pu \
	--bdgdiff=$bd --t1=HrWT_young_TMX_SC_NA_H3K9me2 --t2=HrWT_young_TMX_SC_NA_H3K9me3 --c1=HrWT_young_TMX_SC_NA_NA --c2=HrWT_young_TMX_SC_NA_NA      --GB=$gb --annotateRegion=$ar --callpeak=$cp
#	--bdgdiff=$bd --t1=HrCKO_young_TMX_SC_NA_H3K9me3 --t2=HrWT_young_TMX_SC_NA_H3K9me3 --c1=HrCKO_young_TMX_SC_NA_NA --c2=HrWT_young_TMX_SC_NA_NA    --GB=$gb --annotateRegion=$ar     --callpeak=$cp
#	--bdgdiff=$bd --t1=HrCKO_young_TMX_SC_NA_H3K9me2 --t2=HrWT_young_TMX_SC_NA_H3K9me2 --c1=HrCKO_young_TMX_SC_NA_NA --c2=HrWT_young_TMX_SC_NA_NA    --GB=$gb --annotateRegion=$ar --callpeak=$cp
#	--bdgdiff=$bd --t1=HrCKO_young_TMX_SC_NA_H3K9me2 --t2=HrCKO_young_TMX_SC_NA_H3K9me3 --c1=HrCKO_young_TMX_SC_NA_NA --c2=HrCKO_young_TMX_SC_NA_NA  --GB=$gb --annotateRegion=$ar --callpeak=$cp

#slurmMain PF --datadir=~/Projects/biter_lingliu_DNAdamage/data/Liu_2016/ChIPseq2 --nc=16 --stag=WT --calcsim=$cs --poolreps=$pr --predictd=$pd --filterdup=$fd --pileup=$pu \
#	--bdgdiff=$bd --t1=WT_old_NA_SC_H3K9me2 --t2=WT_old_NA_SC_H3K9me3 --c1=WT_old_NA_SC_NA --c2=WT_old_NA_SC_NA    --GB=$gb --annotateRegion=$ar     --callpeak=$cp
#	--bdgdiff=$bd --t1=WT_old_NA_SC_H3K9me3 --t2=WT_young_NA_SC_H3K9me3 --c1=WT_old_NA_SC_NA --c2=WT_young_NA_SC_NA    --GB=$gb --annotateRegion=$ar     --callpeak=$cp
#	--bdgdiff=$bd --t1=WT_old_NA_SC_H3K9me2 --t2=WT_young_NA_SC_H3K9me2 --c1=WT_old_NA_SC_NA --c2=WT_young_NA_SC_NA    --GB=$gb --annotateRegion=$ar     --callpeak=$cp
#	--bdgdiff=$bd --t1=WT_young_NA_SC_H3K9me2 --t2=WT_young_NA_SC_H3K9me3 --c1=WT_young_NA_SC_NA --c2=WT_young_NA_SC_NA    --GB=$gb --annotateRegion=$ar     --callpeak=$cp

#slurmMain PF --datadir=~/Projects/biter_lingliu_DNAdamage/data/Liu_2015/ChIPseq --stag=WT --calcsim=$cs --poolreps=$pr --predictd=$pd --filterdup=$fd --pileup=$pu \
#	--bdgdiff=$bd --t1=WT_young_NA_SC_NA_H3K9me2 --t2=WT_young_NA_SC_NA_H3K9me3 --c1=WT_young_NA_SC_NA_NA --c2=WT_young_NA_SC_NA_NA                  --GB=$gb --annotateRegion=$ar     --callpeak=$cp

#slurmMain PF --datadir=~/Projects/biter_lingliu_DNAdamage/data/Liu_2015/ChIPseq2 --indir=RM_TR_ada --nc=16 --stag=WT --calcsim=$cs --poolreps=$pr --predictd=$pd --filterdup=$fd --pileup=$pu \
#	--bdgdiff=$bd --t1=HrCKO_young_TMX_SC_NA_H3K9me2 --t2=HrWT_young_TMX_SC_NA_H3K9me2 --c1=HrCKO_young_TMX_SC_NA_NA --c2=HrWT_young_TMX_SC_NA_NA      --GB=$gb --annotateRegion=$ar   --callpeak=$cp
#	--bdgdiff=$bd --t1=HrCKO_young_TMX_SC_NA_H3K9me3 --t2=HrWT_young_TMX_SC_NA_H3K9me3 --c1=HrCKO_young_TMX_SC_NA_NA --c2=HrWT_young_TMX_SC_NA_NA      --GB=$gb --annotateRegion=$ar --callpeak=$cp
#	--bdgdiff=$bd --t1=HrCKO_young_TMX_SC_NA_H3K9me2 --t2=HrCKO_young_TMX_SC_NA_H3K9me3 --c1=HrCKO_young_TMX_SC_NA_NA --c2=HrCKO_young_TMX_SC_NA_NA    --GB=$gb --annotateRegion=$ar --callpeak=$cp 
#	--bdgdiff=$bd --t1=HrWT_young_TMX_SC_NA_H3K9me2 --t2=HrWT_young_TMX_SC_NA_H3K9me3 --c1=HrWT_young_TMX_SC_NA_NA --c2=HrWT_young_TMX_SC_NA_NA      --GB=$gb --annotateRegion=$ar --callpeak=$cp
	 

#slurmMain PF --datadir=~/Projects/biter_lingliu_DNAdamage/data/Liu_2015/ATACseq --nc=16 --indir=RM_TR_ada --stag=WT --otag=yAyQoQ 

#slurmMain PF --datadir=~/Projects/biter_lingliu_DNAdamage/data/Liu_2016/ChIPseq --nc=16 --stag=HrCKO_young_TMX_SC_NA_H3K9me2 --ctag=HrCKO_young_TMX_SC_NA_NA --otag=HrCKOK9me2
#slurmMain PF --datadir=~/Projects/biter_lingliu_DNAdamage/data/Liu_2016/ChIPseq --nc=16 --stag=HrCKO_young_TMX_SC_NA_H3K9me2 --ctag=HrCKO_young_TMX_SC_NA_NA --otag=HrCKOK9me2
#slurmMain PF --datadir=~/Projects/biter_lingliu_DNAdamage/data/Liu_2016/ChIPseq --nc=16 --stag=HrCKO_young_TMX_SC_NA_H3K9me3 --ctag=HrCKO_young_TMX_SC_NA_NA --otag=HrCKOK9me3
#slurmMain PF --datadir=~/Projects/biter_lingliu_DNAdamage/data/Liu_2016/ChIPseq --nc=16 --stag=HrWT_young_TMX_SC_NA_H3K9me2 --ctag=HrWT_young_TMX_SC_NA_NA --otag=HrWTK9me2
#slurmMain PF --datadir=~/Projects/biter_lingliu_DNAdamage/data/Liu_2016/ChIPseq --nc=16 --stag=HrWT_young_TMX_SC_NA_H3K9me3 --ctag=HrWT_young_TMX_SC_NA_NA --otag=HrWTK9me3

exit

slurmMain PF --datadir=~/Projects/biter_lingliu_DNAdamage/data/Liu_2016/ChIPseq --stag=HrCKO_young_TMX_SC_NA_H3K9me2_1 --ctag=HrCKO_young_TMX_SC_NA_NA --otag=HrCKOK9me21
slurmMain PF --datadir=~/Projects/biter_lingliu_DNAdamage/data/Liu_2016/ChIPseq --stag=HrCKO_young_TMX_SC_NA_H3K9me3_1 --ctag=HrCKO_young_TMX_SC_NA_NA --otag=HrCKOK9me31
slurmMain PF --datadir=~/Projects/biter_lingliu_DNAdamage/data/Liu_2016/ChIPseq --stag=HrWT_young_TMX_SC_NA_H3K9me2_1 --ctag=HrWT_young_TMX_SC_NA_NA --otag=HrWTK9me21
slurmMain PF --datadir=~/Projects/biter_lingliu_DNAdamage/data/Liu_2016/ChIPseq --stag=HrWT_young_TMX_SC_NA_H3K9me3_1 --ctag=HrWT_young_TMX_SC_NA_NA --otag=HrWTK9me31
slurmMain PF --datadir=~/Projects/biter_lingliu_DNAdamage/data/Liu_2016/ChIPseq --stag=HrCKO_young_TMX_SC_NA_H3K9me2_2 --ctag=HrCKO_young_TMX_SC_NA_NA --otag=HrCKOK9me22
slurmMain PF --datadir=~/Projects/biter_lingliu_DNAdamage/data/Liu_2016/ChIPseq --stag=HrCKO_young_TMX_SC_NA_H3K9me3_2 --ctag=HrCKO_young_TMX_SC_NA_NA --otag=HrCKOK9me32
slurmMain PF --datadir=~/Projects/biter_lingliu_DNAdamage/data/Liu_2016/ChIPseq --stag=HrWT_young_TMX_SC_NA_H3K9me2_2 --ctag=HrWT_young_TMX_SC_NA_NA --otag=HrWTK9me22
slurmMain PF --datadir=~/Projects/biter_lingliu_DNAdamage/data/Liu_2016/ChIPseq --stag=HrWT_young_TMX_SC_NA_H3K9me3_2 --ctag=HrWT_young_TMX_SC_NA_NA --otag=HrWTK9me32

