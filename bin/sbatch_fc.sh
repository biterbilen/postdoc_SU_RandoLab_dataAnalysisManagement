#!/bin/bash - 
#===============================================================================
#
#          FILE: sbatch_fc.sh
# 
#         USAGE: ./sbatch_fc.sh 
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
source $PI_SCRATCH/Projects/biter_biter_shared/bin/auxslurm.sh
po=${1:-F}

# For Fosl1 project
#slurmMain FC --datadir=$PI_SCRATCH/Projects/biter_demorree_Fosl1/data/Wold_2011/ChIPseq --bamsuf=.ddnsrt.bam --binsize=1000 --plotfrac=0.1 --nc=8 --t=2 --parseOutput=$po

# 4 Ling 2017
#slurmMain FC --datadir=$PI_SCRATCH/Projects/biter_lingliu_DNAdamage/data/Liu_2016/ChIPseqV2 --binsize=10000 --nc=8 --t=1 --parseOutput=$po
#slurmMain FC --datadir=$PI_SCRATCH/Projects/biter_lingliu_DNAdamage/data/Liu_2017/ATACseq --binsize=10000 --maxins=2000 --nc=8 --t=2 --samplerm=T --parseOutput=$po
#slurmMain FC --datadir=$PI_SCRATCH/Projects/biter_lingliu_DNAdamage/data/Liu_2017/ATACseq --binsize=10000 --maxins=150 --nc=8 --t=2 --otag=oc --samplerm=T --parseOutput=$po
#slurmMain FC --datadir=$PI_SCRATCH/Projects/biter_lingliu_DNAdamage/data/Liu_2016/ATACseq --binsize=10000 --maxins=2000 --nc=8 --t=2 --parseOutput=$po
#slurmMain FC --datadir=$PI_SCRATCH/Projects/biter_lingliu_DNAdamage/data/Liu_2016/ATACseq --binsize=10000 --maxins=150 --nc=8 --t=2 --otag=oc --parseOutput=$po
#slurmMain FC --datadir=$PI_SCRATCH/Projects/biter_lingliu_DNAdamage/data/Liu_2016/RNAseq --t=1 --parseOutput=$po 
slurmMain FC --datadir=$PI_SCRATCH/Projects/biter_lingliu_DNAdamage/data/Liu_2015/RNAseq --t=1 --parseOutput=$po 

# -----------------
# 4 Mike 2017
#slurmMain FC --datadir=$PI_SCRATCH/Projects/biter_wosczyna_miRNAsInMuscleCellFate/data/Wosczyna_2015/RNAseq --nc=4 --t=1 --parseOutput=$po 
#slurmMain FC --datadir=$PI_SCRATCH/Projects/biter_wosczyna_miRNAsInMuscleCellFate/data/Wosczyna_2015/PullDown --genome=mm10_ERCC92 --nc=8 --t=1 --parseOutput=$po

# 4 Jamie Luiz Daniel 2017
#slurmMain FC --datadir=$PI_SCRATCH/Projects/biter_jbrett_DNAmethylation/data/Brett_2015/RNAseq --genome=mm10_ERCC92 --nc=16 --parseOutput=$po

# 4 Biter
#slurmMain FC --datadir=$PI_SCRATCH/Projects/biter_biter_noiseInAging/data/Guess_2015/RNAseq
#slurmMain FC --datadir=$PI_SCRATCH/Projects/biter_biter_stemCellFateRegulation/data/Ikeda_2014/RNAseq --nc=16 --t=12 --parseOutput=$po
#slurmMain FC --datadir=$PI_SCRATCH/Projects/biter_biter_stemCellFateRegulation/data/Brett_2015/RNAseq --genome=mm10_ERCC92 --nc=16 --t=12 --parseOutput=$po
#slurmMain FC --datadir=$PI_SCRATCH/Projects/biter_biter_stemCellFateRegulation/data/Liu_2014/RNAseq --nc=16 --parseOutput=$po

# 4 Jamie
#slurmMain FC --datadir=$PI_SCRATCH/Projects/biter_jbrett_DNAmethylation/data/Brett_2015/RNAseq --genome=mm10_ERCC92 --nc=16 --t=12
#slurmMain FC --datadir=$PI_SCRATCH/Projects/biter_jbrett_DNAmethylation/data/Ikeda_2014/RNAseq --genome=mm10_ERCC92 --nc=16 --t=12

# 4 Ling - Hairless paper submission
#slurmMain FC --datadir=$PI_SCRATCH/Projects/biter_lingliu_DNAdamage/data/Liu_2016/RNAseq --genome=mm10 --nc=16 --parseOutput=$po

# 4 Ling
#slurmMain FC --datadir=$PI_SCRATCH/Projects/biter_lingliu_DNAdamage/data/Liu_2014/RNAseq --genome=mm10 --nc=16
#slurmMain FC --datadir=$PI_SCRATCH/Projects/biter_lingliu_DNAdamage/data/Cheung_2013/RNAseq --genome=mm10 --nc=16
#slurmMain FC --datadir=$PI_SCRATCH/Projects/biter_lingliu_DNAdamage/data/Liu_2015/RNAseq --genome=mm10 --nc=16 --parseOutput=$po
#slurmMain FC --datadir=$PI_SCRATCH/Projects/biter_lingliu_DNAdamage/data/Liu_2016/RNAseq --genome=mm10 --nc=16 --parseOutput=$po

# 4 Cindy
#slurmMain FC --datadir=$PI_SCRATCH/Projects/biter_cvanvelt_SCTrxTrl/data/vanVelthoven_2015/PullDownTrx --genome=mm10_ERCC92 --nc=16 --t=12
#slurmMain FC --datadir=$PI_SCRATCH/Projects/biter_cvanvelt_SCTrxTrl/data/vanVelthoven_2015/PullDownTrl --genome=mm10_ERCC92 --nc=16 --t=12
#slurmMain FC --datadir=$PI_SCRATCH/Projects/biter_cvanvelt_SCTrxTrl/data/vanVelthoven_2015/PullDownTrx2 --genome=mm10_ERCC92 --nc=16 --t=12 --indir=RM_TR_ada --parseOutput=$po
#slurmMain FC --datadir=$PI_SCRATCH/Projects/biter_cvanvelt_SCTrxTrl/data/vanVelthoven_2016/PullDownTrx --genome=mm10_ERCC92 --nc=16 --t=12 --indir=RM_TR_ada --parseOutput=$po

# 4 Mike
#slurmMain FC --datadir=$PI_SCRATCH/Projects/biter_wosczyna_miRNAsInMuscleCellFate/data/Wosczyna_2015/RNAseq 
#slurmMain FC --datadir=$PI_SCRATCH/Projects/biter_wosczyna_miRNAsInMuscleCellFate/data/Wosczyna_2015/PullDown --genome=mm10_ERCC92

