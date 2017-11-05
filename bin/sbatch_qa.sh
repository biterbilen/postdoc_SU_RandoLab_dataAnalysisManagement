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

source $PI_SCRATCH/Projects/biter_biter_shared/bin/auxslurm.sh
# 4 Fosl1 project
#slurmMain QA --datadir=$PI_SCRATCH/Projects/biter_demorree_Fosl1/data/Wold_2011/ChIPseq

# 4 Mike 2017/04
#slurmMain QA --datadir=$PI_SCRATCH/Projects/biter_wosczyna_miRNAsInMuscleCellFate/data/Wosczyna_2015/RNAseq
#slurmMain QA --datadir=$PI_SCRATCH/Projects/biter_wosczyna_miRNAsInMuscleCellFate/data/Wosczyna_2015/PullDown

# 4 Ling 2017
#slurmMain QA --datadir=$PI_SCRATCH/Projects/biter_lingliu_DNAdamage/data/Liu_2016/ATACseq 
#slurmMain QA --datadir=$PI_SCRATCH/Projects/biter_lingliu_DNAdamage/data/Liu_2017/ATACseq
#slurmMain QA --datadir=$PI_SCRATCH/Projects/biter_lingliu_DNAdamage/data/Liu_2016/ChIPseqV2 
#slurmMain QA --datadir=$PI_SCRATCH/Projects/biter_lingliu_DNAdamage/data/Liu_2016/RNAseq --nc=16
slurmMain QA --datadir=$PI_SCRATCH/Projects/biter_lingliu_DNAdamage/data/Liu_2015/RNAseq --nc=16

# 4 Cindy 2017
#slurmMain QA --datadir=$PI_SCRATCH/Projects/biter_cvanvelt_SCTrxTrl/data/vanVelthoven_2016/PullDownTrx2

# 4 Jamie Luiz Daniel 2017
#slurmMain QA --datadir=$PI_SCRATCH/Projects/biter_jbrett_DNAmethylation/data/Brett_2015/RNAseq --skipfastqc=T
#slurmMain QA --datadir=$PI_SCRATCH/Projects/biter_jbrett_DNAmethylation/data/Ikeda_2014/RNAseq --skipfastqc=T
#slurmMain QA --datadir=$PI_SCRATCH/Projects/biter_lingliu_DNAdamage/data/Liu_2014/RNAseq --skipfastqc=T

# 4 Jamie
#slurmMain QA --datadir=$PI_SCRATCH/Projects/biter_jbrett_DNAmethylation/data/Brett_2014/WGBseq --nc=8 --t=8 --skipfastqc=T
##slurmMain QA --datadir=$PI_SCRATCH/Projects/biter_jbrett_DNAmethylation/data/Brett_2014/WGBseq --nc=16 --t=24
#slurmMain QA --datadir=$PI_SCRATCH/Projects/biter_jbrett_DNAmethylation/data/Brett_2015/RNAseq 
#slurmMain QA --datadir=$PI_SCRATCH/Projects/biter_jbrett_DNAmethylation/data/Ikeda_2014/RNAseq 

# 4 Jay
#slurmMain QA --datadir=$PI_SCRATCH/Projects/biter_jsalvi_spuriousReplication/data/Salvi_2016/ChIPseq
#slurmMain QA --datadir=$PI_SCRATCH/Projects/biter_jsalvi_spuriousReplication/data/Dellino_2013_GenomeResearch/ChIPseq

# 4 Antoine
#slurmMain QA --datadir=$PI_SCRATCH/Projects/biter_demorree_PAS/data/deMorree_2012/directRNAseq --pattern=.fasta.gz --nc=2 --t=2
#slurmMain QA --datadir=$PI_SCRATCH/Projects/biter_demorree_PAS/data/Soleimani_2012_DevCell/ChIPseq

# 4 Biter
#slurmMain QA --datadir=$PI_SCRATCH/Projects/biter_biter_noiseInAging/data/Guess_2015/RNAseq
#slurmMain QA --datadir=$PI_SCRATCH/Projects/biter_biter_EWSR1/data/Bilen_2015/ChIPseq --limitsfile=$PI_SCRATCH/Projects/biter_biter_EWSR1/results/2015-12-21/limits_Bilen_2015_ChIPseq.txt

#slurmMain QA --datadir=$PI_SCRATCH/Projects/biter_biter_stemCellFateRegulation/data/Liu_2015/ATACseq 
#slurmMain QA --datadir=$PI_SCRATCH/Projects/biter_biter_stemCellFateRegulation/data/DellOrso_2016_CellReports/ATACseq 
#slurmMain QA --datadir=$PI_SCRATCH/Projects/biter_biter_stemCellFateRegulation/data/Buenrostro_2013_NatureMethods/ATACseq 
#slurmMain QA --datadir=$PI_SCRATCH/Projects/biter_biter_stemCellFateRegulation/data/Brett_2015/RNAseq 
#slurmMain QA --datadir=$PI_SCRATCH/Projects/biter_biter_stemCellFateRegulation/data/Ikeda_2014/RNAseq --nc=16
#slurmMain QA --datadir=$PI_SCRATCH/Projects/biter_biter_stemCellFateRegulation/data/Liu_2014/RNAseq 

# 4 Ling - Hairless paper submission
#slurmMain QA --datadir=$PI_SCRATCH/Projects/biter_lingliu_DNAdamage/data/Liu_2016/ATACseq 
#slurmMain QA --datadir=$PI_SCRATCH/Projects/biter_lingliu_DNAdamage/data/Liu_2016/RNAseq 
#slurmMain QA --datadir=$PI_SCRATCH/Projects/biter_lingliu_DNAdamage/data/Liu_2016/ChIPseqV2 
#slurmMain QA --datadir=$PI_SCRATCH/Projects/biter_lingliu_DNAdamage/data/Liu_2016/ChIPseq2 

# 4 Ling
#slurmMain QA --datadir=$PI_SCRATCH/Projects/biter_lingliu_DNAdamage/data/Cheung_2013/RNAseq 
#slurmMain QA --datadir=$PI_SCRATCH/Projects/biter_lingliu_DNAdamage/data/Liu_2014/RNAseq 
#slurmMain QA --datadir=$PI_SCRATCH/Projects/biter_lingliu_DNAdamage/data/Liu_2014/MNaseseq 
#slurmMain QA --datadir=$PI_SCRATCH/Projects/biter_lingliu_DNAdamage/data/Liu_2011/ChIPseq 
#slurmMain QA --datadir=$PI_SCRATCH/Projects/biter_lingliu_DNAdamage/data/Liu_2012/ChIPseq 
#slurmMain QA --datadir=$PI_SCRATCH/Projects/biter_lingliu_DNAdamage/data/Liu_2013/ChIPseq 
##slurmMain QA --datadir=$PI_SCRATCH/Projects/biter_lingliu_DNAdamage/data/Liu_2015/ATACseq 
#slurmMain QA --datadir=$PI_SCRATCH/Projects/biter_lingliu_DNAdamage/data/Wenemoser_2015/ChIPseq
#slurmMain QA --datadir=$PI_SCRATCH/Projects/biter_lingliu_DNAdamage/data/Liu_2016/ChIPseq 
#slurmMain QA --datadir=$PI_SCRATCH/Projects/biter_lingliu_DNAdamage/data/Liu_2015/RNAseq 
#slurmMain QA --datadir=$PI_SCRATCH/Projects/biter_lingliu_DNAdamage/data/Liu_2016/RNAseq 
#slurmMain QA --datadir=$PI_SCRATCH/Projects/biter_lingliu_DNAdamage/data/Liu_2015/ChIPseq 
#slurmMain QA --datadir=$PI_SCRATCH/Projects/biter_lingliu_DNAdamage/data/Zhang_2015_Science/ChIPseq 
#slurmMain QA --datadir=$PI_SCRATCH/Projects/biter_lingliu_DNAdamage/data/Liu_2015/ChIPseq2 
#slurmMain QA --datadir=$PI_SCRATCH/Projects/biter_lingliu_DNAdamage/data/Liu_2016/ChIPseq2 
#slurmMain QA --datadir=$PI_SCRATCH/Projects/biter_lingliu_DNAdamage/data/Liu_2016/ChIPseqV2 
#slurmMain QA --datadir=$PI_SCRATCH/Projects/biter_lingliu_DNAdamage/data/Bernstein_2014/ChIPseq --nc=16

# 4 Cindy
#slurmMain QA --datadir=$PI_SCRATCH/Projects/biter_cvanvelt_SCTrxTrl/data/vanVelthoven_2015/PullDownTrx
#slurmMain QA --datadir=$PI_SCRATCH/Projects/biter_cvanvelt_SCTrxTrl/data/vanVelthoven_2015/PullDownTrl
#slurmMain QA --datadir=$PI_SCRATCH/Projects/biter_cvanvelt_SCTrxTrl/data/Gonzalez_2014_JNeurosci/PullDownTrl
#slurmMain QA --datadir=$PI_SCRATCH/Projects/biter_cvanvelt_SCTrxTrl/data/vanVelthoven_2015/PullDownTrx2
#slurmMain QA --datadir=$PI_SCRATCH/Projects/biter_cvanvelt_SCTrxTrl/data/vanVelthoven_2016/PullDownTrx

# 4 Mike
#slurmMain QA --datadir=$PI_SCRATCH/Projects/biter_wosczyna_miRNAsInMuscleCellFate/data/Wosczyna_2015/RNAseq 
#slurmMain QA --datadir=$PI_SCRATCH/Projects/biter_wosczyna_miRNAsInMuscleCellFate/data/Wosczyna_2015/PullDown


