#!/bin/bash - 
#===============================================================================
#
#          FILE: hpc_qa.sh
# 
#         USAGE: ./hpc_qa.sh
# 
#   DESCRIPTION: 
# 
#       OPTIONS: ---
#  REQUIREMENTS: ---
#          BUGS: ---
#         NOTES: ---
#        AUTHOR: Dr. Biter Bilen (bb), biterbilen@yahoo.com
#  ORGANIZATION: Stanford, Palo Alto, USA
#       CREATED: 11/13/2017 09:51
#      REVISION:  ---
#===============================================================================

set -o nounset # Treat unset variables as an error

SRCDIR=$HOME/Projects/postdoc_SU_RandoLab_dataAnalysisManagement/bin
source $SRCDIR/auxHPC.sh

# 4 Antoine 2017
Main QA --srcdir=$SRCDIR --datadir=$HOME/Projects/biter_demorree_PAS/data/deMorree_2012/directRNAseq --pattern=.fasta.gz --nc=2 --t=2 --scheduler=sge --dryrun=F

# 4 Fosl1 project
#Main QA --srcdir=$SRCDIR --datadir==$HOME/Projects/biter_demorree_Fosl1/data/Wold_2011/ChIPseq

# 4 Mike 2017/04
#Main QA --srcdir=$SRCDIR --datadir==$HOME/Projects/biter_wosczyna_miRNAsInMuscleCellFate/data/Wosczyna_2015/RNAseq
#Main QA --srcdir=$SRCDIR --datadir==$HOME/Projects/biter_wosczyna_miRNAsInMuscleCellFate/data/Wosczyna_2015/PullDown

# 4 Ling 2017
#Main QA --srcdir=$SRCDIR --datadir==$HOME/Projects/biter_lingliu_DNAdamage/data/Liu_2016/ATACseq 
#Main QA --srcdir=$SRCDIR --datadir==$HOME/Projects/biter_lingliu_DNAdamage/data/Liu_2017/ATACseq
#Main QA --srcdir=$SRCDIR --datadir==$HOME/Projects/biter_lingliu_DNAdamage/data/Liu_2016/ChIPseqV2 
#Main QA --srcdir=$SRCDIR --datadir==$HOME/Projects/biter_lingliu_DNAdamage/data/Liu_2016/RNAseq --nc=16
#Main QA --srcdir=$SRCDIR --datadir==$HOME/Projects/biter_lingliu_DNAdamage/data/Liu_2015/RNAseq --nc=16

# 4 Cindy 2017
#Main QA --srcdir=$SRCDIR --datadir==$HOME/Projects/biter_cvanvelt_SCTrxTrl/data/vanVelthoven_2016/PullDownTrx2

# 4 Jamie Luiz Daniel 2017
#Main QA --srcdir=$SRCDIR --datadir==$HOME/Projects/biter_jbrett_DNAmethylation/data/Brett_2015/RNAseq --skipfastqc=T
#Main QA --srcdir=$SRCDIR --datadir==$HOME/Projects/biter_jbrett_DNAmethylation/data/Ikeda_2014/RNAseq --skipfastqc=T
#Main QA --srcdir=$SRCDIR --datadir==$HOME/Projects/biter_lingliu_DNAdamage/data/Liu_2014/RNAseq --skipfastqc=T

# 4 Jamie
#Main QA --srcdir=$SRCDIR --datadir==$HOME/Projects/biter_jbrett_DNAmethylation/data/Brett_2014/WGBseq --nc=8 --t=8 --skipfastqc=T
##Main QA --srcdir=$SRCDIR --datadir==$HOME/Projects/biter_jbrett_DNAmethylation/data/Brett_2014/WGBseq --nc=16 --t=24
#Main QA --srcdir=$SRCDIR --datadir==$HOME/Projects/biter_jbrett_DNAmethylation/data/Brett_2015/RNAseq 
#Main QA --srcdir=$SRCDIR --datadir==$HOME/Projects/biter_jbrett_DNAmethylation/data/Ikeda_2014/RNAseq 

# 4 Jay
#Main QA --srcdir=$SRCDIR --datadir==$HOME/Projects/biter_jsalvi_spuriousReplication/data/Salvi_2016/ChIPseq
#Main QA --srcdir=$SRCDIR --datadir==$HOME/Projects/biter_jsalvi_spuriousReplication/data/Dellino_2013_GenomeResearch/ChIPseq

# 4 Antoine
#Main QA --srcdir=$SRCDIR --datadir==$HOME/Projects/biter_demorree_PAS/data/deMorree_2012/directRNAseq --pattern=.fasta.gz --nc=2 --t=2
#Main QA --srcdir=$SRCDIR --datadir==$HOME/Projects/biter_demorree_PAS/data/Soleimani_2012_DevCell/ChIPseq

# 4 Biter
#Main QA --srcdir=$SRCDIR --datadir==$HOME/Projects/biter_biter_noiseInAging/data/Guess_2015/RNAseq
#Main QA --srcdir=$SRCDIR --datadir==$HOME/Projects/biter_biter_EWSR1/data/Bilen_2015/ChIPseq --limitsfile=$PI_SCRATCH/Projects/biter_biter_EWSR1/results/2015-12-21/limits_Bilen_2015_ChIPseq.txt

#Main QA --srcdir=$SRCDIR --datadir==$HOME/Projects/biter_biter_stemCellFateRegulation/data/Liu_2015/ATACseq 
#Main QA --srcdir=$SRCDIR --datadir==$HOME/Projects/biter_biter_stemCellFateRegulation/data/DellOrso_2016_CellReports/ATACseq 
#Main QA --srcdir=$SRCDIR --datadir==$HOME/Projects/biter_biter_stemCellFateRegulation/data/Buenrostro_2013_NatureMethods/ATACseq 
#Main QA --srcdir=$SRCDIR --datadir==$HOME/Projects/biter_biter_stemCellFateRegulation/data/Brett_2015/RNAseq 
#Main QA --srcdir=$SRCDIR --datadir==$HOME/Projects/biter_biter_stemCellFateRegulation/data/Ikeda_2014/RNAseq --nc=16
#Main QA --srcdir=$SRCDIR --datadir==$HOME/Projects/biter_biter_stemCellFateRegulation/data/Liu_2014/RNAseq 

# 4 Ling - Hairless paper submission
#Main QA --srcdir=$SRCDIR --datadir==$HOME/Projects/biter_lingliu_DNAdamage/data/Liu_2016/ATACseq 
#Main QA --srcdir=$SRCDIR --datadir==$HOME/Projects/biter_lingliu_DNAdamage/data/Liu_2016/RNAseq 
#Main QA --srcdir=$SRCDIR --datadir==$HOME/Projects/biter_lingliu_DNAdamage/data/Liu_2016/ChIPseqV2 
#Main QA --srcdir=$SRCDIR --datadir==$HOME/Projects/biter_lingliu_DNAdamage/data/Liu_2016/ChIPseq2 

# 4 Ling
#Main QA --srcdir=$SRCDIR --datadir==$HOME/Projects/biter_lingliu_DNAdamage/data/Cheung_2013/RNAseq 
#Main QA --srcdir=$SRCDIR --datadir==$HOME/Projects/biter_lingliu_DNAdamage/data/Liu_2014/RNAseq 
#Main QA --srcdir=$SRCDIR --datadir==$HOME/Projects/biter_lingliu_DNAdamage/data/Liu_2014/MNaseseq 
#Main QA --srcdir=$SRCDIR --datadir==$HOME/Projects/biter_lingliu_DNAdamage/data/Liu_2011/ChIPseq 
#Main QA --srcdir=$SRCDIR --datadir==$HOME/Projects/biter_lingliu_DNAdamage/data/Liu_2012/ChIPseq 
#Main QA --srcdir=$SRCDIR --datadir==$HOME/Projects/biter_lingliu_DNAdamage/data/Liu_2013/ChIPseq 
##Main QA --srcdir=$SRCDIR --datadir==$HOME/Projects/biter_lingliu_DNAdamage/data/Liu_2015/ATACseq 
#Main QA --srcdir=$SRCDIR --datadir==$HOME/Projects/biter_lingliu_DNAdamage/data/Wenemoser_2015/ChIPseq
#Main QA --srcdir=$SRCDIR --datadir==$HOME/Projects/biter_lingliu_DNAdamage/data/Liu_2016/ChIPseq 
#Main QA --srcdir=$SRCDIR --datadir==$HOME/Projects/biter_lingliu_DNAdamage/data/Liu_2015/RNAseq 
#Main QA --srcdir=$SRCDIR --datadir==$HOME/Projects/biter_lingliu_DNAdamage/data/Liu_2016/RNAseq 
#Main QA --srcdir=$SRCDIR --datadir==$HOME/Projects/biter_lingliu_DNAdamage/data/Liu_2015/ChIPseq 
#Main QA --srcdir=$SRCDIR --datadir==$HOME/Projects/biter_lingliu_DNAdamage/data/Zhang_2015_Science/ChIPseq 
#Main QA --srcdir=$SRCDIR --datadir==$HOME/Projects/biter_lingliu_DNAdamage/data/Liu_2015/ChIPseq2 
#Main QA --srcdir=$SRCDIR --datadir==$HOME/Projects/biter_lingliu_DNAdamage/data/Liu_2016/ChIPseq2 
#Main QA --srcdir=$SRCDIR --datadir==$HOME/Projects/biter_lingliu_DNAdamage/data/Liu_2016/ChIPseqV2 
#Main QA --srcdir=$SRCDIR --datadir==$HOME/Projects/biter_lingliu_DNAdamage/data/Bernstein_2014/ChIPseq --nc=16

# 4 Cindy
#Main QA --srcdir=$SRCDIR --datadir==$HOME/Projects/biter_cvanvelt_SCTrxTrl/data/vanVelthoven_2015/PullDownTrx
#Main QA --srcdir=$SRCDIR --datadir==$HOME/Projects/biter_cvanvelt_SCTrxTrl/data/vanVelthoven_2015/PullDownTrl
#Main QA --srcdir=$SRCDIR --datadir==$HOME/Projects/biter_cvanvelt_SCTrxTrl/data/Gonzalez_2014_JNeurosci/PullDownTrl
#Main QA --srcdir=$SRCDIR --datadir==$HOME/Projects/biter_cvanvelt_SCTrxTrl/data/vanVelthoven_2015/PullDownTrx2
#Main QA --srcdir=$SRCDIR --datadir==$HOME/Projects/biter_cvanvelt_SCTrxTrl/data/vanVelthoven_2016/PullDownTrx

# 4 Mike
#Main QA --srcdir=$SRCDIR --datadir==$HOME/Projects/biter_wosczyna_miRNAsInMuscleCellFate/data/Wosczyna_2015/RNAseq 
#Main QA --srcdir=$SRCDIR --datadir==$HOME/Projects/biter_wosczyna_miRNAsInMuscleCellFate/data/Wosczyna_2015/PullDown


