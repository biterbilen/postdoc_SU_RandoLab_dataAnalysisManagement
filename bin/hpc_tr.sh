#!/bin/bash - 
#===============================================================================
#
#          FILE: hpc_tr.sh
# 
#         USAGE: ./hpc_tr.sh 
# 
#   DESCRIPTION: 
# 
#       OPTIONS: ---
#  REQUIREMENTS: ---
#          BUGS: ---
#         NOTES: ---
#        AUTHOR: Dr. Biter Bilen (bb), biterbilen@yahoo.com
#  ORGANIZATION: Stanford, Palo Alto, USA
#       CREATED: 11/16/2017 22:20
#      REVISION:  ---
#===============================================================================

set -o nounset # Treat unset variables as an error

SRCDIR=$HOME/Projects/postdoc_SU_RandoLab_dataAnalysisManagement/bin
source $SRCDIR/auxHPC.sh

po=${1:-F}

# 4 Antoine 2017
Main TR --srcdir=$SRCDIR --datadir=$HOME/Projects/biter_demorree_PAS/data/deMorree_2012/directRNAseq --nc=2 --t=2 --scheduler=sge \
  --otag=trim5termT --trim5termN=T --skipcutadapt=T --n=2 --parseOutput=$po

# 4 Fosl1 project
#Main TR --srcdir=$SRCDIR --datadir==$HOME/Projects/biter_demorree_Fosl1/data/Wold_2011/ChIPseq --t=12 \
#	--as="TruSeqIndexedAdaptor:AGATCGGAAGAGC" --As="TruSeqIndexedAdaptor:AGATCGGAAGAGC" --n=2 --parseOutput=$po

# 4 Mike 2017
#Main TR --srcdir=$SRCDIR --datadir==$HOME/Projects/biter_wosczyna_miRNAsInMuscleCellFate/data/Wosczyna_2015/PullDown --t=12 --n=2 \
#	--as="TruSeqIndexedAdaptor:AGATCGGAAGAGC" --As="TruSeqIndexedAdaptor:AGATCGGAAGAGC" --parseOutput=$po #--debug=T
#Main TR --srcdir=$SRCDIR --datadir==$HOME/Projects/biter_wosczyna_miRNAsInMuscleCellFate/data/Wosczyna_2015/RNAseq --t=12 --n=2 \
#	--as="TruSeqIndexedAdaptor:AGATCGGAAGAGC" --As="TruSeqIndexedAdaptor:AGATCGGAAGAGC" --parseOutput=$po #--debug=T

# 4 Ling 2017
#Main TR --srcdir=$SRCDIR --datadir==$HOME/Projects/biter_lingliu_DNAdamage/data/Liu_2017/ATACseq --t=12 --n=2 \
#	--as="NexteraTransposaseSeq:CTGTCTCTTATA" --As="NexteraTransposaseSeq:CTGTCTCTTATA" --n=2 --parseOutput=$po
#Main TR --srcdir=$SRCDIR --datadir==$HOME/Projects/biter_lingliu_DNAdamage/data/Liu_2016/ATACseq --t=12 --n=2 \
#	--as="NexteraTransposaseSeq:CTGTCTCTTATA" --As="NexteraTransposaseSeq:CTGTCTCTTATA" --n=2 --parseOutput=$po
#Main TR --srcdir=$SRCDIR --datadir==$HOME/Projects/biter_lingliu_DNAdamage/data/Liu_2016/RNAseq --t=12 \
#	--as="TruSeqIndexedAdaptor:AGATCGGAAGAGC" --As="TruSeqIndexedAdaptor:AGATCGGAAGAGC" --n=2 --parseOutput=$po
#Main TR --srcdir=$SRCDIR --datadir==$HOME/Projects/biter_lingliu_DNAdamage/data/Liu_2015/RNAseq --t=12 \
#	--as="TruSeqIndexedAdaptor:AGATCGGAAGAGC" --As="TruSeqIndexedAdaptor:AGATCGGAAGAGC" --n=2 --parseOutput=$po

# 4 Cindy 2017
#Main TR --srcdir=$SRCDIR --datadir==$HOME/Projects/biter_cvanvelt_SCTrxTrl/data/vanVelthoven_2016/PullDownTrx2 --t=12 --n=2 \
#	--as="TruSeqIndexedAdaptor:AGATCGGAAGAGC" --As="TruSeqIndexedAdaptor:AGATCGGAAGAGC" --parseOutput=$po #--debug=T

# 4 Jamie Luiz Daniel 2017
#Main TR --srcdir=$SRCDIR --datadir==$HOME/Projects/biter_jbrett_DNAmethylation/data/Brett_2015/RNAseq --t=12 --n=2 \
#	--as="TruSeqIndexedAdaptor:AGATCGGAAGAGC" --As="TruSeqIndexedAdaptor:AGATCGGAAGAGC" --skipfastqc=T --parseOutput=$po #--debug=T

# 4 Antoine
#Main TR --srcdir=$SRCDIR --datadir=$HOME/Projects/biter_demorree_PAS/data/deMorree_2012/directRNAseq --nc=2 --t=2 --otag=trim5termT \
#	--trim5termN=T --skipcutadapt=T --n=2 --parseOutput=$po
#Main TR --srcdir=$SRCDIR --datadir=$HOME/Projects/biter_demorree_PAS/data/deMorree_2012/directRNAseq --nc=2 --t=2 --otag=trim5termTanypolyA \
#	--trim5termN=T --as="polyA:'A{20}'" --n=2 --parseOutput=$po
#Main TR --srcdir=$SRCDIR --datadir=$HOME/Projects/biter_demorree_PAS/data/deMorree_2012/directRNAseq --nc=2 --t=2 --otag=trim5termTany10A \
#	--trim5termN=T --as="polyA:'A{10}'" --n=2 --parseOutput=$po

# 4 Jamie
#Main TR --srcdir=$SRCDIR --datadir==$HOME/Projects/biter_jbrett_DNAmethylation/data/Brett_2014/WGBseq --t=30 --nc=16 \
#	--as="TruSeqIndexedAdaptor:AGATCGGAAGAGC" --As="TruSeqIndexedAdaptor:AGATCGGAAGAGC" --n=2 --parseOutput=$po

# 4 Jay
#Main TR --srcdir=$SRCDIR --datadir==$HOME/Projects/biter_jsalvi_spuriousReplication/data/Salvi_2016/ChIPseq \
#	--as="TruSeqIndexedAdaptor:AGATCGGAAGAGC" --n=2 --parseOutput=$po

# 4 Biter
#Main TR --srcdir=$SRCDIR --datadir==$HOME/Projects/biter_biter_stemCellFateRegulation/data/Liu_2015/ATACseq --t=12 \
#	--as="NexteraTransposaseSeq:CTGTCTCTTATA" --As="NexteraTransposaseSeq:CTGTCTCTTATA" --n=2 --parseOutput=$po
#Main TR --srcdir=$SRCDIR --datadir==$HOME/Projects/biter_biter_stemCellFateRegulation/data/DellOrso_2016_CellReports/ATACseq --t=12 \
#	--otag=ada --as="NexteraTransposaseSeq:CTGTCTCTTATA" --As="NexteraTransposaseSeq:CTGTCTCTTATA" --n=2 --parseOutput=$po

# 4 Ling - Hairless paper submission
#Main TR --srcdir=$SRCDIR --datadir==$HOME/Projects/biter_lingliu_DNAdamage/data/Liu_2016/ATACseq --t=12 \
#	--as="NexteraTransposaseSeq:CTGTCTCTTATA" --As="NexteraTransposaseSeq:CTGTCTCTTATA" --n=2 --parseOutput=$po
#Main TR --srcdir=$SRCDIR --datadir==$HOME/Projects/biter_lingliu_DNAdamage/data/Liu_2016/RNAseq --t=12 \
#	--as="TruSeqIndexedAdaptor:AGATCGGAAGAGC" --As="TruSeqIndexedAdaptor:AGATCGGAAGAGC" --n=2 --parseOutput=$po

# For Ling
#Main TR --srcdir=$SRCDIR --datadir==$HOME/Projects/biter_lingliu_DNAdamage/data/Liu_2014/MNaseseq --t=12 \
#	--otag=ada --as="TruSeqIndexedAdaptor:AGATCGGAAGAGC" --n=2
#Main TR --srcdir=$SRCDIR --datadir==$HOME/Projects/biter_lingliu_DNAdamage/data/Wenemoser_2015/ChIPseq --t=12 \
#	--otag=ada --as="TruSeqIndexedAdaptor:AGATCGGAAGAGC" --n=2
#Main TR --srcdir=$SRCDIR --datadir==$HOME/Projects/biter_lingliu_DNAdamage/data/Liu_2011/ChIPseq --t=12 \
#	--otag=ada --as="TruSeqIndexedAdaptor:AGATCGGAAGAGC" --n=2
#Main TR --srcdir=$SRCDIR --datadir==$HOME/Projects/biter_lingliu_DNAdamage/data/Liu_2012/ChIPseq --t=12 \
#	--otag=ada --as="TruSeqIndexedAdaptor:AGATCGGAAGAGC" --n=2
#Main TR --srcdir=$SRCDIR --datadir==$HOME/Projects/biter_lingliu_DNAdamage/data/Liu_2013/ChIPseq --t=12 \
#	--otag=ada --as="TruSeqIndexedAdaptor:AGATCGGAAGAGC" --n=2
##Main TR --srcdir=$SRCDIR --datadir==$HOME/Projects/biter_lingliu_DNAdamage/data/Liu_2015/ATACseq --t=12 \
##	--otag=ada --as="NexteraTransposaseSeq:CTGTCTCTTATA" --As="NexteraTransposaseSeq:CTGTCTCTTATA" --n=2
#Main TR --srcdir=$SRCDIR --datadir==$HOME/Projects/biter_lingliu_DNAdamage/data/Liu_2015/RNAseq --t=12 \
#	--otag=ada --as="TruSeqIndexedAdaptor:AGATCGGAAGAGC" --As="TruSeqIndexedAdaptor:AGATCGGAAGAGC" --n=2
#Main TR --srcdir=$SRCDIR --datadir==$HOME/Projects/biter_lingliu_DNAdamage/data/Zhang_2015_Science/ChIPseq --t=16 \
#	--otag=ada --as="TruSeqIndexedAdaptor:AGATCGGAAGAGC" --As="TruSeqIndexedAdaptor:AGATCGGAAGAGC" --n=2 --nc=16 --parseOutput=$po
#Main TR --srcdir=$SRCDIR --datadir==$HOME/Projects/biter_lingliu_DNAdamage/data/Liu_2015/ChIPseq2 --t=12 \
#	--otag=ada --as="TruSeqIndexedAdaptor:AGATCGGAAGAGC" --n=2 --parseOutput=$po

# 4 Cindy
#Main TR --srcdir=$SRCDIR --datadir==$HOME/Projects/biter_cvanvelt_SCTrxTrl/data/vanVelthoven_2015/PullDownTrl --t=12 \
#	--otag=ada --as="TruSeqIndexedAdaptor:AGATCGGAAGAGC" --As="TruSeqIndexedAdaptor:AGATCGGAAGAGC" --n=2
#Main TR --srcdir=$SRCDIR --datadir==$HOME/Projects/biter_cvanvelt_SCTrxTrl/data/vanVelthoven_2015/PullDownTrx --t=12 \
#	--otag=ada --as="TruSeqIndexedAdaptor:AGATCGGAAGAGC" --As="TruSeqIndexedAdaptor:AGATCGGAAGAGC" --n=2
#Main TR --srcdir=$SRCDIR --datadir==$HOME/Projects/biter_cvanvelt_SCTrxTrl/data/Gonzalez_2014_JNeurosci/PullDownTrl --t=12 \
#	--otag=ada --as="RiboProfiling3pAdaptor:CTGTAGGCACCAT;A100_1:'A{100}';T100_1:'T{100}'" --n=2
#Main TR --srcdir=$SRCDIR --datadir==$HOME/Projects/biter_cvanvelt_SCTrxTrl/data/vanVelthoven_2015/PullDownTrx2 --t=16 \
#	--otag=ada --as="TruSeqIndexedAdaptor:AGATCGGAAGAGC" --As="TruSeqIndexedAdaptor:AGATCGGAAGAGC" --n=2 --nc=16 --parseOutput=$po
#Main TR --srcdir=$SRCDIR --datadir==$HOME/Projects/biter_cvanvelt_SCTrxTrl/data/vanVelthoven_2016/PullDownTrx --t=16 \
#	--otag=ada --as="TruSeqIndexedAdaptor:AGATCGGAAGAGC" --As="TruSeqIndexedAdaptor:AGATCGGAAGAGC" --n=2 --nc=16 --parseOutput=$po

# 4 Mike
#Main TR --srcdir=$SRCDIR --datadir==$HOME/Projects/biter_wosczyna_miRNAsInMuscleCellFate/data/Wosczyna_2015/PullDown --t=12 \
#	--otag=ada_A100_T100 --as="TruSeqIndexedAdaptor:AGATCGGAAGAGC;A100_1:'A{100}';T100_1:'T{100}'" --As="TruSeqIndexedAdaptor:AGATCGGAAGAGC;A100_2:'A{100}';T100_2:'T{100}'" --n=2
#Main TR --srcdir=$SRCDIR --datadir==$HOME/Projects/biter_wosczyna_miRNAsInMuscleCellFate/data/Wosczyna_2015/RNAseq --t=12 \
#	--otag=ada_A100_T100 --as="TruSeqIndexedAdaptor:AGATCGGAAGAGC;A100_1:'A{100}';T100_1:'T{100}'" --As="TruSeqIndexedAdaptor:AGATCGGAAGAGC;A100_2:'A{100}';T100_2:'T{100}'" --n=2

