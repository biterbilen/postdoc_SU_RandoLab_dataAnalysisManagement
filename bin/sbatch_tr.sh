#!/bin/bash - 
#===============================================================================
#
#          FILE: sbatch_tr.sh
# 
#         USAGE: ./sbatch_tr.sh 
# 
#   DESCRIPTION: 
# 
#       OPTIONS: ---
#  REQUIREMENTS: ---
#          BUGS: ---
#         NOTES: ---
#        AUTHOR: Dr. Biter Bilen (bb), biterbilen@yahoo.com
#  ORGANIZATION: Stanford, Palo Alto, USA
#       CREATED: 08/13/2015 11:19
#      REVISION:  ---
#===============================================================================

source $PI_SCRATCH/Projects/biter_biter_shared/bin/auxslurm.sh
po=${1:-F}

# 4 Fosl1 project
#slurmMain TR --datadir=$PI_SCRATCH/Projects/biter_demorree_Fosl1/data/Wold_2011/ChIPseq --t=12 \
#	--as="TruSeqIndexedAdaptor:AGATCGGAAGAGC" --As="TruSeqIndexedAdaptor:AGATCGGAAGAGC" --n=2 --parseOutput=$po

# 4 Mike 2017
#slurmMain TR --datadir=$PI_SCRATCH/Projects/biter_wosczyna_miRNAsInMuscleCellFate/data/Wosczyna_2015/PullDown --t=12 --n=2 \
#	--as="TruSeqIndexedAdaptor:AGATCGGAAGAGC" --As="TruSeqIndexedAdaptor:AGATCGGAAGAGC" --parseOutput=$po #--debug=T
#slurmMain TR --datadir=$PI_SCRATCH/Projects/biter_wosczyna_miRNAsInMuscleCellFate/data/Wosczyna_2015/RNAseq --t=12 --n=2 \
#	--as="TruSeqIndexedAdaptor:AGATCGGAAGAGC" --As="TruSeqIndexedAdaptor:AGATCGGAAGAGC" --parseOutput=$po #--debug=T

# 4 Ling 2017
#slurmMain TR --datadir=$PI_SCRATCH/Projects/biter_lingliu_DNAdamage/data/Liu_2017/ATACseq --t=12 --n=2 \
#	--as="NexteraTransposaseSeq:CTGTCTCTTATA" --As="NexteraTransposaseSeq:CTGTCTCTTATA" --n=2 --parseOutput=$po
#slurmMain TR --datadir=$PI_SCRATCH/Projects/biter_lingliu_DNAdamage/data/Liu_2016/ATACseq --t=12 --n=2 \
#	--as="NexteraTransposaseSeq:CTGTCTCTTATA" --As="NexteraTransposaseSeq:CTGTCTCTTATA" --n=2 --parseOutput=$po
#slurmMain TR --datadir=$PI_SCRATCH/Projects/biter_lingliu_DNAdamage/data/Liu_2016/RNAseq --t=12 \
#	--as="TruSeqIndexedAdaptor:AGATCGGAAGAGC" --As="TruSeqIndexedAdaptor:AGATCGGAAGAGC" --n=2 --parseOutput=$po
slurmMain TR --datadir=$PI_SCRATCH/Projects/biter_lingliu_DNAdamage/data/Liu_2015/RNAseq --t=12 \
	--as="TruSeqIndexedAdaptor:AGATCGGAAGAGC" --As="TruSeqIndexedAdaptor:AGATCGGAAGAGC" --n=2 --parseOutput=$po

# 4 Cindy 2017
#slurmMain TR --datadir=$PI_SCRATCH/Projects/biter_cvanvelt_SCTrxTrl/data/vanVelthoven_2016/PullDownTrx2 --t=12 --n=2 \
#	--as="TruSeqIndexedAdaptor:AGATCGGAAGAGC" --As="TruSeqIndexedAdaptor:AGATCGGAAGAGC" --parseOutput=$po #--debug=T

# 4 Jamie Luiz Daniel 2017
#slurmMain TR --datadir=$PI_SCRATCH/Projects/biter_jbrett_DNAmethylation/data/Brett_2015/RNAseq --t=12 --n=2 \
#	--as="TruSeqIndexedAdaptor:AGATCGGAAGAGC" --As="TruSeqIndexedAdaptor:AGATCGGAAGAGC" --skipfastqc=T --parseOutput=$po #--debug=T

# 4 Antoine
#slurmMain TR --nc=2 --t=2 --datadir=$PI_SCRATCH/Projects/biter_demorree_PAS/data/deMorree_2012/directRNAseq --otag=trim5termT \
#	--trim5termN=T --skipcutadapt=T --n=2 --parseOutput=$po
#slurmMain TR --nc=2 --t=2 --datadir=$PI_SCRATCH/Projects/biter_demorree_PAS/data/deMorree_2012/directRNAseq --otag=trim5termTanypolyA \
#	--trim5termN=T --as="polyA:'A{20}'" --n=2 --parseOutput=$po
#slurmMain TR --nc=2 --t=2 --datadir=$PI_SCRATCH/Projects/biter_demorree_PAS/data/deMorree_2012/directRNAseq --otag=trim5termTany10A \
#	--trim5termN=T --as="polyA:'A{10}'" --n=2 --parseOutput=$po

# 4 Jamie
#slurmMain TR --datadir=$PI_SCRATCH/Projects/biter_jbrett_DNAmethylation/data/Brett_2014/WGBseq --t=30 --nc=16 \
#	--as="TruSeqIndexedAdaptor:AGATCGGAAGAGC" --As="TruSeqIndexedAdaptor:AGATCGGAAGAGC" --n=2 --parseOutput=$po

# 4 Jay
#slurmMain TR --datadir=$PI_SCRATCH/Projects/biter_jsalvi_spuriousReplication/data/Salvi_2016/ChIPseq \
#	--as="TruSeqIndexedAdaptor:AGATCGGAAGAGC" --n=2 --parseOutput=$po

# 4 Biter
#slurmMain TR --datadir=$PI_SCRATCH/Projects/biter_biter_stemCellFateRegulation/data/Liu_2015/ATACseq --t=12 \
#	--as="NexteraTransposaseSeq:CTGTCTCTTATA" --As="NexteraTransposaseSeq:CTGTCTCTTATA" --n=2 --parseOutput=$po
#slurmMain TR --datadir=$PI_SCRATCH/Projects/biter_biter_stemCellFateRegulation/data/DellOrso_2016_CellReports/ATACseq --t=12 \
#	--otag=ada --as="NexteraTransposaseSeq:CTGTCTCTTATA" --As="NexteraTransposaseSeq:CTGTCTCTTATA" --n=2 --parseOutput=$po

# 4 Ling - Hairless paper submission
#slurmMain TR --datadir=$PI_SCRATCH/Projects/biter_lingliu_DNAdamage/data/Liu_2016/ATACseq --t=12 \
#	--as="NexteraTransposaseSeq:CTGTCTCTTATA" --As="NexteraTransposaseSeq:CTGTCTCTTATA" --n=2 --parseOutput=$po
#slurmMain TR --datadir=$PI_SCRATCH/Projects/biter_lingliu_DNAdamage/data/Liu_2016/RNAseq --t=12 \
#	--as="TruSeqIndexedAdaptor:AGATCGGAAGAGC" --As="TruSeqIndexedAdaptor:AGATCGGAAGAGC" --n=2 --parseOutput=$po

# For Ling
#slurmMain TR --datadir=$PI_SCRATCH/Projects/biter_lingliu_DNAdamage/data/Liu_2014/MNaseseq --t=12 \
#	--otag=ada --as="TruSeqIndexedAdaptor:AGATCGGAAGAGC" --n=2
#slurmMain TR --datadir=$PI_SCRATCH/Projects/biter_lingliu_DNAdamage/data/Wenemoser_2015/ChIPseq --t=12 \
#	--otag=ada --as="TruSeqIndexedAdaptor:AGATCGGAAGAGC" --n=2
#slurmMain TR --datadir=$PI_SCRATCH/Projects/biter_lingliu_DNAdamage/data/Liu_2011/ChIPseq --t=12 \
#	--otag=ada --as="TruSeqIndexedAdaptor:AGATCGGAAGAGC" --n=2
#slurmMain TR --datadir=$PI_SCRATCH/Projects/biter_lingliu_DNAdamage/data/Liu_2012/ChIPseq --t=12 \
#	--otag=ada --as="TruSeqIndexedAdaptor:AGATCGGAAGAGC" --n=2
#slurmMain TR --datadir=$PI_SCRATCH/Projects/biter_lingliu_DNAdamage/data/Liu_2013/ChIPseq --t=12 \
#	--otag=ada --as="TruSeqIndexedAdaptor:AGATCGGAAGAGC" --n=2
##slurmMain TR --datadir=$PI_SCRATCH/Projects/biter_lingliu_DNAdamage/data/Liu_2015/ATACseq --t=12 \
##	--otag=ada --as="NexteraTransposaseSeq:CTGTCTCTTATA" --As="NexteraTransposaseSeq:CTGTCTCTTATA" --n=2
#slurmMain TR --datadir=$PI_SCRATCH/Projects/biter_lingliu_DNAdamage/data/Liu_2015/RNAseq --t=12 \
#	--otag=ada --as="TruSeqIndexedAdaptor:AGATCGGAAGAGC" --As="TruSeqIndexedAdaptor:AGATCGGAAGAGC" --n=2
#slurmMain TR --datadir=$PI_SCRATCH/Projects/biter_lingliu_DNAdamage/data/Zhang_2015_Science/ChIPseq --t=16 \
#	--otag=ada --as="TruSeqIndexedAdaptor:AGATCGGAAGAGC" --As="TruSeqIndexedAdaptor:AGATCGGAAGAGC" --n=2 --nc=16 --parseOutput=$po
#slurmMain TR --datadir=$PI_SCRATCH/Projects/biter_lingliu_DNAdamage/data/Liu_2015/ChIPseq2 --t=12 \
#	--otag=ada --as="TruSeqIndexedAdaptor:AGATCGGAAGAGC" --n=2 --parseOutput=$po

# 4 Cindy
#slurmMain TR --datadir=$PI_SCRATCH/Projects/biter_cvanvelt_SCTrxTrl/data/vanVelthoven_2015/PullDownTrl --t=12 \
#	--otag=ada --as="TruSeqIndexedAdaptor:AGATCGGAAGAGC" --As="TruSeqIndexedAdaptor:AGATCGGAAGAGC" --n=2
#slurmMain TR --datadir=$PI_SCRATCH/Projects/biter_cvanvelt_SCTrxTrl/data/vanVelthoven_2015/PullDownTrx --t=12 \
#	--otag=ada --as="TruSeqIndexedAdaptor:AGATCGGAAGAGC" --As="TruSeqIndexedAdaptor:AGATCGGAAGAGC" --n=2
#slurmMain TR --datadir=$PI_SCRATCH/Projects/biter_cvanvelt_SCTrxTrl/data/Gonzalez_2014_JNeurosci/PullDownTrl --t=12 \
#	--otag=ada --as="RiboProfiling3pAdaptor:CTGTAGGCACCAT;A100_1:'A{100}';T100_1:'T{100}'" --n=2
#slurmMain TR --datadir=$PI_SCRATCH/Projects/biter_cvanvelt_SCTrxTrl/data/vanVelthoven_2015/PullDownTrx2 --t=16 \
#	--otag=ada --as="TruSeqIndexedAdaptor:AGATCGGAAGAGC" --As="TruSeqIndexedAdaptor:AGATCGGAAGAGC" --n=2 --nc=16 --parseOutput=$po
#slurmMain TR --datadir=$PI_SCRATCH/Projects/biter_cvanvelt_SCTrxTrl/data/vanVelthoven_2016/PullDownTrx --t=16 \
#	--otag=ada --as="TruSeqIndexedAdaptor:AGATCGGAAGAGC" --As="TruSeqIndexedAdaptor:AGATCGGAAGAGC" --n=2 --nc=16 --parseOutput=$po

# 4 Mike
#slurmMain TR --datadir=$PI_SCRATCH/Projects/biter_wosczyna_miRNAsInMuscleCellFate/data/Wosczyna_2015/PullDown --t=12 \
#	--otag=ada_A100_T100 --as="TruSeqIndexedAdaptor:AGATCGGAAGAGC;A100_1:'A{100}';T100_1:'T{100}'" --As="TruSeqIndexedAdaptor:AGATCGGAAGAGC;A100_2:'A{100}';T100_2:'T{100}'" --n=2
#slurmMain TR --datadir=$PI_SCRATCH/Projects/biter_wosczyna_miRNAsInMuscleCellFate/data/Wosczyna_2015/RNAseq --t=12 \
#	--otag=ada_A100_T100 --as="TruSeqIndexedAdaptor:AGATCGGAAGAGC;A100_1:'A{100}';T100_1:'T{100}'" --As="TruSeqIndexedAdaptor:AGATCGGAAGAGC;A100_2:'A{100}';T100_2:'T{100}'" --n=2

