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

source $PI_SCRATCH/Projects/biter_biter_shared/bin/auxslurm.sh
po=${1:-F}

# 4 Ling 2017
#slurmMain ER --datadir=$PI_SCRATCH/Projects/biter_lingliu_DNAdamage/data/Liu_2016/ATACseq --t=2 --nc=1 --parseOutput=$po
#slurmMain ER --datadir=$PI_SCRATCH/Projects/biter_lingliu_DNAdamage/data/Liu_2017/ATACseq --t=2 --nc=1 --parseOutput=$po
#slurmMain ER --datadir=$PI_SCRATCH/Projects/biter_lingliu_DNAdamage/data/Liu_2016/ChIPseqV2 --t=2 --nc=1 --parseOutput=$po
#slurmMain ER --skippostproc=T --analysisdir=ER_Liu_2016_ATACseq_FC_FCoc_2016_ChIPseqV2_FC --datadir=$PI_SCRATCH/Projects/biter_lingliu_DNAdamage/data/Liu_2017/ATACseq:/scratch/PI/casco/Projects/biter_lingliu_DNAdamage/data/Liu_2017/ATACseq:/scratch/PI/casco/Projects/biter_lingliu_DNAdamage/data/Liu_2016/ChIPseqV2 --indir=FC:FCoc:FC --nameCol=Name:Name:Name --parseOutput=$po 
#slurmMain ER --datadir=$PI_SCRATCH/Projects/biter_lingliu_DNAdamage/data/Liu_2016/RNAseq --parseOutput=$po
slurmMain ER --datadir=$PI_SCRATCH/Projects/biter_lingliu_DNAdamage/data/Liu_2015/RNAseq --parseOutput=$po

# 4 Mike 2017
#slurmMain ER --datadir=$PI_SCRATCH/Projects/biter_wosczyna_miRNAsInMuscleCellFate/data/Wosczyna_2015/RNAseq --t=2 --parseOutput=$po
#slurmMain ER --datadir=$PI_SCRATCH/Projects/biter_wosczyna_miRNAsInMuscleCellFate/data/Wosczyna_2015/PullDown --removeSpike=T --t=2 --parseOutput=$po

# 4 Jamie Luiz Daniel 2017
#slurmMain ER --datadir=$PI_SCRATCH/Projects/biter_jbrett_DNAmethylation/data/Brett_2015/RNAseq --removeSpike=T --t=2 --parseOutput=$po 
#slurmMain ER --datadir=$PI_SCRATCH/Projects/biter_jbrett_DNAmethylation/data/Brett_2015/RNAseq --removeSpike=T --t=2 --otag=1sr --removeSamplePat=22m_NA_muscle_SC_18hGM_1 --parseOutput=$po 
#slurmMain ER --datadir=$PI_SCRATCH/Projects/biter_jbrett_DNAmethylation/data/Brett_2015/RNAseq --removeSpike=T --t=2 --otag=2sr --removeSamplePat="22m_NA_muscle_SC_18hGM_1|22m_NA_muscleTG_SC_NA_1" --parseOutput=$po 
#slurmMain ER --datadir=$PI_SCRATCH/Projects/biter_jbrett_DNAmethylation/data/Brett_2015/RNAseq --removeSpike=T --t=2 --otag=3sr --removeSamplePat="22m_NA_muscle_SC_18hGM_1|22m_NA_muscleTG_SC_NA_1|04m_NA_muscle_SC_NA_2" --parseOutput=$po 
#slurmMain ER --datadir=$PI_SCRATCH/Projects/biter_jbrett_DNAmethylation/data/Brett_2015/RNAseq --removeSpike=T --t=2 --otag=4sr --removeSamplePat="22m_NA_muscle_SC_18hGM_1|22m_NA_muscleTG_SC_NA_1|04m_NA_muscle_SC_NA_2|04m_NA_muscleTG_SC_NA_1" --parseOutput=$po 


# TODO organize
#slurmMain ER --datadir=$PI_SCRATCH/Projects/biter_jbrett_DNAmethylation/data/Ikeda_2014/RNAseq
#slurmMain ER --datadir=$PI_SCRATCH/Projects/biter_lingliu_DNAdamage/data/Liu_2015/ChIPseq
#slurmMain ER --datadir=$PI_SCRATCH/Projects/biter_lingliu_DNAdamage/data/Liu_2015/ATACseq
#slurmMain ER --datadir=$PI_SCRATCH/Projects/biter_lingliu_DNAdamage/data/Liu_2011/ChIPseq
#slurmMain ER --datadir=$PI_SCRATCH/Projects/biter_lingliu_DNAdamage/data/Liu_2015/RNAseq

# 4 Biter
#slurmMain ER --datadir=$PI_SCRATCH/Projects/biter_biter_stemCellFateRegulation/data/Ikeda_2014/RNAseq --parseOutput=$po 
#slurmMain ER --datadir=$PI_SCRATCH/Projects/biter_biter_stemCellFateRegulation/data/Brett_2015/RNAseq --removeSpike=T --parseOutput=$po 
#slurmMain ER --datadir=$PI_SCRATCH/Projects/biter_biter_stemCellFateRegulation/data/Liu_2014/RNAseq --parseOutput=$po

# 4 Ling - Hairless paper submission
#slurmMain ER --datadir=$PI_SCRATCH/Projects/biter_lingliu_DNAdamage/data/Liu_2016/RNAseq --parseOutput=$po

# 4 Jamie
#slurmMain ER --datadir=$PI_SCRATCH/Projects/biter_jbrett_DNAmethylation/data/Brett_2015/RNAseq 
#slurmMain ER --datadir=$PI_SCRATCH/Projects/biter_jbrett_DNAmethylation/data/Ikeda_2014/RNAseq 
#slurmMain ER --skippostproc=F --analysisdir=ER_Brett_2015_Ikeda_2014_RNAseq --datadir=$PI_SCRATCH/Projects/biter_jbrett_DNAmethylation/data/Brett_2015/RNAseq:$PI_SCRATCH/Projects/biter_jbrett_DNAmethylation/data/Ikeda_2014/RNAseq --headPatRep1=04m_NA_muscle_SC_NA_ --headPatRep2=04m_NA_muscle_SC_

# 4 Ling
#slurmMain ER --datadir=$PI_SCRATCH/Projects/biter_lingliu_DNAdamage/data/Cheung_2013/RNAseq
#slurmMain ER --datadir=$PI_SCRATCH/Projects/biter_lingliu_DNAdamage/data/Liu_2014/RNAseq
#slurmMain ER --skippostproc=F --analysisdir=ER_Liu_2014_Cheung_2013_RNAseq --datadir=$PI_SCRATCH/Projects/biter_lingliu_DNAdamage/data/Liu_2014/RNAseq:$PI_SCRATCH/Projects/biter_lingliu_DNAdamage/data/Cheung_2013/RNAseq --headPatRep1=NA_SC_ --headPatRep2=SC_
#slurmMain ER --skippostproc=F --analysisdir=ER_Liu_2015_Ikeda_2014_RNAseq --datadir=$PI_SCRATCH/Projects/biter_lingliu_DNAdamage/data/Liu_2015/RNAseq:$PI_SCRATCH/Projects/biter_jbrett_DNAmethylation/data/Ikeda_2014/RNAseq --headPatRep1=HrWT_SC_ --headPatRep2=04m_NA_muscle_SC_ --removeSpike=T

#slurmMain ER --datadir=$PI_SCRATCH/Projects/biter_lingliu_DNAdamage/data/Liu_2015/RNAseq --parseOutput=$po
#slurmMain ER --datadir=$PI_SCRATCH/Projects/biter_lingliu_DNAdamage/data/Liu_2016/RNAseq --parseOutput=$po
#slurmMain ER --skippostproc=F --analysisdir=ER_Liu_2015_2016_RNAseq --datadir=$PI_SCRATCH/Projects/biter_lingliu_DNAdamage/data/Liu_2015/RNAseq:$PI_SCRATCH/Projects/biter_lingliu_DNAdamage/data/Liu_2016/RNAseq --parseOutput=$po


# 4 Cindy
#slurmMain ER --datadir=$PI_SCRATCH/Projects/biter_cvanvelt_SCTrxTrl/data/vanVelthoven_2015/PullDownTrl
#slurmMain ER --datadir=$PI_SCRATCH/Projects/biter_cvanvelt_SCTrxTrl/data/vanVelthoven_2015/PullDownTrx
#slurmMain ER --skippostproc=F --analysisdir=ER_vanVelthoven_2015_PullDownTrx_PullDownTrl --datadir=$PI_SCRATCH/Projects/biter_cvanvelt_SCTrxTrl/data/vanVelthoven_2015/PullDownTrx:$PI_SCRATCH/Projects/biter_cvanvelt_SCTrxTrl/data/vanVelthoven_2015/PullDownTrl
#slurmMain ER --datadir=$PI_SCRATCH/Projects/biter_cvanvelt_SCTrxTrl/data/vanVelthoven_2015/PullDownTrx2 --indir=FC_RM_TR_ada --parseOutput=$po
#slurmMain ER --datadir=$PI_SCRATCH/Projects/biter_cvanvelt_SCTrxTrl/data/vanVelthoven_2016/PullDownTrx --indir=FC_RM_TR_ada --parseOutput=$po
#slurmMain ER --skippostproc=F --analysisdir=ER_FC_RM_TR_ada_vanVelthoven_2015_PullDownTrx2_2016_PullDownTrx --datadir=$PI_SCRATCH/Projects/biter_cvanvelt_SCTrxTrl/data/vanVelthoven_2015/PullDownTrx2:$PI_SCRATCH/Projects/biter_cvanvelt_SCTrxTrl/data/vanVelthoven_2016/PullDownTrx --indir=FC_RM_TR_ada --parseOutput=$po

# 4 Mike
#slurmMain ER --datadir=$PI_SCRATCH/Projects/biter_wosczyna_miRNAsInMuscleCellFate/data/Wosczyna_2015/RNAseq 
#slurmMain ER --datadir=$PI_SCRATCH/Projects/biter_wosczyna_miRNAsInMuscleCellFate/data/Wosczyna_2015/PullDown --removeSpike=F 
#slurmMain ER --datadir=$PI_SCRATCH/Projects/biter_wosczyna_miRNAsInMuscleCellFate/data/Wosczyna_2015/PullDown --removeSpike=T 
#slurmMain ER --skippostproc=F --analysisdir=ER_Wosczyna_2015_RNAseq_PullDown --datadir=$PI_SCRATCH/Projects/biter_wosczyna_miRNAsInMuscleCellFate/data/Wosczyna_2015/RNAseq:$PI_SCRATCH/Projects/biter_wosczyna_miRNAsInMuscleCellFate/data/Wosczyna_2015/PullDown --removeSpike=T


