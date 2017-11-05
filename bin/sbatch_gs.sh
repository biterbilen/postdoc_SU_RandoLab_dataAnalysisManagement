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
#slurmMain GS --datadir=$PI_SCRATCH/Projects/biter_lingliu_DNAdamage/data/Liu_2016/RNAseq --t=2 --nc=1 --parseOutput=$po --indir=DEQSC --otag=QSC
#slurmMain GS --datadir=$PI_SCRATCH/Projects/biter_lingliu_DNAdamage/data/Liu_2016/RNAseq --t=2 --nc=1 --parseOutput=$po --indir=DEASC --otag=ASC
#slurmMain GS --datadir=$PI_SCRATCH/Projects/biter_lingliu_DNAdamage/data/Liu_2016/RNAseq --t=2 --nc=1 --parseOutput=$po

# 4 Mike 2017
#slurmMain GS --datadir=$PI_SCRATCH/Projects/biter_wosczyna_miRNAsInMuscleCellFate/data/Wosczyna_2015/RNAseq --t=2 --nc=1 --parseOutput=$po --flecdir=PullDown
slurmMain GS --datadir=$PI_SCRATCH/Projects/biter_wosczyna_miRNAsInMuscleCellFate/data/Wosczyna_2015/PullDown --t=2 --nc=1 --parseOutput=$po --indir=DEtotRNA --otag=totRNA --flecdir=PullDown

#slurmMain GS --datadir=$PI_SCRATCH/Projects/biter_wosczyna_miRNAsInMuscleCellFate/data/Wosczyna_2015/PullDown --t=1 --nc=1 --parseOutput=$po --flecdir=PullDown
#slurmMain GS --datadir=$PI_SCRATCH/Projects/biter_wosczyna_miRNAsInMuscleCellFate/data/Wosczyna_2015/PullDown --t=1 --nc=1 --parseOutput=$po --indir=DEnulltotRNA --otag=nulltotRNA --flecdir=PullDown
#slurmMain GS --datadir=$PI_SCRATCH/Projects/biter_wosczyna_miRNAsInMuscleCellFate/data/Wosczyna_2015/PullDown --t=1 --nc=1 --parseOutput=$po --indir=DEpulRNA --otag=pulRNA --flecdir=PullDown # --dryrun=T

exit


# 4 Jamie Luiz Daniel 2017
#slurmMain GS --datadir=$PI_SCRATCH/Projects/biter_jbrett_DNAmethylation/data/Brett_2015/RNAseq --removeSpike=T --t=2 --parseOutput=$po 
#slurmMain GS --datadir=$PI_SCRATCH/Projects/biter_jbrett_DNAmethylation/data/Brett_2015/RNAseq --removeSpike=T --t=2 --otag=1sr --removeSamplePat=22m_NA_muscle_SC_18hGM_1 --parseOutput=$po 
#slurmMain GS --datadir=$PI_SCRATCH/Projects/biter_jbrett_DNAmethylation/data/Brett_2015/RNAseq --removeSpike=T --t=2 --otag=2sr --removeSamplePat="22m_NA_muscle_SC_18hGM_1|22m_NA_muscleTG_SC_NA_1" --parseOutput=$po 
#slurmMain GS --datadir=$PI_SCRATCH/Projects/biter_jbrett_DNAmethylation/data/Brett_2015/RNAseq --removeSpike=T --t=2 --otag=3sr --removeSamplePat="22m_NA_muscle_SC_18hGM_1|22m_NA_muscleTG_SC_NA_1|04m_NA_muscle_SC_NA_2" --parseOutput=$po 
#slurmMain GS --datadir=$PI_SCRATCH/Projects/biter_jbrett_DNAmethylation/data/Brett_2015/RNAseq --removeSpike=T --t=2 --otag=4sr --removeSamplePat="22m_NA_muscle_SC_18hGM_1|22m_NA_muscleTG_SC_NA_1|04m_NA_muscle_SC_NA_2|04m_NA_muscleTG_SC_NA_1" --parseOutput=$po 


# TODO organize
#slurmMain GS --datadir=$PI_SCRATCH/Projects/biter_jbrett_DNAmethylation/data/Ikeda_2014/RNAseq
#slurmMain GS --datadir=$PI_SCRATCH/Projects/biter_lingliu_DNAdamage/data/Liu_2015/ChIPseq
#slurmMain GS --datadir=$PI_SCRATCH/Projects/biter_lingliu_DNAdamage/data/Liu_2015/ATACseq
#slurmMain GS --datadir=$PI_SCRATCH/Projects/biter_lingliu_DNAdamage/data/Liu_2011/ChIPseq
#slurmMain GS --datadir=$PI_SCRATCH/Projects/biter_lingliu_DNAdamage/data/Liu_2015/RNAseq

# 4 Biter
#slurmMain GS --datadir=$PI_SCRATCH/Projects/biter_biter_stemCellFateRegulation/data/Ikeda_2014/RNAseq --parseOutput=$po 
#slurmMain GS --datadir=$PI_SCRATCH/Projects/biter_biter_stemCellFateRegulation/data/Brett_2015/RNAseq --removeSpike=T --parseOutput=$po 
#slurmMain GS --datadir=$PI_SCRATCH/Projects/biter_biter_stemCellFateRegulation/data/Liu_2014/RNAseq --parseOutput=$po

# 4 Ling - Hairless paper submission
#slurmMain GS --datadir=$PI_SCRATCH/Projects/biter_lingliu_DNAdamage/data/Liu_2016/RNAseq --parseOutput=$po

# 4 Jamie
#slurmMain GS --datadir=$PI_SCRATCH/Projects/biter_jbrett_DNAmethylation/data/Brett_2015/RNAseq 
#slurmMain GS --datadir=$PI_SCRATCH/Projects/biter_jbrett_DNAmethylation/data/Ikeda_2014/RNAseq 
#slurmMain GS --skippostproc=F --analysisdir=ER_Brett_2015_Ikeda_2014_RNAseq --datadir=$PI_SCRATCH/Projects/biter_jbrett_DNAmethylation/data/Brett_2015/RNAseq:$PI_SCRATCH/Projects/biter_jbrett_DNAmethylation/data/Ikeda_2014/RNAseq --headPatRep1=04m_NA_muscle_SC_NA_ --headPatRep2=04m_NA_muscle_SC_

# 4 Ling
#slurmMain GS --datadir=$PI_SCRATCH/Projects/biter_lingliu_DNAdamage/data/Cheung_2013/RNAseq
#slurmMain GS --datadir=$PI_SCRATCH/Projects/biter_lingliu_DNAdamage/data/Liu_2014/RNAseq
#slurmMain GS --skippostproc=F --analysisdir=ER_Liu_2014_Cheung_2013_RNAseq --datadir=$PI_SCRATCH/Projects/biter_lingliu_DNAdamage/data/Liu_2014/RNAseq:$PI_SCRATCH/Projects/biter_lingliu_DNAdamage/data/Cheung_2013/RNAseq --headPatRep1=NA_SC_ --headPatRep2=SC_
#slurmMain GS --skippostproc=F --analysisdir=ER_Liu_2015_Ikeda_2014_RNAseq --datadir=$PI_SCRATCH/Projects/biter_lingliu_DNAdamage/data/Liu_2015/RNAseq:$PI_SCRATCH/Projects/biter_jbrett_DNAmethylation/data/Ikeda_2014/RNAseq --headPatRep1=HrWT_SC_ --headPatRep2=04m_NA_muscle_SC_ --removeSpike=T

#slurmMain GS --datadir=$PI_SCRATCH/Projects/biter_lingliu_DNAdamage/data/Liu_2015/RNAseq --parseOutput=$po
#slurmMain GS --datadir=$PI_SCRATCH/Projects/biter_lingliu_DNAdamage/data/Liu_2016/RNAseq --parseOutput=$po
#slurmMain GS --skippostproc=F --analysisdir=ER_Liu_2015_2016_RNAseq --datadir=$PI_SCRATCH/Projects/biter_lingliu_DNAdamage/data/Liu_2015/RNAseq:$PI_SCRATCH/Projects/biter_lingliu_DNAdamage/data/Liu_2016/RNAseq --parseOutput=$po


# 4 Cindy
#slurmMain GS --datadir=$PI_SCRATCH/Projects/biter_cvanvelt_SCTrxTrl/data/vanVelthoven_2015/PullDownTrl
#slurmMain GS --datadir=$PI_SCRATCH/Projects/biter_cvanvelt_SCTrxTrl/data/vanVelthoven_2015/PullDownTrx
#slurmMain GS --skippostproc=F --analysisdir=ER_vanVelthoven_2015_PullDownTrx_PullDownTrl --datadir=$PI_SCRATCH/Projects/biter_cvanvelt_SCTrxTrl/data/vanVelthoven_2015/PullDownTrx:$PI_SCRATCH/Projects/biter_cvanvelt_SCTrxTrl/data/vanVelthoven_2015/PullDownTrl
#slurmMain GS --datadir=$PI_SCRATCH/Projects/biter_cvanvelt_SCTrxTrl/data/vanVelthoven_2015/PullDownTrx2 --indir=FC_RM_TR_ada --parseOutput=$po
#slurmMain GS --datadir=$PI_SCRATCH/Projects/biter_cvanvelt_SCTrxTrl/data/vanVelthoven_2016/PullDownTrx --indir=FC_RM_TR_ada --parseOutput=$po
#slurmMain GS --skippostproc=F --analysisdir=ER_FC_RM_TR_ada_vanVelthoven_2015_PullDownTrx2_2016_PullDownTrx --datadir=$PI_SCRATCH/Projects/biter_cvanvelt_SCTrxTrl/data/vanVelthoven_2015/PullDownTrx2:$PI_SCRATCH/Projects/biter_cvanvelt_SCTrxTrl/data/vanVelthoven_2016/PullDownTrx --indir=FC_RM_TR_ada --parseOutput=$po

# 4 Mike
#slurmMain GS --datadir=$PI_SCRATCH/Projects/biter_wosczyna_miRNAsInMuscleCellFate/data/Wosczyna_2015/RNAseq 
#slurmMain GS --datadir=$PI_SCRATCH/Projects/biter_wosczyna_miRNAsInMuscleCellFate/data/Wosczyna_2015/PullDown --removeSpike=F 
#slurmMain GS --datadir=$PI_SCRATCH/Projects/biter_wosczyna_miRNAsInMuscleCellFate/data/Wosczyna_2015/PullDown --removeSpike=T 
#slurmMain GS --skippostproc=F --analysisdir=ER_Wosczyna_2015_RNAseq_PullDown --datadir=$PI_SCRATCH/Projects/biter_wosczyna_miRNAsInMuscleCellFate/data/Wosczyna_2015/RNAseq:$PI_SCRATCH/Projects/biter_wosczyna_miRNAsInMuscleCellFate/data/Wosczyna_2015/PullDown --removeSpike=T


