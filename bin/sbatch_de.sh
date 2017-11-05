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
po=${1:-F}

# 4 Mike 2017
#slurmMain DE --datadir=$PI_SCRATCH/Projects/biter_wosczyna_miRNAsInMuscleCellFate/data/Wosczyna_2015/RNAseq --nc=2 --t=1 --contrast='c(1,-1)' --parseOutput=$po
slurmMain DE --datadir=$PI_SCRATCH/Projects/biter_wosczyna_miRNAsInMuscleCellFate/data/Wosczyna_2015/PullDown --nc=2 --t=1 --contrast='c(1,0,-1,0)' --otag=nulltotRNA --parseOutput=$po
#slurmMain DE --datadir=$PI_SCRATCH/Projects/biter_wosczyna_miRNAsInMuscleCellFate/data/Wosczyna_2015/PullDown --nc=2 --t=1 --contrast='c(1,-1,-1,1)' --parseOutput=$po
#slurmMain DE --datadir=$PI_SCRATCH/Projects/biter_wosczyna_miRNAsInMuscleCellFate/data/Wosczyna_2015/PullDown --nc=2 --t=1 --contrast='c(1,-1)' --otag=totRNA --colselpat=QiagenRNEasy --bamsuf=.nsrt.bam --parseOutput=$po
#slurmMain DE --datadir=$PI_SCRATCH/Projects/biter_wosczyna_miRNAsInMuscleCellFate/data/Wosczyna_2015/PullDown --nc=2 --t=1 --contrast='c(1,-1)' --otag=pulRNA --colselpat=MWmodPD --bamsuf=.nsrt.bam --parseOutput=$po

# 4 Ling 2017
#slurmMain DE --datadir=$PI_SCRATCH/Projects/biter_lingliu_DNAdamage/data/Liu_2016/RNAseq --nc=2 --t=1 --contrast='c(1,-1)' --otag=ASC --colselpat=_TMX05hdPI_ --bamsuf=.nsrt.bam --parseOutput=$po
#slurmMain DE --datadir=$PI_SCRATCH/Projects/biter_lingliu_DNAdamage/data/Liu_2016/RNAseq --nc=2 --t=1 --contrast='c(1,-1)' --otag=QSC --colselpat=_TMX_ --bamsuf=.nsrt.bam  --parseOutput=$po
#slurmMain DE --datadir=$PI_SCRATCH/Projects/biter_lingliu_DNAdamage/data/Liu_2015/RNAseq --nc=2 --t=1 --contrast='c(1,-1)' --parseOutput=$po

exit
# Integrate
#slurmMain DE --datadir=$PI_SCRATCH/Projects/biter_jbrett_DNAmethylation/data/Brett_2015/RNAseq 
#slurmMain DE --datadir=. --analysisdir=TEMP --skippostproc=F
#slurmMain DE --datadir=. --analysisdir=GMMfixed --skippostproc=F
#slurmMain DE --datadir=. --analysisdir=OUT.none --skippostproc=F
#slurmMain DE --datadir=. --analysisdir=OUT.scale --skippostproc=F

#for scale in none upper; do
for scale in upper; do

	# 4 Jamie Luiz Daniel 2017
	k=1; #k=1; #k=3
	# order yASC yRSC yEASC yQSC_WM yQSC oASC oRSC oEASC oQSC_WM oQSC
	#TODO check the contrasts
	contrasts='c(0,0,0,0,-1,0,0,0,1)'; tags="'OvYquiescence'" # 5 9 
	contrasts+=',c(0,0,0,-1,0,0,0,1,0)'; tags+=",'OvYquiescenceWM'" # 4 8 
	contrasts+=',c(-1,0,0,0,0,1,0,0,0)'; tags+=",'OvYactivation'" # 6 1
	contrasts+=',c(0,0,-1,0,0,0,1,0,0)'; tags+=",'OvYearlyactivation'" # 7 3
	contrasts+=',c(0,1,0,0,-1,0,0,0,0)'; tags+=",'cautionRvQyoung'" # 2 5 
	contrasts+=',c(0,0,0,0,0,1,0,0,-1)'; tags+=",'cautionAvQold'" # 6 9 
	contrasts+=',c(1,0,0,0,-1,0,0,0,0)'; tags+=",'cautionAvQyoung'" # 1 5 
	contrasts+=',c(0,0,1,-1,0,0,0,0,0)'; tags+=",'cautionEAvQyoung'" # 1 5 
	contrasts+=',c(1,0,0,0,0,1,0,-1,0)'; tags+=",'cautionEAvQold'" # 1 5 
	slurmMain DE --datadir=$PI_SCRATCH/Projects/biter_jbrett_DNAmethylation/data/Brett_2015/RNAseq --indir=EG1sr --scale=$scale --k=$k --otag=k${k}_sc${scale} --contrast="$contrasts" --tag="$tags" --nc=2

	# 4 Biter
	# TODO change contrast for Brett_2015
	#slurmMain DE --datadir=$PI_SCRATCH/Projects/biter_biter_stemCellFateRegulation/data/Brett_2015/RNAseq  --scale=$scale --otag=$scale --contrast='c(1,0,0,-1,0,0,0,0,0)'
	#slurmMain DE --datadir=$PI_SCRATCH/Projects/biter_biter_stemCellFateRegulation/data/Ikeda_2014/RNAseq  --scale=$scale --otag=$scale --contrast='c(0,0,-1,0,1)'

	# 4 Ling
	#slurmMain DE --datadir=$PI_SCRATCH/Projects/biter_lingliu_DNAdamage/results/2015-10-12 --scale=$scale --analysisdir=DE_$scale\_Liu_2014_Cheung_2013_RNAseq --indir=EG_Liu_2014_Cheung_2013_RNAseq --contrast='c(-1,1,0,0,0)' --skippostproc=F 
	#slurmMain DE --datadir=$PI_SCRATCH/Projects/biter_lingliu_DNAdamage/results/2015-10-12 --scale=$scale --analysisdir=DE_$scale\_Liu_2015_Ikeda_2014_RNAseq --indir=EG_Liu_2015_Ikeda_2014_RNAseq --contrast='c(1,0,0,-1,0,0)' --skippostproc=F 
	
	#slurmMain DE --datadir=$PI_SCRATCH/Projects/biter_lingliu_DNAdamage/results/2015-10-12 --scale=$scale --analysisdir=DE_$scale\_Liu_2015_2016_RNAseq --indir=EG_Liu_2015_2016_RNAseq --contrast='c(1,0,-1,0),c(0,1,0,-1)' --skippostproc=F 
	#slurmMain DE --datadir=$PI_SCRATCH/Projects/biter_lingliu_DNAdamage/data/Liu_2015/RNAseq  --scale=$scale --otag=$scale
	#slurmMain DE --datadir=$PI_SCRATCH/Projects/biter_lingliu_DNAdamage/data/Liu_2016/RNAseq  --scale=$scale --otag=$scale --contrast='c(1,0,-1,0),c(0,1,0,-1)'

	# 4 Mike
	#slurmMain DE --datadir=$PI_SCRATCH/Projects/biter_wosczyna_miRNAsInMuscleCellFate/data/Wosczyna_2015/RNAseq  --scale=$scale --otag=$scale
	#slurmMain DE --datadir=$PI_SCRATCH/Projects/biter_wosczyna_miRNAsInMuscleCellFate/data/Wosczyna_2015/PullDown --scale=$scale --otag=$scale --contrast='c(1,-1,-1,1)'

	# 4 Jamie
#	slurmMain DE --datadir=$PI_SCRATCH/Projects/biter_jbrett_DNAmethylation/data/Brett_2015/RNAseq  --scale=$scale --otag=$scale --contrast='c(1,0,0,-1,0,0,0,0,0)'
#	slurmMain DE --datadir=$PI_SCRATCH/Projects/biter_jbrett_DNAmethylation/data/Ikeda_2014/RNAseq  --scale=$scale --otag=$scale --contrast='c(0,0,-1,0,1)'
#	slurmMain DE --datadir=$PI_SCRATCH/Projects/biter_jbrett_DNAmethylation/results/2015-09-14 --scale=$scale --analysisdir=DE_$scale\_Brett_2015_Ikeda_2014_RNAseq --indir=EG_Brett_2015_Ikeda_2014_RNAseq --contrast='c(0,0,0,0,-1,0,1,0,0,0,0,0,0)' --skippostproc=F 

	# 4 Cindy
	#slurmMain DE --datadir=$PI_SCRATCH/Projects/biter_cvanvelt_SCTrxTrl/data/vanVelthoven_2015/PullDownTrl --scale=$scale --otag=$scale --contrast='c(-1,1)'
	#slurmMain DE --datadir=$PI_SCRATCH/Projects/biter_cvanvelt_SCTrxTrl/data/vanVelthoven_2015/PullDownTrx --scale=$scale --otag=$scale --contrast='c(0,0,0,1,-1)'
	#slurmMain DE --datadir=$PI_SCRATCH/Projects/biter_cvanvelt_SCTrxTrl/results/2015-11-06 --scale=$scale --analysisdir=DE_$scale\_vanVelthoven_2015_PullDownTrx_PullDownTrl --indir=EG_vanVelthoven_2015_PullDownTrx_PullDownTrl --contrast='c(0,0,0,1,0,-1,0)' --skippostproc=F
	#slurmMain DE --datadir=$PI_SCRATCH/Projects/biter_cvanvelt_SCTrxTrl/results/2016-03-15 --scale=$scale --analysisdir=DE_EG_FC_RM_TR_ada_$scale\_vanVelthoven_2015_PullDownTrx2_2016_PullDownTrx --indir=EG_FC_RM_TR_ada_vanVelthoven_2015_PullDownTrx2_2016_PullDownTrx --contrast='c(0,1,0,-1,0,0)' --skippostproc=F --removeSpike=T
#	slurmMain DE --datadir=$PI_SCRATCH/Projects/biter_cvanvelt_SCTrxTrl/data/vanVelthoven_2015/PullDownTrx2 --scale=$scale --otag=$scale --contrast='c(1,-1,0)' --removeSpike=T --indir=EG_FC_RM_TR_ada
#	slurmMain DE --datadir=$PI_SCRATCH/Projects/biter_cvanvelt_SCTrxTrl/data/vanVelthoven_2016/PullDownTrx --scale=$scale --otag=$scale --contrast='c(0,1,0,-1,0,0)' --removeSpike=T --indir=EG_FC_RM_TR_ada

done


