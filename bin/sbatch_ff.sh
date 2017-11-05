#!/bin/bash - 
#===============================================================================
#
#          FILE: sbatch_ff.sh
# 
#         USAGE: ./sbatch_ff.sh
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

##slurmMain FF --datadir=$PI_SCRATCH/Projects/biter_wosczyna_miRNAsInMuscleCellFate/data/Wosczyna_2015/RNAseq --regionsStr=promoters --t=3 --nc=2 --parseOutput=$po
slurmMain FF --datadir=$PI_SCRATCH/Projects/biter_wosczyna_miRNAsInMuscleCellFate/data/Wosczyna_2015/RNAseq --regionsStr=promoters --t=2 --nc=6 --minScore="90%" --otag=PWMms90Prom --parseOutput=$po
regionfile=$PI_SCRATCH/Projects/biter_wosczyna_miRNAsInMuscleCellFate/data/Wosczyna_2015/RNAseq/DSWu2017BBAstateStrongActiveEnh/self.fa.gz; otag=PWMms90DSWu2017BBAstateStrongActiveEnh #45' run time for 21K regions
slurmMain FF --datadir=$PI_SCRATCH/Projects/biter_wosczyna_miRNAsInMuscleCellFate/data/Wosczyna_2015/RNAseq --regionsStr=self --regionfile=$regionfile --t=2 --nc=6 --minScore="90%" --otag=$otag --parseOutput=$po
##slurmMain FF --datadir=$PI_SCRATCH/Projects/biter_wosczyna_miRNAsInMuscleCellFate/data/Wosczyna_2015/RNAseq --regionsStr=promoters --t=12 --nc=6 --minScore="80%" --parseOutput=$po
##slurmMain FF --datadir=$PI_SCRATCH/Projects/biter_wosczyna_miRNAsInMuscleCellFate/data/Wosczyna_2015/RNAseq --regionsStr=promoters --t=12 --nc=6 --minScore="85%" --otag=minScore85 --parseOutput=$po

miRExpfile=$PI_SCRATCH/Projects/biter_wosczyna_miRNAsInMuscleCellFate/data/Wosczyna_2013/selected.txt.gz
#slurmMain FF --datadir=$PI_SCRATCH/Projects/biter_wosczyna_miRNAsInMuscleCellFate/data/Wosczyna_2015/PullDown --regionsStr=threeUTRsByTranscript --t=1 --nc=1 --otag=miRseed7mer3UTR --miRExpfile=$miRExpfile --parseOutput=$po
#slurmMain FF --datadir=$PI_SCRATCH/Projects/biter_wosczyna_miRNAsInMuscleCellFate/data/Wosczyna_2015/PullDown --regionsStr=cdsBy --t=1 --nc=1									--otag=miRseed7merCDS --miRExpfile=$miRExpfile --parseOutput=$po
##slurmMain FF --datadir=$PI_SCRATCH/Projects/biter_wosczyna_miRNAsInMuscleCellFate/data/Wosczyna_2015/PullDown --regionsStr=threeUTRsByTranscript --t=1 --nc=1             --otag=adipogenesisMiRseed7mer --miRExpfile=$miRExpfile --parseOutput=$po
##slurmMain FF --datadir=$PI_SCRATCH/Projects/biter_wosczyna_miRNAsInMuscleCellFate/data/Wosczyna_2015/PullDown --regionsStr=threeUTRsByTranscript --t=1 --nc=1 --seedEnd=7 --otag=adipogenesisMiRseed6mer --miRExpfile=$miRExpfile --parseOutput=$po
##slurmMain FF --datadir=$PI_SCRATCH/Projects/biter_wosczyna_miRNAsInMuscleCellFate/data/Wosczyna_2015/PullDown --regionsStr=threeUTRsByTranscript --t=1 --nc=1             --otag=seed7mer --parseOutput=$po
##slurmMain FF --datadir=$PI_SCRATCH/Projects/biter_wosczyna_miRNAsInMuscleCellFate/data/Wosczyna_2015/PullDown --regionsStr=threeUTRsByTranscript --t=1 --nc=1 --seedEnd=7 --otag=seed6mer --parseOutput=$po

exit

# 4 Ling 
slurmMain FF --datadir=~/Projects/biter_lingliu_DNAdamage/data/Liu_2016/RNAseq --indir=FC --infile=featureCounts.txt.gz --analysisdir=FF_final --skippostproc=F --regionsStr=promoters,threeUTRsByTranscript

# 4 Mike
#slurmMain FF --datadir=~/Projects/biter_wosczyna_miRNAsInMuscleCellFate/results/2015-07-17 --indir=EG_Wosczyna_2015_RNAseq_PullDown --analysisdir=FF_final --skippostproc=F --miRexpfile=~/Projects/biter_wosczyna_miRNAsInMuscleCellFate/data/Wosczyna_2013/selected3



