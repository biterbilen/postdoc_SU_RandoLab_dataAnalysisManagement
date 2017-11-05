#!/bin/bash - 
#===============================================================================
#
#          FILE: sbatch_ds.sh
# 
#         USAGE: ./sbatch_ds.sh
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
regionfile=$PI_HOME/Data/casco/SRA/Wu_2017_BBA/bmsc_d00_15_histone_rx2_ctcf_smc1a_all_time_mm10_dense.bed.gz; linesSkipped=1
slurmMain DS --datadir=$PI_SCRATCH/Projects/biter_wosczyna_miRNAsInMuscleCellFate/data/Wosczyna_2015/RNAseq --regionsStr=self --t=3 --nc=2 --regionfile=$regionfile --linesSkipped=$linesSkipped --filterStr='state==4' --colsStr='chr,start,end,state' --posStrand=T --otag=Wu2017BBAstateStrongActiveEnh --parseOutput=$po #--debug=T --dryrun=T
slurmMain DS --datadir=$PI_SCRATCH/Projects/biter_wosczyna_miRNAsInMuscleCellFate/data/Wosczyna_2015/RNAseq --regionsStr=self --t=3 --nc=2 --regionfile=$regionfile --linesSkipped=$linesSkipped --filterStr='state==5' --colsStr='chr,start,end,state' --posStrand=T --otag=Wu2017BBAstateActiveEnh --parseOutput=$po #--debug=T --dryrun=T

#slurmMain DS --datadir=$PI_SCRATCH/Projects/biter_wosczyna_miRNAsInMuscleCellFate/data/Wosczyna_2015/RNAseq --regionsStr=self --t=3 --nc=2 --regionfile=$regionfile --linesSkipped=$linesSkipped --filterStr='state==7' --colsStr='chr,start,end,state' --posStrand=T --otag=Wu2017BBAstateCTCFmedEnh --parseOutput=$po #--debug=T --dryrun=T
#slurmMain DS --datadir=$PI_SCRATCH/Projects/biter_wosczyna_miRNAsInMuscleCellFate/data/Wosczyna_2015/RNAseq --regionsStr=self --t=3 --nc=2 --regionfile=$regionfile --linesSkipped=$linesSkipped --filterStr='state==8' --colsStr='chr,start,end,state' --posStrand=T --otag=Wu2017BBAstateRunx2medEnh --parseOutput=$po #--debug=T --dryrun=T
#slurmMain DS --datadir=$PI_SCRATCH/Projects/biter_wosczyna_miRNAsInMuscleCellFate/data/Wosczyna_2015/RNAseq --regionsStr=self --t=3 --nc=2 --regionfile=$regionfile --linesSkipped=$linesSkipped --filterStr='state==9' --colsStr='chr,start,end,state' --posStrand=T --otag=Wu2017BBAstateCTCFCohesinIns --parseOutput=$po #--debug=T --dryrun=T
