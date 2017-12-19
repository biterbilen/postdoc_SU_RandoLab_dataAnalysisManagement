#!/bin/bash - 
#==============================================================================
#
#          FILE: hpc_rm.sh
# 
#         USAGE: ./hpc_rm.sh 
# 
#   DESCRIPTION: 
# 
#       OPTIONS: ---
#  REQUIREMENTS: ---
#          BUGS: ---
#         NOTES: ---
#        AUTHOR: Dr. Biter Bilen (bb), biterbilen@yahoo.com
#  ORGANIZATION: Stanford, Palo Alto, USA
#       CREATED: 11/16/2017 22:30
#      REVISION:  ---
#==============================================================================

set -o nounset # Treat unset variables as an error

SRCDIR=$HOME/Projects/postdoc_SU_RandoLab_dataAnalysisManagement/bin
GENOME=mm10
GENIND=$HOME/Data/HISAT2_indexes/$GENOME/genome
source $SRCDIR/auxHPC.sh

po=${1:-F}

# 4 Antoine 2017
# TODO: spliced=T
Main RM --genomeindex=$GENIND --srcdir=$SRCDIR --datadir=$HOME/Projects/biter_demorree_PAS/data/deMorree_2012/directRNAseq --scheduler=sge --scoremin=L,0.0,-0.6 --indir=TRtrim5termT --parseOutput=$po 

# 4 Fosl1 project
#Main RM --genomeindex=$GENIND --srcdir=$SRCDIR --datadir=$HOME/Projects/biter_demorree_Fosl1/data/Wold_2011/ChIPseq --indir=TR --rmdup=T --parseOutput=$po 

# 4 Mike 2017
#Main RM --genomeindex=$GENIND --srcdir=$SRCDIR --datadir=$HOME/Projects/biter_wosczyna_miRNAsInMuscleCellFate/data/Wosczyna_2015/RNAseq --spliced=T --indir=TR --parseOutput=$po
#Main RM --genomeindex=$GENIND --srcdir=$SRCDIR --datadir=$HOME/Projects/biter_wosczyna_miRNAsInMuscleCellFate/data/Wosczyna_2015/PullDown --spliced=T --indir=TR --genome=mm10_ERCC92 --parseOutput=$po

# 4 Ling 2017
#Main RM --genomeindex=$GENIND --srcdir=$SRCDIR --datadir=$HOME/Projects/biter_lingliu_DNAdamage/data/Liu_2017/ATACseq --maxins=2000 --indir=TR --parseOutput=$po
#Main RM --genomeindex=$GENIND --srcdir=$SRCDIR --datadir=$HOME/Projects/biter_lingliu_DNAdamage/data/Liu_2016/ATACseq --maxins=2000 --indir=TR --parseOutput=$po
#Main RM --genomeindex=$GENIND --srcdir=$SRCDIR --datadir=$HOME/Projects/biter_lingliu_DNAdamage/data/Liu_2016/RNAseq --spliced=T --indir=TR --parseOutput=$po
#Main RM --genomeindex=$GENIND --srcdir=$SRCDIR --datadir=$HOME/Projects/biter_lingliu_DNAdamage/data/Liu_2015/RNAseq --spliced=T --indir=TR --parseOutput=$po
#Main RM --genomeindex=$GENIND --srcdir=$SRCDIR --datadir=$HOME/Projects/biter_lingliu_DNAdamage/data/Liu_2016/ChIPseqV2 --rmdup=T --parseOutput=$po 

# 4 Cindy 2017
#Main RM --genomeindex=$GENIND --srcdir=$SRCDIR --datadir=$HOME/Projects/biter_cvanvelt_SCTrxTrl/data/vanVelthoven_2016/PullDownTrx2 --indir=TR --spliced=T --nc=16 --parseOutput=$po

# 4 Jamie Luiz Daniel 2017
#Main RM --genomeindex=$GENIND --srcdir=$SRCDIR --datadir=$HOME/Projects/biter_jbrett_DNAmethylation/data/Brett_2015/RNAseq --indir=TR --spliced=T --genome=mm10_ERCC92 --nc=16 --parseOutput=$po #--debug=T #--dryrun=T

# 4 Jay
#Main RM --genomeindex=$GENIND --srcdir=$SRCDIR --datadir=$HOME/Projects/biter_jsalvi_spuriousReplication/data/Salvi_2016/ChIPseq --parseOutput=$po 
#Main RM --genomeindex=$GENIND --srcdir=$SRCDIR --datadir=$HOME/Projects/biter_jsalvi_spuriousReplication/data/Salvi_2016/ChIPseq --indir=TR --parseOutput=$po
#Main RM --genomeindex=$GENIND --srcdir=$SRCDIR --datadir=$HOME/Projects/biter_jsalvi_spuriousReplication/data/Dellino_2013_GenomeResearch/ChIPseq --genome=hg38 --parseOutput=$po

# 4 Antoine
#Main RM --genomeindex=$GENIND --srcdir=$SRCDIR --datadir=$HOME/Projects/biter_demorree_PAS/data/Soleimani_2012_DevCell/ChIPseq --parseOutput=$po 
##Main RM --genomeindex=$GENIND --srcdir=$SRCDIR --datadir=$HOME/Projects/biter_demorree_PAS/data/deMorree_2012/directRNAseq --scoremin=L,0.0,-0.4 --otag=L04minscore --indir=TR_trim5termT --parseOutput=$po 
#Main RM --genomeindex=$GENIND --srcdir=$SRCDIR --datadir=$HOME/Projects/biter_demorree_PAS/data/deMorree_2012/directRNAseq --scoremin=L,0.0,-0.4 --otag=L04minscoreNSC --indir=TR_trim5termT --nosoftclipping=T --parseOutput=$po 
#Main RM --genomeindex=$GENIND --srcdir=$SRCDIR --datadir=$HOME/Projects/biter_demorree_PAS/data/deMorree_2012/directRNAseq --scoremin=L,0.0,-0.6 --otag=L06minscore --indir=TR_trim5termT --parseOutput=$po 
#Main RM --genomeindex=$GENIND --srcdir=$SRCDIR --datadir=$HOME/Projects/biter_demorree_PAS/data/deMorree_2012/directRNAseq --scoremin=L,0.0,-0.6 --otag=L06minscore --indir=TR_trim5termTanypolyA --parseOutput=$po 
#Main RM --genomeindex=$GENIND --srcdir=$SRCDIR --datadir=$HOME/Projects/biter_demorree_PAS/data/deMorree_2012/directRNAseq --scoremin=L,0.0,-0.6 --otag=L06minscore --indir=TR_trim5termTany10A --parseOutput=$po 
#Main RM --genomeindex=$GENIND --srcdir=$SRCDIR --datadir=$HOME/Projects/biter_demorree_PAS/data/deMorree_2012/directRNAseq --scoremin=L,0.0,-0.6 --otag=L06minscoreNSC --indir=TR_trim5termT --nosoftclipping=T --parseOutput=$po 
#Main RM --genomeindex=$GENIND --srcdir=$SRCDIR --datadir=$HOME/Projects/biter_demorree_PAS/data/deMorree_2012/directRNAseq --scoremin=L,0.0,-0.6 --otag=L06minscoreNSC --indir=TR_trim5termTanypolyA --nosoftclipping=T --parseOutput=$po 
#Main RM --genomeindex=$GENIND --srcdir=$SRCDIR --datadir=$HOME/Projects/biter_demorree_PAS/data/deMorree_2012/directRNAseq --scoremin=L,0.0,-0.6 --otag=L06minscoreNSC --indir=TR_trim5termTany10A --nosoftclipping=T --parseOutput=$po 
#Main RM --genomeindex=$GENIND --srcdir=$SRCDIR --datadir=$HOME/Projects/biter_demorree_PAS/data/deMorree_2012/directRNAseq --scoremin=L,0.0,-0.6 --otag=L06minscore --parseOutput=$po 
#Main RM --genomeindex=$GENIND --srcdir=$SRCDIR --datadir=$HOME/Projects/biter_demorree_PAS/data/deMorree_2012/directRNAseq --scoremin=L,0.0,-0.6 --otag=L06minscoreNSC --nosoftclipping=T --parseOutput=$po 
#Main RM --genomeindex=$GENIND --srcdir=$SRCDIR --datadir=$HOME/Projects/biter_demorree_PAS/data/deMorree_2012/directRNAseq --scoremin=L,0.0,-0.6 --otag=L06minscoreNSC --indir=TR_gspolyT12 --nosoftclipping=T --parseOutput=$po 
#Main RM --genomeindex=$GENIND --srcdir=$SRCDIR --datadir=$HOME/Projects/biter_demorree_PAS/data/deMorree_2012/directRNAseq --scoremin=L,0.0,-0.6 --otag=L06minscoreNSC --nosoftclipping=T --indir=TR_gspolyT12 --parseOutput=$po 

# 4 Biter
#Main RM --genomeindex=$GENIND --srcdir=$SRCDIR --datadir=$HOME/Projects/biter_biter_stemCellFateRegulation/data/Liu_2015/ATACseq --nc=16 --indir=TR --parseOutput=$po
#Main RM --genomeindex=$GENIND --srcdir=$SRCDIR --datadir=$HOME/Projects/biter_biter_noiseInAging/data/Guess_2015/RNAseq --spliced=T
#Main RM --genomeindex=$GENIND --srcdir=$SRCDIR --datadir=$HOME/Projects/biter_biter_EWSR1/data/Bilen_2015/ChIPseq --nc=16 --genome=hg38 --parseOutput=$po
#Main RM --genomeindex=$GENIND --srcdir=$SRCDIR --datadir=$HOME/Projects/biter_biter_stemCellFateRegulation/data/DellOrso_2016_CellReports/ATACseq --nc=16 --parseOutput=$po
#Main RM --genomeindex=$GENIND --srcdir=$SRCDIR --datadir=$HOME/Projects/biter_biter_stemCellFateRegulation/data/DellOrso_2016_CellReports/ATACseq --nc=16 --indir=TR_ada --parseOutput=$po
#Main RM --genomeindex=$GENIND --srcdir=$SRCDIR --datadir=$HOME/Projects/biter_biter_stemCellFateRegulation/data/Buenrostro_2013_NatureMethods/ATACseq --nc=16 --genome=hg38 --parseOutput=$po
#Main RM --genomeindex=$GENIND --srcdir=$SRCDIR --datadir=$HOME/Projects/biter_biter_stemCellFateRegulation/data/Ikeda_2014/RNAseq --spliced=T --srt=T --nc=16 --parseOutput=$po
#Main RM --genomeindex=$GENIND --srcdir=$SRCDIR --datadir=$HOME/Projects/biter_biter_stemCellFateRegulation/data/Brett_2015/RNAseq --spliced=T --genome=mm10_ERCC92 --srt=T --nc=16 --parseOutput=$po
#Main RM --genomeindex=$GENIND --srcdir=$SRCDIR --datadir=$HOME/Projects/biter_biter_stemCellFateRegulation/data/Liu_2014/RNAseq --spliced=T --srt=T --nc=16 --parseOutput=$po

# 4 Jamie
#Main RM --genomeindex=$GENIND --srcdir=$SRCDIR --datadir=$HOME/Projects/biter_jbrett_DNAmethylation/data/Brett_2014/WGBseq --saveall=T --indir=TR --t=16 --nc=8 --parseOutput=$po 
#Main RM --genomeindex=$GENIND --srcdir=$SRCDIR --datadir=$HOME/Projects/biter_jbrett_DNAmethylation/data/Brett_2014/WGBseq --unique=T --saveall=T --indir=TR --otag=10K --t=0 --tmin=30 --nc=2 --parseOutput=$po --debug=T --skipfastqc=T

#Main RM --genomeindex=$GENIND --srcdir=$SRCDIR --datadir=$HOME/Projects/biter_jbrett_DNAmethylation/data/Brett_2014/WGBseq --t=16 --nc=16 --indir=TR --parseOutput=$po --dryrun=T --skippostproc=T
#Main RM --genomeindex=$GENIND --srcdir=$SRCDIR --datadir=$HOME/Projects/biter_jbrett_DNAmethylation/data/Brett_2015/RNAseq --spliced=T --genome=mm10_ERCC92
#Main RM --genomeindex=$GENIND --srcdir=$SRCDIR --datadir=$HOME/Projects/biter_jbrett_DNAmethylation/data/Ikeda_2014/RNAseq --spliced=T

# 4 Ling - Hairless paper submission
#Main RM --genomeindex=$GENIND --srcdir=$SRCDIR --datadir=$HOME/Projects/biter_lingliu_DNAdamage/data/Liu_2016/ATACseq --nc=16 --indir=TR --parseOutput=$po
#Main RM --genomeindex=$GENIND --srcdir=$SRCDIR --datadir=$HOME/Projects/biter_lingliu_DNAdamage/data/Liu_2016/RNAseq --spliced=T --indir=TR --parseOutput=$po 
#Main RM --genomeindex=$GENIND --srcdir=$SRCDIR --datadir=$HOME/Projects/biter_lingliu_DNAdamage/data/Liu_2016/ChIPseqV2 --nc=16 --parseOutput=$po 
#Main RM --genomeindex=$GENIND --srcdir=$SRCDIR --datadir=$HOME/Projects/biter_lingliu_DNAdamage/data/Liu_2016/ChIPseq2 --nc=16 --parseOutput=$po


# 4 Ling
#Main RM --genomeindex=$GENIND --srcdir=$SRCDIR --datadir=$HOME/Projects/biter_lingliu_DNAdamage/data/Cheung_2013/RNAseq --spliced=T --nc=16
#Main RM --genomeindex=$GENIND --srcdir=$SRCDIR --datadir=$HOME/Projects/biter_lingliu_DNAdamage/data/Liu_2014/RNAseq --spliced=T --nc=16

##Main RM --genomeindex=$GENIND --srcdir=$SRCDIR --datadir=$HOME/Projects/biter_lingliu_DNAdamage/data/Liu_2014/MNaseseq --nc=16
#Main RM --genomeindex=$GENIND --srcdir=$SRCDIR --datadir=$HOME/Projects/biter_lingliu_DNAdamage/data/Liu_2011/ChIPseq --nc=16
#Main RM --genomeindex=$GENIND --srcdir=$SRCDIR --datadir=$HOME/Projects/biter_lingliu_DNAdamage/data/Liu_2012/ChIPseq --nc=16
#Main RM --genomeindex=$GENIND --srcdir=$SRCDIR --datadir=$HOME/Projects/biter_lingliu_DNAdamage/data/Liu_2013/ChIPseq --nc=16
##Main RM --genomeindex=$GENIND --srcdir=$SRCDIR --datadir=$HOME/Projects/biter_lingliu_DNAdamage/data/Liu_2015/ATACseq --nc=16
##Main RM --genomeindex=$GENIND --srcdir=$SRCDIR --datadir=$HOME/Projects/biter_lingliu_DNAdamage/data/Wenemoser_2015/ChIPseq --nc=16
#Main RM --genomeindex=$GENIND --srcdir=$SRCDIR --datadir=$HOME/Projects/biter_lingliu_DNAdamage/data/Liu_2014/MNaseseq --nc=16 --indir=TR_ada
##Main RM --genomeindex=$GENIND --srcdir=$SRCDIR --datadir=$HOME/Projects/biter_lingliu_DNAdamage/data/Liu_2011/ChIPseq --nc=16 --indir=TR_ada
##Main RM --genomeindex=$GENIND --srcdir=$SRCDIR --datadir=$HOME/Projects/biter_lingliu_DNAdamage/data/Liu_2012/ChIPseq --nc=16 --indir=TR_ada
##Main RM --genomeindex=$GENIND --srcdir=$SRCDIR --datadir=$HOME/Projects/biter_lingliu_DNAdamage/data/Liu_2013/ChIPseq --nc=16 --indir=TR_ada
#Main RM --genomeindex=$GENIND --srcdir=$SRCDIR --datadir=$HOME/Projects/biter_lingliu_DNAdamage/data/Liu_2015/ChIPseq --nc=16 --indir=TR_ada
#Main RM --genomeindex=$GENIND --srcdir=$SRCDIR --datadir=$HOME/Projects/biter_lingliu_DNAdamage/data/Wenemoser_2015/ChIPseq --nc=16 --indir=TR_ada

#Main RM --genomeindex=$GENIND --srcdir=$SRCDIR --datadir=$HOME/Projects/biter_lingliu_DNAdamage/data/Liu_2015/RNAseq --spliced=T --srt=T --indir=TR_ada --parseOutput=$po
#Main RM --genomeindex=$GENIND --srcdir=$SRCDIR --datadir=$HOME/Projects/biter_lingliu_DNAdamage/data/Liu_2016/RNAseq --spliced=T --srt=T --indir=TR_ada --parseOutput=$po 
#Main RM --genomeindex=$GENIND --srcdir=$SRCDIR --datadir=$HOME/Projects/biter_lingliu_DNAdamage/data/Liu_2015/RNAseq --spliced=T --srt=T --parseOutput=$po
#Main RM --genomeindex=$GENIND --srcdir=$SRCDIR --datadir=$HOME/Projects/biter_lingliu_DNAdamage/data/Liu_2016/RNAseq --spliced=T --srt=T --parseOutput=$po
#Main RM --genomeindex=$GENIND --srcdir=$SRCDIR --datadir=$HOME/Projects/biter_lingliu_DNAdamage/data/Liu_2016/ChIPseq --nc=16 --parseOutput=$po 
##Main RM --genomeindex=$GENIND --srcdir=$SRCDIR --datadir=$HOME/Projects/biter_lingliu_DNAdamage/data/Liu_2015/ATACseq --nc=16 --indir=TR_ada --parseOutput=$po
#Main RM --genomeindex=$GENIND --srcdir=$SRCDIR --datadir=$HOME/Projects/biter_lingliu_DNAdamage/data/Liu_2015/ChIPseq --nc=16 --saveallm=T --parseOutput=$po
#Main RM --genomeindex=$GENIND --srcdir=$SRCDIR --datadir=$HOME/Projects/biter_lingliu_DNAdamage/data/Zhang_2015_Science/ChIPseq --nc=16 --indir=TR_ada --saveallm=T --genome=hg38 --parseOutput=$po
#Main RM --genomeindex=$GENIND --srcdir=$SRCDIR --datadir=$HOME/Projects/biter_lingliu_DNAdamage/data/Liu_2015/ChIPseq2 --indir=TR_ada --saveallm=T --parseOutput=$po
#Main RM --genomeindex=$GENIND --srcdir=$SRCDIR --datadir=$HOME/Projects/biter_lingliu_DNAdamage/data/Liu_2016/ChIPseq2 --nc=16 --parseOutput=$po 
#Main RM --genomeindex=$GENIND --srcdir=$SRCDIR --datadir=$HOME/Projects/biter_lingliu_DNAdamage/data/Liu_2016/ChIPseqV2 --nc=16 --parseOutput=$po 

# 4 Cindy
## wo/adaptor trimming
#Main RM --genomeindex=$GENIND --srcdir=$SRCDIR --datadir=$HOME/Projects/biter_cvanvelt_SCTrxTrl/data/vanVelthoven_2015/PullDownTrx --spliced=T --genome=mm10_ERCC92
#Main RM --genomeindex=$GENIND --srcdir=$SRCDIR --datadir=$HOME/Projects/biter_cvanvelt_SCTrxTrl/data/vanVelthoven_2015/PullDownTrl --spliced=T --genome=mm10_ERCC92
#Main RM --genomeindex=$GENIND --srcdir=$SRCDIR --datadir=$HOME/Projects/biter_cvanvelt_SCTrxTrl/data/Gonzalez_2014_JNeurosci/PullDownTrl --spliced=T --genome=mm10_ERCC92
## w/adaptor trimming
#Main RM --genomeindex=$GENIND --srcdir=$SRCDIR --datadir=$HOME/Projects/biter_cvanvelt_SCTrxTrl/data/vanVelthoven_2015/PullDownTrx --spliced=T --genome=mm10_ERCC92 --indir=TR_ada
#Main RM --genomeindex=$GENIND --srcdir=$SRCDIR --datadir=$HOME/Projects/biter_cvanvelt_SCTrxTrl/data/vanVelthoven_2015/PullDownTrl --spliced=T --genome=mm10_ERCC92 --indir=TR_ada
#Main RM --genomeindex=$GENIND --srcdir=$SRCDIR --datadir=$HOME/Projects/biter_cvanvelt_SCTrxTrl/data/Gonzalez_2014_JNeurosci/PullDownTrl --spliced=T --genome=mm10_ERCC92 --indir=TR_ada
#Main RM --genomeindex=$GENIND --srcdir=$SRCDIR --datadir=$HOME/Projects/biter_cvanvelt_SCTrxTrl/data/vanVelthoven_2015/PullDownTrx2 --nc=16 --spliced=T --genome=mm10_ERCC92 --indir=TR_ada --srt=T --parseOutput=$po
#Main RM --genomeindex=$GENIND --srcdir=$SRCDIR --datadir=$HOME/Projects/biter_cvanvelt_SCTrxTrl/data/vanVelthoven_2016/PullDownTrx --nc=16 --spliced=T --genome=mm10_ERCC92 --indir=TR_ada --srt=T --parseOutput=$po

# 4 Mike
#Main RM --genomeindex=$GENIND --srcdir=$SRCDIR --datadir=$HOME/Projects/biter_wosczyna_miRNAsInMuscleCellFate/data/Wosczyna_2015/RNAseq --spliced=T
#Main RM --genomeindex=$GENIND --srcdir=$SRCDIR --datadir=$HOME/Projects/biter_wosczyna_miRNAsInMuscleCellFate/data/Wosczyna_2015/PullDown --spliced=T --genome=mm10_ERCC92
#Main RM --genomeindex=$GENIND --srcdir=$SRCDIR --datadir=$HOME/Projects/biter_wosczyna_miRNAsInMuscleCellFate/data/Wosczyna_2015/RNAseq --spliced=T --indir=TR
#Main RM --genomeindex=$GENIND --srcdir=$SRCDIR --datadir=$HOME/Projects/biter_wosczyna_miRNAsInMuscleCellFate/data/Wosczyna_2015/PullDown --spliced=T --genome=mm10_ERCC92 --indir=TR

#Main RM --genomeindex=$GENIND --srcdir=$SRCDIR --datadir=$HOME/Projects/biter_wosczyna_miRNAsInMuscleCellFate/data/Wosczyna_2015/RNAseq --spliced=T --indir=TR_ada_A100_T100
#Main RM --genomeindex=$GENIND --srcdir=$SRCDIR --datadir=$HOME/Projects/biter_wosczyna_miRNAsInMuscleCellFate/data/Wosczyna_2015/PullDown --spliced=T --genome=mm10_ERCC92 --indir=TR_ada_A100_T100


