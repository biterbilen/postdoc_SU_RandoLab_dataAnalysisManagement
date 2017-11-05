#!/bin/bash - 
#===============================================================================
#
#          FILE: cmds.sh
# 
#         USAGE: ./cmds.sh 
# 
#   DESCRIPTION: 
# 
#       OPTIONS: ---
#  REQUIREMENTS: ---
#          BUGS: ---
#         NOTES: ---
#        AUTHOR: Dr. Biter Bilen (bb), biterbilen@yahoo.com
#  ORGANIZATION: Stanford, Palo Alto, USA
#       CREATED: 01/20/2016 17:01
#      REVISION:  ---
#===============================================================================

set -o nounset                              # Treat unset variables as an error
nc=8; t=8
nc=16; t=8

odir=LM
mkdir -p $odir; pushd $odir

# test
regDir=~/PI_HOME/Data/casco/UCSC_tracks/mm10
reg=a

regDir=~/PI_HOME/Data/casco/UCSC_tracks/mm10/canonical
reg=promoters

fle=$regDir/$reg.bed.gz; regions=self
# prep TFBS table
motifSrc=JASPAR2014; motifSrcName=JASPAR2014;                          matrixtype="PWM"; all_versions=T;
motifSrc=~/PI_HOME/Data/casco/Homer/custom.motifs; motifSrcName=Homer; matrixtype="PWM"; ii=2; byrow=T;
motifSrc=~/PI_HOME/Data/casco/JASPAR/JASPAR_CORE_nonredundant_pfm_vertebrates.txt; motifSrcName=JASPARcoreVert; matrixtype="PFM"; ii=1; byrow=F;
tag=$motifSrcName\_$reg
sbatch -J $tag -p normal -t $t:0:0 -c $nc -e slurm_$tag.out -o slurm_$tag.out <<- END
#!/bin/bash -l
Rscript -e 'source("~/Projects/biter_biter_shared/bin/auxML.R")' \
	-e 'caller.getMotifCountsInRegion("$fle", as.fasta=F, expressedTx.colInd=1, regions="$regions",  src="$motifSrc", matrixtype="$matrixtype", all_versions=$all_versions, ii=$ii, byrow=$byrow,   ofleTag="$tag", mc.cores=$nc)';
Rscript --version;
END

exit


# prep DM table 
#TODO
samplenamepat=(A4 Q4 R5)
samplenamepat=(Q4 Q10 Q16 Q22)

tag=$reg$motifSrc
sbatch -J $tag -p normal -t $t:0:0 -c $nc -e slurm_$tag.out -o slurm_$tag.out <<- END
#!/bin/bash -l
Rscript -e 'source("~/Projects/biter_biter_shared/bin/auxMethylSeekR.R")' \
	-e 'dta <- caller.estimate.combine.callDMR("$fle", "$msrFilesStr","$msrFilesnamesStr", regionfile.wheader=$regionfilewheader, nCpG.cutoff=3, min.CpGcover=4, coveredCG.ratio=0.25, cols2save="$cols2save", samplename.pat="$samplenamepatStr", b samplename.pat =$nc, outfle.tag="$outfleTag",save.image=$saveimage, debug=$debug)'
Rscript --version;
END


popd $odir
