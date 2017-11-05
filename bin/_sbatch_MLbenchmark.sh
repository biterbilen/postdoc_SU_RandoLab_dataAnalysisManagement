#!/bin/bash -l 
# sbatch -J stack -p normal -t 4:0:0 -c 8 -e stack.out -o stack.out _sbatch_MLbenchmark.sh
# sbatch -J bag   -p normal -t 4:0:0 -c 8 -e bag.out -o bag.out _sbatch_MLbenchmark.sh
# sbatch -J boost -p normal -t 4:0:0 -c 8 -e boost.out -o boost.out _sbatch_MLbenchmark.sh
date
nseq=10000
nseq=2000
#nseq=500
nc=8
ii=5
Rscript -e "source('~/Projects/biter_biter_shared/bin/auxHighLevelGenomicFeatures.R'); test.caller.getSeq('stack',$nseq,$nc,$ii); sessionInfo();"
#Rscript -e "source('~/Projects/biter_biter_shared/bin/auxHighLevelGenomicFeatures.R'); test.caller.getSeq('bag',$nseq,$nc,$ii); sessionInfo();"
#Rscript -e "source('~/Projects/biter_biter_shared/bin/auxHighLevelGenomicFeatures.R'); test.caller.getSeq('boost',$nseq,$nc,$ii); sessionInfo();"
date
echo DONE
