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

nc=1; t=1
gnm="mouse.mm10"

# caller for TSS and TES region extraction
for tag in canonical; do
#for tag in proximal distal; do
	trxfile=$PI_HOME/Data/casco/UCSC_tracks/mm10/mm10_refGene_$tag.v2.gtf.gz
	mkdir -p $tag; pushd $tag
	reg=transcripts
	oflet=$reg
	Rscript -e 'source("~/Projects/biter_biter_shared/bin/auxGenomicFeatures.R")' \
		-e "caller.expressedTxRegionsOrSeqs('$trxfile', regions='$reg', regions.only=T, ofleTag='', mc.cores=$nc)"; 

	less $reg.bed.gz | bedtools flank -s -l 1 -r 0 -g ~/aux/genomes/$gnm.genome -i stdin | \
		gzip -c > ${reg}_TSS.bed.gz;

	less $reg.bed.gz | bedtools flank -s -l 0 -r 1 -g ~/aux/genomes/$gnm.genome -i stdin | \
		gzip -c > ${reg}_TES.bed.gz;

	popd;
done
Rscript --version;

exit

nc=8; t=8
gnm="mouse.mm10"

# caller for canonical promoter regions extraction
# TODO set
for tag in canonical; do
#for tag in proximal distal; do
	trxfile=../mm10_refGene_$tag.v2.gtf.gz
	#TODO set
	reg=TES # TODO run first reg=transcripts for this option
	reg=promoters
	reg=transcripts
	oflet=$reg
	random=$RANDOM
	sbatch -J $reg -p normal -t $t:0:0 -c $nc -e slurm_${tag}_$reg.out -o slurm_${tag}_$reg.out <<- END
	#!/bin/bash -l
	mkdir -p $tag; pushd $tag
	[[ $reg == TES ]] && 
	ln -sf transcripts.bed.gz $random$reg.bed.gz ||
		Rscript -e 'source("~/Projects/biter_biter_shared/bin/auxGenomicFeatures.R")' \
		-e "caller.expressedTxRegionsOrSeqs('$trxfile', regions='$reg', regions.only=T, ofleTag='$random', mc.cores=$nc)"; 

	zless $random* | cut -f 1-6 | gzip -c > $oflet.bed.gz;
	[[ $reg == transcripts ]] &&
		less $oflet.bed.gz | bedtools slop -s -l -1000 -r 0 -i stdin -g ~/aux/genomes/$gnm.genome | \
			awk '\$3>\$2{ print }' | gzip -c > $oflet\_wopromoters.bed.gz 

	[[ $reg == TES ]] &&
		less $oflet.bed.gz | bedtools flank -s -l 0 -r 1000 -g ~/aux/genomes/$gnm.genome -i stdin | \
			bedtools slop -s -l 2000 -r 0 -i stdin -g ~/aux/genomes/$gnm.genome | \
			awk '\$3>\$2{ print }' | gzip -c > tmp$oflet.bed.gz; mv tmp$oflet.bed.gz $oflet.bed.gz; 

	rm $random*;
	Rscript --version;
	popd;
	END
done

exit

# caller for regions extraction
# TODO set
reg=promoters; oflet=mm10_refGene_$reg
reg=transcripts; oflet=mm10_refGene_$reg
random=$RANDOM
sbatch -J $reg -p normal -t $t:0:0 -c $nc -e slurm_$reg.out -o slurm_$reg.out <<- END
#!/bin/bash -l
Rscript -e 'source("~/Projects/biter_biter_shared/bin/auxGenomicFeatures.R")' \
	-e 'caller.expressedTxRegionsOrSeqs("mm10_refGene.gtf.gz", regions="$reg", regions.only=T, ofleTag="$random", mc.cores=$nc)'; 
Rscript --version;
zless $random* | cut -f 1-6 | gzip -c > $oflet.bed.gz;
if [[ $reg -eq transcripts ]]; then
	less $oflet.bed.gz | bedtools slop -s -l -1000 -r 0 -i stdin -g ~/aux/genomes/$gnm.genome | \
		awk '$3>$2{ print }' | gzip -c > $oflet\_wopromoters.bed.gz
fi
rm $random*;
END


