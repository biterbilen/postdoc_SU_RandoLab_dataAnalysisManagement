#!/bin/bash - 
#===============================================================================
#
#          FILE: sbatchLike_getDataMetadataFromENCODE.sh
# 
#         USAGE: ./sbatchLike_getDataMetadataFromENCODE.sh
# 
#   DESCRIPTION: 
# 
#       OPTIONS: ---
#  REQUIREMENTS: ---
#          BUGS: ---
#         NOTES: ---
#        AUTHOR: Dr. Biter Bilen (bb), biterbilen@yahoo.com
#  ORGANIZATION: Stanford, Palo Alto, USA
#       CREATED: 12/05/2016 16:47
#      REVISION:  ---
#===============================================================================

set -o nounset                              # Treat unset variables as an error
source ~/Projects/biter_biter_shared/bin/auxslurm.sh

declare -A o=( [jobnames]=getENCODEData [queue]=normal [t]=1 [nc]=1 [infle]='Bernstein_2014/urls.txt' [outdir]=Bernstein_2014 [grep]=H3K9me2 )
declare -A o=( [jobnames]=getENCODEMetadata [queue]=normal [t]=1 [nc]=1 [query]="target='H3',file_format='fastq'" [outdir]=Bernstein_2014 [grep]=H3K9me2 )
declare -A o=( [jobnames]=getENCODEMetadata [queue]=normal [t]=1 [nc]=1 [query]="biosample='C2C12',file_format='bed'" [outdir]=Wold_2011 )
declare -A o=( [jobnames]=getENCODEData [queue]=normal [t]=1 [nc]=1 [infle]='Wold_2011/urls.txt' [outdir]=Wold_2011)
declare -A o=( [jobnames]=getENCODEMetadata [queue]=normal [t]=1 [nc]=1 [query]="biosample='GM12878',file_format='bed',assay='ChIP-seq'" [outdir]=GM12878 )
declare -A o=( [jobnames]=getENCODEData [queue]=normal [t]=1 [nc]=1 [infle]='GM12878/urls.txt' [outdir]=GM12878 [over]=~/PI_HOME/Data/casco/UCSC_tracks/hg38/hg19ToHg38.over.chain.gz )

declare -A qcmds=()

#	curl -O -L www.encodeproject.org/{URLT} \
#urlts=(\$(less metadata.target\:H3_file_format\:fastq | grep H3K9me2 | cut -f 3 | while read p; do grep \$p metadata.target\:H3_file_format\:fastq; done  | cut -f 1 | sort | uniq))
# 1. TODO incomplete for generalization
# gets H3K9me2 associated bio-samples, which will include controls and treatments of biosamples
jobname=getENCODEMetadata; qcmd=$(cat <<- END
Rscript \
	-e "source('~/Projects/biter_biter_shared/bin/auxENCODExplorer.R')" \
	-e "caller.getENCODETreatControl(${o[query]:-file_format='fastq'}, nc=${o[nc]})"
f=\$(ls metadata\.*:*)
head -n 1 \$f > metadata
# URLS
less \$f | grep -v href | \
	if [ ${o[grep]:-X} != X ]; then
		grep H3K9me2 
	else
		less
	fi | cut -f 5 | \
	while read p; do grep \$p \$f; done | \
	cut -f 3 | sort | uniq | \
  while read u; do echo www.encodeproject.org\$u; done > urls.txt 
# metadata
xargs -n 1 basename < urls.txt | \
	while read u; do grep \$u \$f; done >> metadata 
cp metadata metadata.back
END
)
if [[ ,${o[jobnames]}, =~ ,$jobname, ]]; then qcmds+=( [$jobname]="$qcmd" ); fi

jobname=getENCODEData; 
while read u; do
	b=$( echo ${u##*/} )     # basename w suffix
	#suffix=$( echo ${b#*.} ) # suffix # TODO generalize
	b=$( echo ${b%%.*} )     # basename wo suffix
	qcmd=$(cat <<- END
	curl -O -L $u
	if [ ${o[over]:-NA} != NA ]; then
		# TODO generalize
		# TODO narrow.bed
		less $b.bed.gz | awk 'BEGIN{ OFS="\t";} { print \$1,\$2,\$3,\$4,\$7,\$6 }' > $b\_DEL.bed

		liftOver $b\_DEL.bed $(auxDirNameTune ${o[over]}) $b.bed $b.unmapped

		mv $b.bed.gz $b.bed.gz.orginal
		gzip $b.bed 

		if [ ! -s $b.unmapped ]; then rm $b.unmapped; fi
		rm $b\_DEL.bed
	fi
	END
	)
	if [[ ,${o[jobnames]}, =~ ,$jobname, ]]; then qcmds+=( [$jobname$b]="$qcmd" ); fi
done < ${o[infle]:-urls.txt}


mkdir -p ${o[outdir]:-outdir}; pushd ${o[outdir]:-outdir};
mkdir -p bin
# TODO Integrate
#slurmSubmitToQueue

#for jobname in $(echo ${o[jobnames]} | sed 's/,/ /g'); do
for jobname in ${!qcmds[@]}; do
	echo -e "Doing $jobname"
	qcmd="${qcmds[$jobname]}"
	wqcmd="sbatch -J $jobname -p ${o[queue]} -t ${o[t]}:0:0 -c ${o[nc]} \
		-e slurm_$jobname.out -o slurm_$jobname.out bin/sbatch_$jobname.sh"
	cat > bin/sbatch_$jobname.sh <<- END
	#!/bin/bash -l 
	# sbatch command: $wqcmd
	$qcmd
	echo DONE
	END
	$wqcmd
done

popd

