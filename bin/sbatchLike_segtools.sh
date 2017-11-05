#!/bin/bash - 
#===============================================================================
#
#          FILE: sbatchLike_w_segtools.sh
# 
#         USAGE: ./sbatchLike_w_segtools.sh 
# 
#   DESCRIPTION: 
# 
#       OPTIONS: ---
#  REQUIREMENTS: ---
#          BUGS: ---
#         NOTES: ---
#        AUTHOR: Dr. Biter Bilen (bb), biterbilen@yahoo.com
#  ORGANIZATION: Stanford, Palo Alto, USA
#       CREATED: 04/20/2016 10:21
#      REVISION:  ---
#===============================================================================

set -o nounset                              # Treat unset variables as an error

t=1
nc=8; nc=8

jobname=segtools
ANNOTATION="$HOME/PI_HOME/Data/casco/UCSC_tracks/mm10/canonical/transcripts.bed.gz"
ANNOTATIONGTF="/home/biter/PI_HOME/Data/casco/UCSC_tracks/mm10/mm10_refGene_canonical.gtf.gz"

#jobname=genomedata

sbatch -J $jobname -p normal -t $t:0:0 -c $nc -e slurm_$jobname.out -o slurm_$jobname.out <<- BATCH
#!/bin/bash -l

# for jobname=genomedata
#genomedata-load -s ~/PI_HOME/Data/casco/UCSC_tracks/mm10/assembly/all/chr.fa genomedata.mm10

wget http://pmgenomics.ca/hoffmanlab/proj/segway/2011/test.genomedata
genomedata-info tracknames_continuous test.genomedata
genomedata-info contigs test.genomedata
segway --num-labels=4 train test.genomedata traindir
segway identify test.genomedata traindir identifydir


echo DONE
exit

echo $jobname
#cat PF_Liu_2016_ChIPseq/bdgdiff/diff_HrCKO_young_TMX_SC_NA_H3K9me2_vs_HrWT_young_TMX_SC_NA_H3K9me2_c3.0_co* | grep "^chr" | \
	#chrsed 's/diff_HrCKO_young_TMX_SC_NA_H3K9me[23]_vs_HrWT_young_TMX_SC_NA_//g' | sed 's/_[1-9][0-9]*//g' | \
	#sedsort -k 1,1 -k 2,2g | bedtools cluster -i stdin | \
	#stdinuniq -f 5 -c | awk 'BEGIN{OFS="\t"}{ if(\$1==1) print \$2,\$3,\$4,\$5,\$6}' > tmp.sorted.bed
##tmp.sorted.bedsort -k 4,4 > tmp.sorted.bed
##tmp.sorted.bedsort -k 1,1 -k 2,2g > tmp.sorted.bed

## Most likely requires ungapped segments like in a wig file?
##segtools-transition --clobber tmp.sorted.bed 

segtools-preprocess --clobber tmp.sorted.bed 
# ? segtools-aggregation --clobber --mode=region --group --normalize tmp.sorted.bed.pkl.gz $ANNOTATIONGTF

segtools-length-distribution --clobber tmp.sorted.bed.pkl.gz 
segtools-nucleotide-frequency --clobber tmp.sorted.bed.pkl.gz genomedata.mm10
segtools-feature-distance --clobber -s -p -a tmp.sorted.bed.pkl.gz $ANNOTATION > feature_distance/feature_distance.all.tab
segtools-overlap --clobber -b segments tmp.sorted.bed.pkl.gz $ANNOTATIONGTF

segtools-html-report --clobber tmp.sorted.bed.pkl.gz

echo DONE
BATCH

