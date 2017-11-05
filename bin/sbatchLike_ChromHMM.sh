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
#       CREATED: 04/26/2016 15:37
#      REVISION:  ---
#===============================================================================

set -o nounset                              # Treat unset variables as an error

# annotate
# TODO
# check w anova expression variance
less GB_PF_RM_TR_ada_Liu_2015_ATACseq_ChromHMM/filterdup_dense.bed | grep "^chr" | cut -f 1-6 | sort -k 1,1 -k 2,2g | bedtools closest -s -a stdin -b ~/PI_HOME/Data/casco/UCSC_tracks/mm10/canonical/transcripts.bed.gz

exit

clf=~/aux/genomes/mouse.mm10.genome
cmft=meta
ibd=filterdup

# metadata
ls $ibd/*bed.gz | awk -F _ '{ OFS="\t"; print $4, $2"_"$3, $0}' |\
	sed "s/$ibd\///g" | sed 's/.gz//g'	> $cmft

# train FIXME did not work on the clusters!
out=$ibd\_training
for f in $ibd/*bed.gz; do 
	less $f | grep ^chr19 > $out/`basename $f .gz`
done
java -mx4000M -jar ./ChromHMM/ChromHMM.jar BinarizeBed $clf $out $cmft $out
# this did not work in the sdev environment
java -mx4000M -jar ./ChromHMM/ChromHMM.jar LearnModel $out $out 8 mm10

# identify
out2=$ibd\_identification
for f in $ibd/*bed.gz; do 
	less $f > $out2/`basename $f .gz`
done
java -mx4000M -jar ./ChromHMM/ChromHMM.jar BinarizeBed $clf $out2 $cmft $out2
java -mx4000M -jar ./ChromHMM/ChromHMM.jar MakeSegmentation $out/model_*.txt $out2 $out2 


# GB
name=DNAaccAyQyOo
out3=GB_PF_RM_TR_ada_Liu_2015_ATACseq_ChromHMM; 
opref=$out3/$ibd
mkdir -p $out3
java -mx4000M -jar ./ChromHMM/ChromHMM.jar MakeBrowserFiles $out2/*_segments.bed $name $opref
bigDataUrl=http://www-dev.stanford.edu/group/rando_lab/biter
cat > $out3/hub.txt <<- END1
hub $name
shortLabel $name
longLabel $name
genomesFile genomes.txt
email biter@stanford.edu 
END1

cat > $out3/genomes.txt <<- END1
trackDb tracks.txt
genome mm10 
END1

# https://genome.ucsc.edu/goldenpath/help/bigBed.html Example Three
wget https://genome.ucsc.edu/goldenpath/help/examples/bedExample2.as
less bedExample2.as | grep -vP "geneSymbol|spID" > bedExample2_4bed9.as
echo browser position chr10:212000-250000 > $out3/tracks.txt
for s in $(find $out3 -name *bed 2> /dev/null); do
	f=$(basename $s .bed)
	less $s | grep ^chr | sort -k 1,1 -k 2,2g > tmp
	bedToBigBed -as=bedExample2_4bed9.as -type=bed9 -extraIndex=name tmp $clf $s.bb
	rm tmp
cat >> $out3/tracks.txt <<- END1

track $f
shortLabel $f
longLabel $f
type bigBed
bigDataUrl $bigDataUrl/$s.bb
visibility full
itemRgb on
nextItemButton on
END1
done

rsync -e 'ssh -p 22' -Lazp $out3 biter@cardinal.stanford.edu:/afs/ir/group/rando_lab/WWW/biter



