#!/bin/bash - 
#===============================================================================
#
#          FILE: sbatch_getfastqDump.sh
# 
#         USAGE: ./sbatch_getfastqDump.sh 
# 
#   DESCRIPTION: 
# 
#       OPTIONS: ---
#  REQUIREMENTS: ---
#          BUGS: ---
#         NOTES: ---
#        AUTHOR: Dr. Biter Bilen (bb), biterbilen@yahoo.com
#  ORGANIZATION: Stanford, Palo Alto, USA
#       CREATED: 02/08/2015 20:07
#      REVISION:  ---
#===============================================================================

# TODO Prep: Download sra_result.csv like
# http://www.sthda.com/english/wiki/wiki.php?id_contents=7273
# 


# Alternative to wget download of sra files
#echo Solve ncbi-cache issue or download SRR files first
#cache-mgr --clear
#prefetch $i $i.sra did not work; some remnants from previous direct download? / Repository/data related?

set -o nounset                              # Treat unset variables as an error
# ENA runs
tag=BS-Seq
tag=J1
urlInfo=URLs.txt; urli=1
SraRunInfo=NULL
sra_result=NULL
FS=" "

# GEO runs
tag=RNA-Seqx
urlInfo=SraRunInfo.csv; urli=12;
SraRunInfo=SraRunInfo.csv; lli=19;
FS=","
#urlInfo=sra_result

# GEO runs
# TODO set tag, tag2, and go 
wgetComplete=${1:-F}
prjn=PRJNA223060; tag=OTHER; tag2="Mus musculus"; lli=18; urli=12; srxi=13; go="-Pw" #Gonzalez_2014
prjn=PRJNA205299; tag=ChIP-Seq; tag2=""; go="-Pw" #lli=19; urli=12; srxi=13; #Liu_2014
prjn=SRP066889; tag=ChIP-Seq; tag2=""; go="-Pw" #EWSR1
prjn=SRP004081; tag=ChIP-Seq; tag2=""; go="-Pw" #Soleimani PAX3 and PAX7
prjn=SRP032942; tag=ChIP-Seq; tag2=""; go="-Pw" # Zhang 2015 H3K9me3
prjn=SRP067391; tag=OTHER; tag2=""; go="-Pw" # DellOrso 2016 CellReports ATAC-seq
prjn=PRJNA207663; tag=OTHER; tag2=""; go="-Pw" # Buenrostro_2013_NatureMethods ATAC-seq
prjn=PRJNA162295; tag="Orc1_ChIPseq|High-density|input DNA|low-density fractions"; tag2=""; go="-Pw" # Dellino_2013_GenomeResearch GSE37583 
prjn=SRP012231; tag="36bp,"; tag2=""; go="-P" # Wold_2011 GSE36024
urlInfo=SraRunInfo.csv; 
SraRunInfo=SraRunInfo.csv; wget -q -O $SraRunInfo "http://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?save=efetch&db=sra&rettype=runinfo&term=$prjn"
sra_result=sra_result.csv; # TODO find rettype for sra_result.csv file; ask in forum because cra.cgi code didn't give much clue; wget -O $SraRunInfo "http://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?save=efetch&db=sra&term=SRP066889&rettype=sra_access"
FS=","

#source ~/Projects/biter_biter_shared/bin/auxslurm.sh
#urli=$(auxGetIndexFromColName $urlInfo download_path $FS)
#lli=$(auxGetIndexFromColName $urlInfo LibraryLayout $FS)
#srxi=$(auxGetIndexFromColName $urlInfo Experiment $FS)
urli=$(head -n 1 $urlInfo | sed "s/$FS/\n/g"  | awk '$1=="download_path"{ print NR}')
lli=$(head -n 1 $urlInfo | sed "s/$FS/\n/g"  | awk '$1=="LibraryLayout"{ print NR}')
srxi=$(head -n 1 $urlInfo | sed "s/$FS/\n/g"  | awk '$1=="Experiment"{ print NR}')

echo $urli $lli $srxi

if [ $wgetComplete == F ]; then
	#less $SraRunInfo | grep -wP "$tag" | grep -wP "$tag2" | \
	less $SraRunInfo | grep $go "$tag" | grep $go "$tag2" | \
		awk -v F="$FS" 'BEGIN{FS=F;}{ print $1}' | \
		while read i; do 
			sraURL=`less $urlInfo | grep -w $i | awk -F "," '{ print $j}' j=$urli`
			libraryLayout=`less $SraRunInfo | grep -w $i | awk -v F="$FS" 'BEGIN{ FS=F;} { print $j}' j=$lli`
			sbatch -J fqDump_$i -p normal -t 16:0:0 -c 4 <<- END
			#!/bin/bash -l
			date 

			wget -nv -c $sraURL -O $i.sra

			if [ $libraryLayout == SINGLE ]; then
				fastq-dump --gzip $i.sra
			else 
				fastq-dump --gzip --split-files $i.sra
			fi

			rm $i.sra

			echo DONE
			date
			END
			echo $i; 
		done 
else
	echo Alias Alias2 Replicate Remove LibType Pair Cell Age Genotype Tissue Molecule Strain Gender Owner Antibody Date invivoGrowthProt invivoTreatProt exVivoGrowthProt exVivoTreatProt ExtrProt LibPrepProt SeqInst | sed 's/ /\t/g' > metadata
	echo Alias Alias2 Replicate Remove LibType Pair SeqInst | sed 's/ /\t/g' >> metadata
	#ls *fastq.gz | while read i; do 
	less $SraRunInfo | grep $go "$tag" | grep $go "$tag2" | \
		bioawk -v j=$srxi -F "$FS" '{ print $j}' | sort | uniq | \
		while read srx; do 
			lli=$(head -n 1 $SraRunInfo | sed "s/$FS/\n/g"  | awk '$1=="LibraryLayout"{ print NR}')
			libraryLayout=$(grep -w $srx $SraRunInfo | awk -v F="$FS" 'BEGIN{ FS=F;} { print $j}' j=$lli)
			sii=$(head -n 1 $sra_result | sed "s/$FS/\n/g"  | awk '$1=="Instrument"{ print NR}')
			seqInst=$(grep -w $srx $sra_result | awk -v F="$FS" 'BEGIN{ FS=F;} { print $j}' j=$sii | \
				sed 's/ //g' | sed 's/"//g')
			a2i=$(head -n 1 $sra_result | sed "s/$FS/\n/g"  | awk '$1=="Experiment Title"{ print NR}')
			alias2=$(grep -w $srx $sra_result | awk -v j=$a2i -F "$FS" '{ print $j}' | \
				awk -F ";" '{ print $1 }' | cut -d " " -f 2- | sed 's/ /_/g' | sed 's/"//g') 
			#srrs=$(grep -w $srx $SraRunInfo | awk -F "$FS" '{ print $1".fastq.gz"}');
			srrs=($(grep -w $srx $SraRunInfo | awk -F "$FS" '{ print $1".fastq.gz"}')); # bug corrected 20160416 - might effect sample sets with multiple srr files per one srx
			if [ $libraryLayout == PAIRED ]; then
				cat $(echo ${srrs[@]} | sed 's/.fastq.gz/_1.fastq.gz/') > $srx\_1.fastq.gz
				cat $(echo ${srrs[@]} | sed 's/.fastq.gz/_2.fastq.gz/') > $srx\_2.fastq.gz
				srrs=($(echo ${srrs[@]} | sed 's/.fastq.gz/_1.fastq.gz/') $(echo ${srrs[@]} | sed 's/.fastq.gz/_2.fastq.gz/')) #added 20160627
				echo $srx | awk '{ OFS="\t"; print $1"_"1, al, 0, 0, "Xunstranded", 1, si }' al=$alias2 si=$seqInst 
				echo $srx | awk '{ OFS="\t"; print $1"_"2, al, 0, 0, "Xunstranded", 2, si }' al=$alias2 si=$seqInst
			else 
				cat ${srrs[@]} > $srx.fastq.gz
				echo $srx | awk '{ OFS="\t"; print $1, al, 1, 0, "Xunstranded", "NA", si }' al=$alias2 si=$seqInst
			fi
			rm ${srrs[@]}
		done | sort -k 2,2 | uniq >> metadata
fi

if [ $SraRunInfo == NULL ]; then
	echo Library Alias Group Remove LibraryLayout LibraryType | sed 's/ /\t/g' > sample.info
	less ena.xml  | grep -P "<TITLE>|<ID>h" | \
		awk 'BEGIN{ RS="</ID>\n"; FS="\n"; OFS=";"; } { print $1,$2 }' | sed 's/<\/*[A-Z]*>//g' | \
		awk 'BEGIN{ FS="/"; OFS=";" }{ print $1, $NF}' | \
		awk -F ";" 'BEGIN{OFS="\t";}{ lo="SINGLE"; if ($1~ "paired") lo="PAIRED"; split($2,a, " "); s=a[1]"_"a[3]; print $NF,s,a[1],0,lo,a[2] }' | \
		sed 's/_experiment//g' >> sample.info
	SraRunInfo=sample.info; lli=5
fi



exit
