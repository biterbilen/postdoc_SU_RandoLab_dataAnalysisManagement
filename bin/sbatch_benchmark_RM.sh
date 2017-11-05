#!/bin/bash - 
#===============================================================================
#
#          FILE: sbatch_benchmark_HISAT.sh
# 
#         USAGE: ./sbatch_benchmark_HISAT.sh 
# 
#   DESCRIPTION: 
# 
#       OPTIONS: ---
#  REQUIREMENTS: ---
#          BUGS: ---
#         NOTES: ---
#        AUTHOR: Dr. Biter Bilen (bb), biterbilen@yahoo.com
#  ORGANIZATION: Stanford, Palo Alto, USA
#       CREATED: 01/15/2016 15:59
#      REVISION:  ---
#===============================================================================

set -o nounset                              # Treat unset variables as an error

#SE-----------------------------------------------------------------------------------------------------
# single end
t=3
nc=4 #16
indir=~/Projects/biter_lingliu_DNAdamage/data/Liu_2015/ChIPseq/TR_ada
file=HrCKO_young_NA_SC_NA_NA_1.fastq.gz
indir=~/Projects/biter_lingliu_DNAdamage/data/Liu_2016/ChIPseq/RM/
file=HrWT_young_TMX_SC_NA_H3K9me2_1_un.fastq.gz
indir=~/Projects/biter_lingliu_DNAdamage/data/Liu_2015/ChIPseq/QA
file=WT_old_60hPI_SC_NA_H3K27me3_1.fastq.gz
indir=~/Projects/biter_biter_shared/results/2016-01-15/WT_old_60hPI_SC_NA_H3K27me3_1.w_hisat2
file=SEshort_woT_un.fastq.gz
tag=`basename $file .fastq.gz`

#nc=16; tag2=default; tool=STAR; commonParams="$tool --runThreadN $nc --genomeDir ~/PI_HOME/Data/casco/STAR_indexes/mm10 --readFilesIn $indir/$file --readFilesCommand zless --alignIntronMax 1 --outSAMtype BAM SortedByCoordinate --outReadsUnmapped fastq --outFileNamePrefix"
#tag2=local; tool=hisat; commonParams="$tool -t -p $nc --sensitive-local -x ~/PI_HOME/Data/casco/HISAT_indexes/mm10/chr --no-spliced-alignment -U $indir/$file -S" 
#tag2=endtoend; tool=hisat; commonParams="$tool -t -p $nc --end-to-end -x ~/PI_HOME/Data/casco/HISAT_indexes/mm10/chr --no-spliced-alignment -U $indir/$file -S" 
#tag2=local; tool=bowtie2; commonParams="$tool --sensitive-local --omit-sec-seq -p $nc -x ~/PI_HOME/Data/casco/bowtie2_indexes/mm10/chr  -U $indir/$file -S "
#tag2=endtoend; tool=bowtie2; commonParams="$tool --omit-sec-seq -p $nc -x ~/PI_HOME/Data/casco/bowtie2_indexes/mm10/chr  -U $indir/$file -S "
#tag2=k100; tool=hisat2; commonParams="$tool --omit-sec-seq -k 100 -5 0 -3 0 -t -p $nc -x ~/PI_HOME/Data/casco/HISAT2_indexes/mm10/chr --no-spliced-alignment -U $indir/$file -S "  #finds 46381  
#tag2=scoreminCm100; tool=hisat2; commonParams="$tool --omit-sec-seq --score-min=C,-100,0 -5 0 -3 0 -t -p $nc -x ~/PI_HOME/Data/casco/HISAT2_indexes/mm10/chr --no-spliced-alignment -U $indir/$file -S "
tag2=scoreminLm06m06; tool=hisat2; commonParams="$tool --omit-sec-seq --score-min=L,-0.6,-0.06 -5 0 -3 0 -t -p $nc -x ~/PI_HOME/Data/casco/HISAT2_indexes/mm10/chr --no-spliced-alignment -U $indir/$file -S "

declare -A qcmds
qcmds+=( [$tag\_$tool\_$tag2]="$commonParams" )
#qcmds+=( [SEshort_5_5_3_0T]="$commonParams	-5 5 -3 0" )
#qcmds+=( [SEshort_5_5_3_0T_qcfilter]="$commonParams	-5 5 -3 0 --qc-filter" )


#	${qcmds[$cmd]} $cmd.sam --un-gz $cmd\_un.fastq.gz;
for cmd in ${!qcmds[@]}; do
	mkdir -p $cmd; pushd $cmd;
	sbatch -J $cmd -p normal -t $t:0:0 -c $nc -e slurm_$cmd.out -o slurm_$cmd.out <<- END
	#!/bin/bash -l
	
	date;
	if [ $tool == STAR ]; then
		${qcmds[$cmd]} $cmd;
	else
		${qcmds[$cmd]} $cmd.sam --un-gz $cmd\_un.fastq.gz;
		if [ $tool == bowtie2 ]; then
			grep -v "XS:i:" $cmd.sam | samtools view -@ $nc -b -o $cmd.bam -;
		else
			grep -P "^@|NH:i:1$|NH:i:1\t" $cmd.sam | samtools view -@ $nc -b -o $cmd.bam -;
		fi;
		samtools view $cmd.bam | wc -l;
		#rm $cmd.sam
	fi;

	date;

	END
	popd
done
exit




#PE-----------------------------------------------------------------------------------------------------
# paired end
t=3
nc=4 #16
indir=~/Projects/biter_lingliu_DNAdamage/data/Liu_2015/ATACseq/TR_ada;
file1=$indir/WT_old_NA_SC_1_1.fastq.gz
file2=$indir/WT_old_NA_SC_1_2.fastq.gz
name=WT_old_NA_SC_1
tag2=default; tool=hisat2; cmd="$tool --omit-sec-seq -5 0 -3 0 -t -p $nc -x ~/PI_HOME/Data/casco/HISAT2_indexes/mm10/chr --no-spliced-alignment -1 $file1 -2 $file2 --un-conc-gz $name\_un.fastq.gz -S $name.sam"  #finds none
tag2=X2000; tool=hisat2; cmd="$tool -X 2000 --omit-sec-seq -5 0 -3 0 -t -p $nc -x ~/PI_HOME/Data/casco/HISAT2_indexes/mm10/chr --no-spliced-alignment -1 $file1 -2 $file2 --un-conc-gz $name\_un.fastq.gz -S $name.sam"  #finds none
tag2=scoreminLm06m06; tool=hisat2; cmd="$tool --score-min L,-0.6,-0.6 --omit-sec-seq -5 0 -3 0 -t -p $nc -x ~/PI_HOME/Data/casco/HISAT2_indexes/mm10/chr --no-spliced-alignment -1 $file1 -2 $file2 --un-conc-gz $name\_un.fastq.gz -S $name.sam"  #finds none
tag2=defaultNamelyscoreminLm06m06; tool=bowtie2; cmd="$tool --omit-sec-seq -t -p $nc -x ~/PI_HOME/Data/casco/bowtie2_indexes/mm10/chr -1 $file1 -2 $file2 --un-conc-gz $name\_un.fastq.gz -S $name.sam"  #finds none

#tag2=endtoend; tool=bowtie2; commonParams="$tool --omit-sec-seq -p $nc -x ~/PI_HOME/Data/casco/bowtie2_indexes/mm10/chr  -U $indir/$file -S "

odir=$name\_w_$tool\_$tag2
mkdir -p $odir; pushd $odir
ln -sf $file1 $file2 .
sbatch -J $tool -p normal -t $t:0:0 -c $nc -e slurm_$name.out -o slurm_$name.out <<- END
	#!/bin/bash -l
	$cmd
END
popd
exit


