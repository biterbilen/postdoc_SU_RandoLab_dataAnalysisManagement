#!/bin/bash - 
#===============================================================================
#
#          FILE: _sbatch_benchmark_RM_of_WGBseq.sh
# 
#         USAGE: ./_sbatch_benchmark_RM_of_WGBseq.sh 
# 
#   DESCRIPTION: 
# 
#       OPTIONS: ---
#  REQUIREMENTS: ---
#          BUGS: ---
#         NOTES: ---
#        AUTHOR: Biter Bilen (), biterbilen@yahoo.com
#  ORGANIZATION: 
#       CREATED: 09/16/2016 05:56:36 PM
#      REVISION:  ---
#===============================================================================

set -o nounset                              # Treat unset variables as an error
source ~/Projects/biter_biter_shared/bin/auxslurm.sh
po=${1:-F}

RM() {
	local analysisdirbase=$1
	local genome=${o[genome]:-mm10}

	#local param="${o[param]:-'--sensitive-local'} -p ${o[nc]} -5 ${o[trim5]:-0} -3 ${o[trim3]:-0} \
	#hisat2		
	local param="-p ${o[nc]} -5 ${o[trim5]:-0} -3 ${o[trim3]:-0} \
		-x ${o[genomeindex]:-/home/biter/PI_HOME/Data/casco/HISAT2_indexes/$genome/chr}"

	if [ ${o[scoremin]:-X} != X ]; then
		param="$param --score-min ${o[scoremin]:-X}"
	fi	
	if [ ${o[rdg]:-X} != X ]; then
		param="$param --rdg ${o[rdg]-:X}"
	fi	
	if [ ${o[rfg]:-X} != X ]; then
		param="$param --rfg ${o[rfg]:-X}"
	fi	

	if [[ "${o[spliced]:-F}" == T ]]; then
		# TODO prep in advance
		local allTXfile=${o[allTXfile]:-/home/biter/PI_HOME/Data/casco/UCSC_tracks/$genome/${genome}_refGene.gtf.gz}
		ln -sf $allTXfile .
		less $allTXfile | hisat2_extract_splice_sites.py - > spliceSites.txt
		#local tmp=$(basename $allTXfile .gz)
		#less $allTXfile > $tmp;
		#hisat2_extract_splice_sites.py $tmp > spliceSites.txt;
		#rm $tmp
		param="$param --known-splicesite-infile spliceSites.txt"
	else 
		param="$param --no-spliced-alignment"
	fi

	local pairi=$(auxGetIndexFromColName metadata ${o[pair]:-Pair})
	local extrProti=$(auxGetIndexFromColName metadata ${o[extrProt]:-ExtrProt})
	local lti=$(auxGetIndexFromColName metadata ${o[extrProt]:-LibType})
	local moli=$(auxGetIndexFromColName metadata ${o[molecule]:-Molecule})

	local extrProti=$(auxGetIndexFromColName metadata ${o[extrProt]:-ExtrProt})
	local namei=$(auxGetIndexFromColName metadata ${o[name]:-Name})
	local -a names=($(awk -v namei=$namei 'NR>1{print $namei}' metadata | sort | uniq))
	local genomefafile=${o[genomefafile]:-/home/biter/PI_HOME/Data/casco/HISAT2_indexes/$genome/all/chr.fa}
	local indir=${o[indir]:-QA}; mkdir -p $indir;
	local name
	# TODO revisit single-end mapping after checking log file sanity
	for name in ${names[@]}; do
		ln -sf ${o[datadir]}/$indir/$name*.fastq.gz $indir/

		local allparam="$param" 
		local libType="$(grep -w -m 1 $name metadata | cut -f $lti)"
		if [ $libType != unstranded ]; then allparam="$allparam --rna-strandness $libType"; fi

		if [ ${o[X]:-NA} != NA ]; then allparam="$allparam -X ${o[X]}"; fi

		local mol=$(grep -w -m 1 $name metadata | cut -f $moli)
		local pair=$(grep -w -m 1 $name metadata | cut -f $pairi)
		local cmd
		# select unique (and concordant) mappers
		if [[ $pair == NA ]]; then
			cmd=$(cat <<- END
				if [ $mol == bsDNA ]; then
					bioawk -c fastx '{ gsub("C", "T", \$seq); print "@"\$name" "\$comment; print \$seq; print "+"; print \$qual;}' $indir/$name.fastq.gz > tmp_$name.fastq
					hisat2 $allparam -U tmp_$name.fastq --un-gz $name\_un.fastq.gz -S tmp$name.sam;
				else 
					hisat2 $allparam -U $indir/$name.fastq.gz --un-gz $name\_un.fastq.gz -S tmp$name.sam;
				fi

				samtools view -H tmp$name.sam > $name.sam;
				grep -w NH:i:1 tmp$name.sam >> $name.sam;
				samtools view -@ ${o[nc]} -F 0x4 tmp$name.sam | grep -vw NH:i:1 | cut -f 1 | sort -T . | uniq > tmpmulti$name.ids;

				END
				)
		else
			# correct for ATAC-seq compatible fragment length selection from paired end reads
			# patch for hisat -un-conc-gz option, which gives redundant output; 
			cmd=$(cat <<- END
				if [ $mol == bsDNA ]; then
					bioawk -c fastx '{ gsub("C", "T", \$seq); print "@"\$name" "\$comment; print \$seq; print "+"; print \$qual;}' $indir/${name}_1.fastq.gz > tmp_${name}_1.fastq
					bioawk -c fastx '{ gsub("C", "T", \$seq); print "@"\$name" "\$comment; print \$seq; print "+"; print \$qual;}' $indir/${name}_2.fastq.gz > tmp_${name}_2.fastq
					hisat2 $allparam -1 tmp_${name}_1.fastq -2 tmp_${name}_2.fastq -S tmp$name.sam; 
				else
					hisat2 $allparam -1 $indir/${name}_1.fastq.gz -2 $indir/${name}_2.fastq.gz -S tmp$name.sam; 
				fi
				samtools fastq -F 0x2 -1 ${name}_un_1.fastq -2 ${name}_un_2.fastq tmp$name.sam;
				gzip -f ${name}_un_[12].fastq 

				samtools view -H tmp$name.sam > $name.sam;
				samtools view -@ ${o[nc]} -f 0x2 tmp$name.sam | grep -w NH:i:1 >> $name.sam;
				samtools view -@ ${o[nc]} -f 0x2 tmp$name.sam | grep -vw NH:i:1 | cut -f 1 | sort -T . | uniq > tmpmulti$name.ids;

				END
				)
		fi

		local qcmd=$(cat <<- END
			$cmd

			samtools view -b -@ ${o[nc]} -o $name.bam $name.sam;
			samtools sort -@ ${o[nc]} -o $name.sorted.bam $name.bam;
			samtools index $name.sorted.bam; 

			if [[ ${o[srt]:-F} == T ]]; then
				samtools sort -n -@ ${o[nc]} -o $name.srt.bam $name.bam;
			fi

			if [[ ${o[bamtobed]:-F} == T ]]; then
				if [[ $pair == NA ]]; then
					bedtools bamtobed -i $name.sorted.bam | gzip -c > $name.bed.gz;
				else
					bedtools bamtobed -bedpe -i $name.sorted.bam | gzip -c > $name.bedpe.gz;
					if [[ ${o[bamtobedpe]:-${o[bamtobed]:-F}} == F ]]; then
						zless $name.bedpe.gz | awk '\$1==\$4 && \$6-\$2 < 2000{ OFS="\t"; print \$1,\$2,\$3,\$7,\$8,\$9; print \$4,\$5,\$6,\$7,\$8,\$10; }' | gzip -c > $name.bed.gz;
						rm $name.bedpe.gz;
					fi
				fi
			fi

			if [[ ${o[savemm]:-F} == T ]]; then
				samtools view -b -@ ${o[nc]} -o all_$name.bam tmp$name.sam;
			fi

			# statistics
			echo single- and multi-locus mapping statistics:;
			less $name.sam | grep -v "^@" | cut -f 1 | sort -T . | uniq | wc -l;
			less tmpmulti$name.ids | cut -f 1 | sort -T . | uniq | wc -l;

			#cleaning
			rm -rf $name.{b,s}am tmpmulti$name.ids tmp$name.{b,s}am;


			bedtools --version;
			hisat2 --version;
			samtools --version;
		END
		)
		qcmds+=( [$analysisdirbase\_$name]="$qcmd" )

		# dependent
		qcmd=$(cat <<- END
			sodir=fastqc; mkdir -p \$sodir;
			fastqc -q -t ${o[nc]} -o \$sodir -d \$sodir $name.sorted.bam
			fastqc --version
		END
		)
		dqcmds+=( [QC$analysisdirbase\_$name]="$qcmd" )

	done 
}

slurmMain RM --datadir=~/Projects/biter_jbrett_DNAmethylation/data/Brett_2014/WGBseq --nc=16 --t=12 --indir=TRdelme --genome=mm10_G2A --otag=G2A --skippostproc=F --parseOutput=$po 
slurmMain RM --datadir=~/Projects/biter_jbrett_DNAmethylation/data/Brett_2014/WGBseq --nc=16 --t=12 --indir=TRdelme --genome=mm10 --skippostproc=F --parseOutput=$po 
slurmMain RM --datadir=~/Projects/biter_jbrett_DNAmethylation/data/Brett_2014/WGBseq --nc=16 --t=12 --indir=TRdelme --genome=mm10_C2T --otag=C2T --skippostproc=F --parseOutput=$po 

