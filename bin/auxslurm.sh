#!/bin/bash - 
#===============================================================================
#
#          FILE: auxslurm.sh
# 
#         USAGE: ./auxslurm.sh
# 
#   DESCRIPTION: SLURM submission scripts for deep sequencing data analysis 
# 
#       OPTIONS: ---
#  REQUIREMENTS: ---
#          BUGS: ---
#         NOTES: ---
#        AUTHOR: Dr. Biter Bilen (bb), biterbilen@yahoo.com
#  ORGANIZATION: Stanford, Palo Alto, USA
#       CREATED: 17/07/2015 18:06
#      REVISION:  ---
#===============================================================================

set -o nounset # Treat unset variables as an error
#set -o verbose # Echoes all commands before executing 
#set -o xtrace  # Echo commands after command-line processing

auxhelp() {
	if [[ -z "${@}" ]]; then 
		echo -e "USAGE:\t\tsh auxslurm.sh <SUB-COMMANDS>"
		echo -e "\nSUB-COMMANDS:"
		echo -e "QA\t\tRaw data quality assessment"
		echo -e "TR\t\tRead trimming"
		echo -e "RM\t\tRead to genome alignment with HISAT2"
		echo -e "GB\t\tUCSC Genome Browser hub upload"
		echo -e "FC\t\tFeature Count"
		echo -e "ER\t\tExpressed Genes with 2-component-GMM"
		echo -e "DE\t\tNormalization and differential comparison"
		echo -e "GS\t\tGSEA"
		echo -e "FF\t\tFind features such as TF and miR binding site counts"
		echo -e "DS\t\tDownload sequence"
		echo -e "?|-h|--help\tPrint this help menu"
	else
		subcommand=$1
		echo -e "USAGE:\t\tsh auxslurm.sh $subcommand [OPTIONS]"
		echo -e "\nOPTIONS:"
		case $subcommand in
			QA)
				echo -e "--datadir=<full/path/to/rawData/folder/w/fastq/gz/and/metadata/file/folder>"
				;;
			TR)
				echo -e "--datadir=<full/path/to/rawData/folder/w/fastq/gz/and/metadata/file/folder>"
				;;
			GB)
				echo -e "--datadir=<full/path/to/rawData/folder/w/fastq/gz/and/metadata/file/folder>"
				;;
			RM)
				echo -e "--datadir=<full/path/to/rawData/folder/w/fastq/gz/and/metadata/file/folder>"
				echo -e "--trim3=<read/3'/end/trim/length>"
				echo -e "--trim5=<read/5'/end/trim/length>"
				;;
			FC)
				echo -e "--datadir=<full/path/to/rawData/folder/w/fastq/gz/and/metadata/file/folder>"
				;;
			ER)
				echo -e "--datadir=<full/path/to/rawData/folder/w/fastq/gz/and/metadata/file/folder>"
				;;
			DE)
				echo -e "--datadir=<full/path/to/rawData/folder/w/fastq/gz/and/metadata/file/folder>"
				;;
			GS)
				echo -e "--datadir=<full/path/to/rawData/folder/w/fastq/gz/and/metadata/file/folder>"
				;;
			DS)
				echo -e "--datadir=<full/path/to/rawData/folder/w/fastq/gz/and/metadata/file/folder>"
				;;
			FF)
				echo -e "--datadir=<full/path/to/rawData/folder/w/fastq/gz/and/metadata/file/folder>"
				;;
			PF)
				echo -e "--datadir=<full/path/to/rawData/folder/w/fastq/gz/and/metadata/file/folder>"
				echo -e "--stag=<treatment/sample/name/tag>"
				;;
			*)
				echo -e "\nERROR: Unknown SUB-COMMAND $subcommand"
				exit
				;;
		esac
		[ $# -eq 2 ] && echo -e "\n$2"  # errormessage
	fi
	exit
}

# Removes tailing backslashes and replaces heading ~ with $HOME
auxDirNameTune() {
	echo $(echo $1 | sed 's|/*$||g' | sed "s|~|$HOME|g")
}

auxPrintObject () {
	echo "#OPTIONS"
	local i
	for i in ${!o[@]}; do
		echo -e " [$i]=${o[$i]}"
	done
	echo "#-------"
}

slurmSetObject() {
	#http://mywiki.wooledge.org/BashFAQ/035
	local i
	for i in ${@}; do
		case $i in
			-h|?|--help)
				auxhelp ${o[analysis]}
				;;
			--*=*)
				local varname="${1#--}"
				varname="${varname%=*}"
				varname="${varname%%=*}" #added on 20170729 for filter=state=2
				local value="${1#*=}"
				if [[ $varname == "analysis" ]]; then
					auxhelp ${o[analysis]} "ERROR: OPTION $varname is invalid"
				elif [[ $varname =~ dir$ ]] || [[ $varname =~ file$ ]]; then
					value="$( auxDirNameTune "$value" )"
				fi
				o["$varname"]="$value"
				shift;;
			*)
				auxhelp ${o[analysis]} "ERROR: Unknown OPTION $i";
				;;
		esac
	done

	# Check obligatory params and set default values if not set before
	[ ${o[debug]:-F} == T ] && o[verbose]=T
	[[ -z "${o[datadir]+1}" ]] && auxhelp ${o[analysis]} "ERROR: Set OPTION --datadir"
	if [ ${o[scheduler]:-slurm} == slurm ]; then
		o[stiv]=SLURM_ARRAY_TASK_ID 
		o[stimv]=SLURM_ARRAY_TASK_MAX
	else
		echo ERROR: unknown scheduler ${o[scheduler]}
	fi

	# Set SUB-COMMAND independent params
	o[analysislnk]=${o[analysis]}${o[otag]:-""}
	o[analysisdir]=${o[analysisdir]:-"$(pwd)/${o[analysislnk]}_$(awk -F '/' '{ print $(NF-1)"_"$NF; }' <<< ${o[datadir]} )"}

}

fa2fq() {
	echo FUNCTION $FUNCNAME $* $(date):  
	local ifle=$1
	local ofle=$2

	bioawk -c fastx '{ qual=sprintf("%*s",length($seq),""); gsub(/ /,"I",qual); print "@"$name$comment"\n"$seq"\n+\n"qual; }' $ifle | gzip -c > $ofle
}	

# ShortREAD
QA() {
	echo FUNCTION $FUNCNAME $* $(date):  

	local pattern=${o[pattern]:-.fastq.gz}
	local names=$(bioawk -tc hdr 'NR>1{ ORS=" "; if ($Pair == "NA" ) print $Name; else print $Name"_"$Pair; }' metadata)
	local aliases=$(bioawk -tc hdr 'NR>1{ ORS=" "; print $Alias; }' metadata)

	export -f fa2fq
	parallel --delay 1 --plus --tmpdir . -j${o[nc]} --header : --tag --xapply \
		"[ $pattern == .fasta.gz ] && fa2fq indir/{a}$pattern {n}.fastq.gz || ln -sf $(pwd)/indir/{a}$pattern {n}$pattern" \
		::: n $names ::: a $aliases

	Rscript	-e "source('$PI_SCRATCH/Projects/biter_biter_shared/bin/auxShortRead.R')" \
					-e "caller.qa(dir(pattern='.fastq.gz'), outDir='${o[analysisdir]}', nc=${o[nc]})"

	parallel --version
	bioawk --version
}


QC() {
	echo FUNCTION $FUNCNAME $* $(date):  

	local sodir=${o[sodir]:-fastqc}
	local param="-q -t ${o[nc]} -o $sodir -d $sodir" 
	[ -e "${o[limitsfile]:-SETME}" ] && param+=" -l ${o[limitsfile]}"
	local name=$(bioawk -tc hdr -v i=$((${o[stiv]}+1)) 'NR>1 && $Pair != 2 { k+=1; if (k==i) print $Name; }' metadata)
	local files=$(ls $(pwd)/$name*${o[QCpattern]})

	echo NAME $name files $files
	mkdir -p $sodir
	parallel --delay 1 --plus --tmpdir . -j0 --header : --xapply fastqc $param {f} ::: f $files

	parallel --version
	fastqc --version
}

fq2trimseqr() {
	echo FUNCTION $FUNCNAME $* $(date):  
	local ifle=$1
	local ofle=$2
	local len=$3

	mkdir -p ${ofle%/*} # just in case

	bioawk -c fastx -v len=$len '{ l=len<length($seq)?len:length($seq); print "@"$name$comment"\n"substr($seq,1,l)"\n+\n"substr($qual,1,l); }' $ifle | \
		gzip -c > $ofle

}	

fq2trimseq5terminalN() {
	echo FUNCTION $FUNCNAME $* $(date):  
	local ifle=$1
	local ofle=$2
	local nuc=$3
	local minlen=${4:-15}

	mkdir -p ${ofle%/*} # just in case

	bioawk -c fastx -v ofle=$ofle -v nuc=$nuc -v minlen=$minlen -v tag=$FUNCNAME 'BEGIN{k=N=0;}{ l1=length($seq); gsub("^"nuc"*","",$seq); l=length($seq); if(l1>l) k++; if (l >= minlen) { N++; print "@"$name$comment"\n"$seq"\n+\n"substr($qual,1+(l1-l)); } } END{ ofles=ofle".stat"; OFS="\t"; print "#"tag > ofles; print "TotalReads", "EffectiveReads", "TrimmedReads"> ofles; print NR,N,k > ofles; }' $ifle | \
		gzip -c > $ofle 
	cat $ofle.stat #log to stdout 
}	

fq2len() {
	echo FUNCTION $FUNCNAME $* $(date):  
	for f in $1*.fastq.gz; do
		base=$(basename $f .fastq.gz)
		bioawk -c fastx '{ print(length($seq)) }' $f | sort | uniq -c | sort -k 2,2g | \
			awk -v name=$base 'BEGIN{OFS="\t"; print "Len", "Count", "Name"; }{ print $2,$1,name; }' 
	done
}

TR() {
	echo FUNCTION $FUNCNAME $* $(date):  

	local name=$(bioawk -tc hdr -v i=$((${o[stiv]}+1)) 'NR>1 && $Pair != 2 { k+=1; if (k==i) print $Name; }' metadata)

	echo NAME $name
	export -f fq2trimseqr fq2trimseq5terminalN 
	# fq2trimseqr 4 pairs
	local len=${o[len]:-0} # default:no trim
	local files=($(ls $(pwd)/indir/$name*.fastq.gz)) 
	local tmpdir=tmp_fq2trimseqr_$name
	mkdir -p $tmpdir
	parallel --delay 1 --plus --tmpdir . -j0 --header : --xapply \
		"[ $len -gt 0 ] && fq2trimseqr {f} $tmpdir/{t} $len || ln -sf {f} $tmpdir/{t}" \
		::: t ${files[@]##*/} ::: f ${files[@]}

	# fq2trimseq5terminalN 4 pairs
	local trim5termN=${o[trim5termN]:-X} #default:no trim
	local minlen=${o[m]:-15}
	files=($(ls $(pwd)/$tmpdir/$name*.fastq.gz)) 
	tmpdir=tmp_fq2trimseq5terminalN_$name
	mkdir -p $tmpdir
	parallel --delay 1 --plus --tmpdir . -j0 --header : --xapply \
		"[ $trim5termN == X ] && ln -sf {f} $tmpdir/{t} || fq2trimseq5terminalN {f} $tmpdir/{t} $trim5termN $minlen" \
		::: t ${files[@]##*/} ::: f ${files[@]}

	# cutadapt
	local skipcutadapt=${o[skipcutadapt]:-F}
	indir=tmp_fq2trimseq5terminalN_$name
	if [ $skipcutadapt == F ]; then
		local pair=$(bioawk -tc hdr -v i=$((${o[stiv]}+1)) 'NR>1 && $Pair != 2 { k+=1; if (k==i) print $Pair; }' metadata)
		local as=$(sed 's/;/ -a /g' <<< ";${o[as]:-TruSeqIndexedAdaptor:AGATCGGAAGAGC;A6read1:'A{6}$';T6read1:'T{6}$'}" |  sed 's/:/=/g' )
		local gs=$(sed 's/;/ -g /g' <<< ";${o[gs]:-A6read1:'^A{6}';T6read1:'^T{6}'}" | sed 's/:/=/g' )
		local As=$(sed 's/;/ -A /g' <<< ";${o[As]:-TruSeqIndexedAdaptor:AGATCGGAAGAGC;A6read2:'A{6}$';T6read2:'T{6}$'}" |  sed 's/:/=/g' )
		local param="-q 33,33 --trim-n -m ${o[m]:-15} -n ${o[n]:-10} $([ ${o[gs]:-F} == F ] && echo $as || echo $gs)"
		[ $pair == NA ] && 
			param+=" -o $name.fastq.gz $indir/$name.fastq.gz" ||
			param+=" $As -o ${name}_1.fastq.gz -p ${name}_2.fastq.gz $indir/${name}_1.fastq.gz $indir/${name}_2.fastq.gz"
		source $PI_HOME/PythonSandbox/cutadapt/bin/activate
		cutadapt $param

		echo cutadapt-version $(cutadapt --version)
		deactivate
	else
		mv indir/$name*.fastq.gz .
	fi

	fq2len $name | gzip -c > $name.read_length_stat.txt.gz

	#cleaning
	rm -rf tmp*$name

	bioawk --version
}

# reduced representation
fq2rrfq() {
	echo FUNCTION $FUNCNAME $* $(date):  
	local indir=$1
	local name=$2
	local debug=${3:-F}
	local files=($(ls $indir/$name*.fastq.gz | sort))

	# first pair
	base=$(basename ${files[0]} .gz)
	if [ $debug == T ]; then
		less ${files[0]} | head -n 40000 |	bioawk -c fastx '{ gsub("C","T", $seq); print "@"$name" "$comment; print $seq; print "+"; print $qual;}' -
	else
		bioawk -c fastx '{ gsub("C","T", $seq); print "@"$name" "$comment; print $seq; print "+"; print $qual;}' ${files[0]}
	fi > $base

	# second pair
	if [[ ${#files[@]} -gt 1 ]]; then
		base=$(basename ${files[1]} .gz)
		if [ $debug == T ]; then
			less ${files[1]} | head -n 40000 | bioawk -c fastx '{ gsub("G","A", $seq); print "@"$name" "$comment; print $seq; print "+"; print $qual;}' -
		else
			bioawk -c fastx '{ gsub("G","A", $seq); print "@"$name" "$comment; print $seq; print "+"; print $qual;}' ${files[1]}
		fi > $base
	fi
}

#examples: count|count,[S|M|U|]sam|count,[S|M|U|]id
samSplitMapperTypes() {
	echo FUNCTION $FUNCNAME $* $(date):  
	local insam=$1;
	local type=${2:-count};
	local gz=${3:-T}

	if [[ ! $insam =~ srt.sam ]]; then
		echo Warning: samSplitMapperTypes requires name sorted sam input for correct counting
	fi
	local tag=$(basename $insam .sam)

	bioawk -c sam -v type=,$type -v tag=$tag '
	BEGIN{ OFS="\t"; prev=""; ctype=""; }
	{
		if (and($flag,4)) {
			ctype="U";
		} else {
			if ($0 ~ /NH:i:1$/) {
				ctype="S";
			} else {
				ctype="M";
			} 
		}
		# id and count
		file=tag"."ctype"id";
		pattern=","ctype"id"
		if (prev != $qname) {
			count[ctype]++;
			if (type ~ /,id/ || type ~ pattern ) 
			print $qname > file;
		}

		# sam
		file=tag"."ctype".sam";
		pattern=","ctype"sam"
		if (type ~ /,sam/ || type ~ pattern )
			print $0 > file; 

			prev=$qname;
	} END { if (type ~ /count/) { print tag" read mapping statistics:"; print "singlemapping","multimapping","unmapped"; print count["S"],count["M"],count["U"]; } }' $insam

	idfiles=($( ls $tag.{M,U,S}id 2> /dev/null ))
	if [[ $gz == T ]] && [[ ${#idfiles[@]} -gt 0 ]]; then
		echo ${idfiles[@]} | xargs gzip
	fi
}

bssam2merge () {
	echo FUNCTION $FUNCNAME $* $(date):  
	local insam1=$1
	local insam2=$2
	local otag=$3

	# alphabetic sorting for join operation
	local itag=$otag.asrt

	if [[ ! $insam1 =~ asrt.sam ]]; then
		sort -T . -k 1,1 $insam1 > 1$itag.sam
	else
		[ ! -e 1$itag.sam] && ln -sf $insam2 2$itag.sam
	fi
	if [[ ! $insam2 =~ asrt.sam ]]; then
		sort -T . -k 1,1 $insam2 > 2$itag.sam
	else
		[ ! -e 2$itag.sam] && ln -sf $insam2 2$itag.sam
	fi

	echo $(date) id types;
	samSplitMapperTypes 1$itag.sam count,id,sam F
	samSplitMapperTypes 2$itag.sam count,id,sam F

	echo $(date) id intersections;
	join 1$itag.Sid 2$itag.Sid > $itag.SSid
	join 1$itag.Sid 2$itag.Mid > $itag.SMid
	join 1$itag.Sid 2$itag.Uid > $itag.SUid
	join 1$itag.Mid 2$itag.Sid > $itag.MSid
	join 1$itag.Mid 2$itag.Mid > $itag.MMid
	join 1$itag.Mid 2$itag.Uid > $itag.MUid
	join 1$itag.Uid 2$itag.Sid > $itag.USid
	join 1$itag.Uid 2$itag.Mid > $itag.UMid
	join 1$itag.Uid 2$itag.Uid > $itag.UUid

	# header 
	echo $(date) sam join
	join -t $'\t' 1$itag.S.sam $itag.SSid | sed 's/NH:i:1$/NH:i:2/' >> $otag.sam
	join -t $'\t' 2$itag.S.sam $itag.SSid | sed 's/NH:i:1$/NH:i:2/' >> $otag.sam 
	join -t $'\t' 1$itag.S.sam $itag.SMid | sed 's/NH:i:/NH:i:x/' >> $otag.sam
	join -t $'\t' 2$itag.M.sam $itag.SMid | sed 's/NH:i:/NH:i:x/' >> $otag.sam 
	join -t $'\t' 1$itag.S.sam $itag.SUid >> $otag.sam
	join -t $'\t' 1$itag.M.sam $itag.MSid | sed 's/NH:i:/NH:i:x/' >> $otag.sam
	join -t $'\t' 2$itag.S.sam $itag.MSid | sed 's/NH:i:/NH:i:x/' >> $otag.sam
	join -t $'\t' 1$itag.M.sam $itag.MMid | sed 's/NH:i:/NH:i:x/' >> $otag.sam
	join -t $'\t' 2$itag.M.sam $itag.MMid | sed 's/NH:i:/NH:i:x/' >> $otag.sam
	join -t $'\t' 1$itag.M.sam $itag.MUid >> $otag.sam
	join -t $'\t' 2$itag.S.sam $itag.USid >> $otag.sam
	join -t $'\t' 2$itag.M.sam $itag.UMid >> $otag.sam
	join -t $'\t' 1$itag.U.sam $itag.UUid >> $otag.sam

	# cleaning
	rm {1,2}$itag.* $itag.{S,M,U}{S,M,U}id

	echo END $FUNCNAME $(date)
}

hisat2header() {
	echo FUNCTION $FUNCNAME $* $(date):  
	local name=$1
	local param="$2"
	local ofile=$3
	parallel --tmpdir . 'head -n 4 {} > tmp{}' ::: $name*.fastq
	parallel --tmpdir . 'mv tmp{} {}' ::: $name*.fastq
	hisat2 $param 2> /dev/null | samtools view -H - | grep ^@ > $ofile
	hisat2 $(echo $param | sed "s/_C2T/_G2A/") 2> /dev/null | samtools view -H - | grep ^@PG >> $ofile

	# cleaning
	rm $name*.fastq
}

RM() {
	echo FUNCTION $FUNCNAME $* $(date):  

	local genome=${o[genome]:-mm10}
	local rmdup=${o[rmdup]:-F}

	# HISAT2 - select unique (and concordant) mappers
	local param="-t -p ${o[nc]} -5 ${o[trim5]:-0} -3 ${o[trim3]:-0} --no-discordant --no-mixed \
		-x ${o[genomeindex]:-$PI_HOME/Data/casco/HISAT2_indexes/$genome/chr}"

	param+="$([ ${o[scoremin]:-X} == X ] && echo "" || echo " --score-min ${o[scoremin]:-X}" )"
	param+="$([ ${o[rdg]:-X} == X ] && echo "" || echo " --rdg ${o[rdg]:-X}" )"
	param+="$([ ${o[rfg]:-X} == X ] && echo "" || echo " --rfg ${o[rfg]:-X}" )"
	param+="$([ ${o[nosoftclipping]:-F} == T ] && echo " --sp 10,20" || echo "")"
	
	local name=$(bioawk -tc hdr -v i=$((${o[stiv]}+1)) 'NR>1 && $Pair != 2 { k+=1; if (k==i) print $Name; }' metadata)
	local pair=$(bioawk -tc hdr -v i=$((${o[stiv]}+1)) 'NR>1 && $Pair != 2 { k+=1; if (k==i) print $Pair; }' metadata)
	local mol=$(bioawk -tc hdr -v i=$((${o[stiv]}+1)) 'NR>1 && $Pair != 2 { k+=1; if (k==i) print $Molecule; }' metadata)
	local lib=$(bioawk -tc hdr -v i=$((${o[stiv]}+1)) 'NR>1 && $Pair != 2 { k+=1; if (k==i) print $LibType; }' metadata)

	if [ ${o[spliced]:-F} == T ]; then
		local hisat2version=$(hisat2 --version | grep version -m 1 | awk  '{ print $NF}')
		local spliceSitesfle=${o[spliceSites]:-$PI_HOME/Data/casco/UCSC_tracks/$genome/hisat2v${hisat2version}_spliceSites.txt}
		[[ ! -e ${spliceSitesfle##*/} ]] && ln -sf $spliceSitesfle .
		param+=" --known-splicesite-infile $spliceSitesfle"
	else 
		param+=" --no-spliced-alignment --maxins ${o[maxins]:-500}"
	fi

	echo NAME $name
	export -f fq2rrfq samSplitMapperTypes bssam2merge bssam2merge hisat2header

	local debug=${o[debug]:-F}
	local unique=${o[unique]:-T}
	local saveall=$([ $unique == F ] && echo T || echo ${o[saveall]:-F})
	param+="$([ $lib == unstranded ] && echo "" || echo " --rna-strandness $lib ")"

	if [ $mol == bsDNA ]; then
		fq2rrfq indir $name $debug
		param+=" --no-head $([ $pair == NA ] && echo " -U $name.fastq" || echo " -1 ${name}_1.fastq -2 ${name}_2.fastq")"
		param=$(sed "s/$genome/${genome}_C2T/" <<< $param)
		echo hisat2 $param -S 1$name.sam
		hisat2 $param -S 1$name.sam
		echo hisat2 $(sed "s/_C2T/_G2A/" <<< $param) -S 2$name.sam
		hisat2header $name "$param" $name.sam
		bssam2merge 1$name.sam 2$name.sam $name

		# cleaning
		rm -rf {1,2}$name.sam
	else 
		param+="$([ $pair == NA ] && echo " -U indir/$name.fastq.gz" || echo " -1 indir/${name}_1.fastq.gz -2 indir/${name}_2.fastq.gz")"
		echo hisat2 $param -S $name.sam
		hisat2 $param -S $name.sam
	fi

	# convert to bam for efficient sorting
	samtools view -b -@ ${o[nc]} -o $name.bam $name.sam
	samtools sort -n -@ ${o[nc]} -o $name.nsrt.bam $name.bam
	rm $name.bam

	# convert to sam for samSplitMapperTypes
	samtools view -@ ${o[nc]} -o $name.nsrt.sam $name.nsrt.bam
	#samSplitMapperTypes $name.nsrt.sam count,Ssam,Uid 
	samSplitMapperTypes $name.nsrt.sam count,Ssam

	if [ $unique == T ]; then
		[ $saveall == T ] && mv $name.nsrt.bam $name.bam
		samtools view -H $name.sam > $name.nsrt.sam
		cat $name.nsrt.S.sam >> $name.nsrt.sam
		samtools view -b -@ ${o[nc]} -o $name.nsrt.bam $name.nsrt.sam
		rm $name.nsrt.S.sam
	else
		# $name.bam has the same content as $name.nsrt.bam
		[ $saveall == T ] && ln -sf $name.nsrt.bam $name.bam
	fi

	samtools sort -@ ${o[nc]} -o $name.sorted.bam $name.nsrt.bam
	samtools index $name.sorted.bam 

	if [ $rmdup == T ] && [ $pair == NA ]; then
		samtools rmdup -s $name.sorted.bam $name.ddsorted.bam
		samtools index $name.ddsorted.bam 
		samtools sort -n -@ ${o[nc]} -o $name.ddnsrt.bam $name.ddsorted.bam
		# TODO clean .nsrt.bam* .sorted.bam* if rmdup-if is executed?
	fi

	#cleaning
	rm -rf *$name*.sam 

	bioawk --version
	hisat2 --version
	samtools --version

}

bam2flselbam() {
	local bam=$1
	local obam=$2
	local minins=${3:-0}
	local maxins=${4:-500}
	local tsam=${5:-$obam.sam}
	local nc=${6:-1}
	
	samtools view -H $bam > $tsam.sam
	samtools view -@ $nc $bam > $tsam
	bioawk -tc sam -v minins=$minins -v maxins=$maxins '($tlen > 0 && $tlen <= maxins && $tlen >= minins) || ($tlen < 0 && -$tlen <= maxins && -$tlen >= minins) { print }' $tsam >> $tsam.sam
	samtools view -@ $nc -b -o $obam $tsam.sam
}

FC() {
	echo FUNCTION $FUNCNAME $* $(date):  

	local ftag=${o[ftag]:-_}
	local names=($(bioawk -tc hdr -v sr=${o[samplerm]:-F} 'NR>1 && $Pair != 2 && ($Remove == 0 || sr == "F") { print $Name; }' metadata))
	local pairs=($(bioawk -tc hdr -v sr=${o[samplerm]:-F} 'NR>1 && $Pair != 2 && ($Remove == 0 || sr == "F") { print $Pair; }' metadata))
	local libs=($(bioawk -tc hdr -v sr=${o[samplerm]:-F} 'NR>1 && $Pair != 2 && ($Remove == 0 || sr == "F") { print $LibType; }' metadata))
	local genome=${o[genome]:-mm10}
	local canTXfle=$PI_HOME/Data/casco/UCSC_tracks/$genome/${genome}_refGene_canonical.v2.gtf.gz
	local removeSpike=${o[removeSpike]:-F}
	local spikeId=${o[spikeId]:-ERCC}
	local binsize=${o[binsize]:-NA}
	local bamsuf=${o[bamsuf]:-.nsrt.bam}
	#local ignoreDup=${o[ignoreDup]:-F} #does not work with HISAT2 alignments
	local plotfrac=${o[plotfrac]:-1}

	local param="-T ${o[nc]} -s $(sed -e 's/unstranded/0/' -e 's/RF/2/' -e 's/FR/1/' <<< ${libs[0]} )"
	#param+=$([ ${pairs[0]} == NA ] && echo "" || echo " -p -B -C -d ${o[minins]:-0} -D ${o[maxins]:-500} ")
	# patch for -d -D options of featureCounts
	if [ ${pairs[0]} == NA ]; then
		ln -sf indir tmp
	else
		param+=" -p -B -C"
		if [ ${o[minins]:-0} -eq 0 ] && [ ${o[maxins]:-500} -eq 500 ]; then
			ln -sf indir tmp
		else
			mkdir -p tmp
			# select reads based on fragment length
			export -f bam2flselbam
			parallel --delay 1 --plus --tmpdir . -j${o[nc]} --header : --tag --xapply \
				bam2flselbam indir/{fle} tmp/{fle} ${o[minins]:-0} ${o[maxins]:-500} tmp/{fle}.sam \
				::: fle ${names[@]/%/$bamsuf}
		fi
	fi
	if [ $binsize == NA ]; then
		less $canTXfle | grep exon > features.gtf
		param+=" -t exon -g transcript_id -a features.gtf"
	else
		canTXfle=$PI_HOME/Data/casco/UCSC_tracks/$genome/${genome}_bin${binsize}windows.saf.gz
		less $canTXfle > features.saf
		param+=" -F SAF -a features.saf"
	fi
	#does not work with HISAT2 alignments
	#if [ $ignoreDup == T ]; then
	#	param+=" --ignoreDup"
	#fi

	ln -sf $canTXfle .
	param+=" -o featureCounts.txt $(ls tmp/*$ftag*$bamsuf)"

	echo featureCounts $param
	featureCounts $param
	gzip -f featureCounts.txt

	Rscript \
		-e "source('/scratch/PI/casco/Projects/biter_biter_shared/bin/auxslurm.R')" \
		-e "get.EDA4featureCounts('featureCounts.txt.gz', plotfrac=$plotfrac, ddata_cols='ERCC', col.sel.pat='$bamsuf', ofle.tag='featureCounts', plotRatio=T)"

	# clean
	rm -rf features.{gtf,saf} tmp

	featureCounts -v
	bioawk --version

}

ER() {
	echo FUNCTION $FUNCNAME $* $(date):  

	local indirs=($(echo ${o[indir]:-FC} | sed 's/:/ /g')) #Merge if neccessary
	local fle=featureCounts.txt.gz
	local datadirs=($(echo ${o[datadir]} | sed 's/:/ /g')) #Merge if neccessary
	local namecols=($(echo ${o[nameCol]:-Name} | sed 's/:/ /g')) #Merge if neccessary
	local removeSpike=${o[removeSpike]:-F}
	local spikeId=${o[spikeId]:-ERCC}
	local removeSamplePat=${o[removeSamplePat]:-X}
	local bamsuf=${o[bamsuf]:-.nsrt.bam}
	local plotfrac=${o[plotfrac]:-1}

	for i in ${!datadirs[@]}; do 
		# link to save disk space
		ln -sf ${datadirs[$i]}/${indirs[$i]}/$fle ${i}_$fle
		ln -sf ${datadirs[$i]}/${indirs[$i]}/metadata ${i}_metadata
		# rename samples if colName!=Name but Name2
		local header=$(less ${i}_$fle | grep -m 1 ^Gene)
		if [ ${namecols[$i]} != Name ]; then
			local colRepS=$(bioawk -c hdr 'NR>1{ ORS=";"; print "s/"$Name"/"$Name2"/" }' ${i}_metadata);
			header=$(echo $header | sed "$colRepS")
		fi
		# remove or rename Geneid of spikeins
		echo $header > tmp$i
		less ${i}_$fle | grep -v "^#" | bioawk -tc hdr -F '\t' -v rs=$removeSpike -v si=$spikeId \
			'rs=="F" || $Chr !~ si { if ($Chr ~ si) { $Geneid=$Chr; } if ($Chr != "Chr") print }' >> tmp$i
	done
	paste tmp* | gzip -c > $fle;

	Rscript \
		-e "source('/scratch/PI/casco/Projects/biter_biter_shared/bin/auxslurm.R')" \
		-e "get.EDA4featureCounts('$fle', plotfrac=$plotfrac, ddata_cols='ERCC', col.sel.pat='$bamsuf', ofle.tag='expressedFeatureCounts', GMM2expressed=T)" \
		-e "get.EDA4featureCounts('$fle', plotfrac=$plotfrac, ddata_cols='ERCC',col.sel.pat='$bamsuf', ofle.tag='TMscaledExpressedFeatureCounts', GMM2expressed=T, trimmedMeanScaled=T)"

	# clean 
	rm tmp*
}

writeGmt() {
	local otag=$1
	local FDRcut=$2
	local tag=$3

	less $otag.txt.gz | grep -v "^#" | \
		bioawk -c hdr -v tag="custom:$tag down" -v cut=$FDRcut \
		'BEGIN{ ORS=""; print tag"\t"tag; } $6<cut && $2<0 { print "\t"$ENTREZID } END{ print "\n"; }' > $otag.gmt

	less $otag.txt.gz | grep -v "^#" | \
		bioawk -c hdr -v tag="custom:$tag up" -v cut=$FDRcut \
		'BEGIN{ ORS=""; print tag"\t"tag; } $6<cut && $2>0 { print "\t"$ENTREZID } END{ print "\n"; }' >> $otag.gmt 

	less $otag.txt.gz | grep -v "^#" | \
		bioawk -c hdr -v tag="custom:$tag updown" -v cut=$FDRcut \
		'BEGIN{ ORS=""; print tag"\t"tag; } $6<cut         { print "\t"$ENTREZID } END{ print "\n"; }' >> $otag.gmt 

	gzip $otag.gmt
}

DE() {
	echo FUNCTION $FUNCNAME $* $(date):  

	local fle=${o[datadir]}"/${o[indir]}/featureCounts.txt.gz"
	local efle=${o[datadir]}"/${o[indir]}/expressedFeatureCounts.id"
	local species=${o[species]:-Mm}
	local contrast=${o[contrast]:-c(1,-1)}
	local bamsuf=${o[bamsuf]:-.nsrt.bam}
	local colselpat=${o[colselpat]:-$bamsuf}
	local FDRcut=${o[FDRcut]:-0.01}
	local GSEApcut=${o[GSEApcut]:-0.01}
	local otag=edgeR_DE

	ln -sf $fle .
	ln -sf $efle .
	Rscript \
		-e "source('$PI_SCRATCH/Projects/biter_biter_shared/bin/auxslurm.R')" \
		-e "param <- list(out='top', scale='upperquartile', FDRcut=$FDRcut, species='$species', contrast=$contrast)" \
		-e "top <- caller.DE.edgeR(fle='$fle', efle='$efle', col.sel.pat='$colselpat', suf.pat='$bamsuf', param=param, tag='$otag', plot2file=T)" 

	writeGmt $otag $FDRcut ${o[analysisdir]##*/}

}

GS() {
	echo FUNCTION $FUNCNAME $* $(date):  

	local fle=${o[datadir]}"/${o[indir]}/edgeR_DE.txt.gz"
	local flecs=($(ls ${o[datadir]}/../${o[flecdir]}/DE*/edgeR_DE.gmt.gz))
	local flets=($(ls ${o[datadir]}/FF*/promoters.gmt.gz ${o[datadir]}/FF*/self.gmt.gz ))
	local flems=($(ls ${o[datadir]}/FF*/threeUTRsByTranscript.gmt.gz ${o[datadir]}/FF*/cdsBy.gmt.gz ))
	local flecS="$(echo ${flecs[@]} | sed "s/ /','/g")"
	local fletS="$(echo ${flets[@]} | sed "s/ /','/g")"
	local flemS="$(echo ${flems[@]} | sed "s/ /','/g")"
	local species=${o[species]:-Mm}
	local FDRcut=${o[FDRcut]:-0.01}
	local stat=${o[stat]:-fisher}
	local GSEApcut=${o[GSEApcut]:-0.05}

	ln -sf $fle .

	Rscript \
		-e "source('$PI_SCRATCH/Projects/biter_biter_shared/bin/auxslurm.R')" \
		-e "param <- list(FDRcut=$FDRcut, GSEApcut=$GSEApcut, species='$species', stat='$stat')" \
		-e "flec <- c('$flecS')" \
		-e "flet <- c('$fletS')" \
		-e "flem <- c('$flemS')" \
		-e "tmp <- caller.GSEA('$fle', plot2file=T, param=param, types=c('goana'), tag='GSEA_GO')" \
		-e "tmp <- caller.GSEA('$fle', plot2file=T, param=param, types=c('kegga'), tag='GSEA_KEGG')" \
		-e "if (all(file.exists(flec))) tmp <- caller.GSEA('$fle', flec=flec, param=param, plot2file=T, types=c('customa'), tag='GSEA_CUSTOMDE')" \
		-e "if (all(file.exists(flet))) for (i in 1:length(flet)) { tag <- paste0('GSEA_CUSTOM', basename(dirname(flet[i])),gsub('.gmt.gz', '', basename(flet[i]))); print(tag); print(flet[i]); tmp <- caller.GSEA('$fle', flec=flet[i], param=param, plot2file=T, types=c('customa'), tag=tag) }" \
		-e "if (all(file.exists(flem))) for (i in 1:length(flem)) { tag <- paste0('GSEA_CUSTOM', basename(dirname(flem[i])),gsub('.gmt.gz', '', basename(flem[i]))); print(tag); print(flem[i]); tmp <- caller.GSEA('$fle', flec=flem[i], param=param, plot2file=T, types=c('customa'), tag=tag) }"
#		-e "flec <- Sys.glob(paste0('$flecdir','/DE*/*gmt.gz'))" \
}

# TODO move the caller.expressedTxRegionsOrSeqs to auxslurm.R
DS() {
	echo FUNCTION $FUNCNAME $* $(date):  

	local fle=${o[datadir]}"/${o[indir]}/featureCounts.txt.gz"
	local regionfle=${o[regionfile]:-NA}
	local regionsStr=${o[regionsStr]:-self}
	local genome=${o[genome]:-mm10}
	local linesSkipped=${o[linesSkipped]:-0}
	local filterStr="${o[filterStr]:-NA}"
	local posStrand=${o[posStrand]:-F}
	local colsStr="${o[colsStr]:-NULL}"

	if [ $regionfle != NA ]; then
		ln -sf $regionfle .
	else
		ln -sf $fle .
	fi

	Rscript \
		-e "source('$PI_SCRATCH/Projects/biter_biter_shared/bin/auxGenomicFeatures.R')" \
		-e "fle <- ifelse('$regionfle' == 'NA', '$fle', '$regionfle')" \
		-e "filterStr <- ifelse('$filterStr' == 'NA', NA, '$filterStr')" \
		-e "if ('$colsStr' == 'NULL') col.names <- c('chr','start','end','name') else col.names <- strsplit('$colsStr', ',')[[1]]" \
		-e "ofleTag <- ''" \
		-e "regionSeqs <- caller.expressedTxRegionsOrSeqs(expressedTx.file=fle, regions='$regionsStr', filterStr=filterStr, col.names=col.names, gnm='$genome', skip=$linesSkipped, posStrand=$posStrand, ofleTag=ofleTag, debug=${o[debug]}, verbose=T)"
}


# run regions one by one
FF() {
	echo FUNCTION $FUNCNAME $* $(date):  

	local fle=${o[datadir]}"/${o[indir]}/edgeR_DE.txt.gz"
	local minScore=${o[minScore]:-'90%'}
	local miRSeqfle=NA
	local miRExpfle=${o[miRExpfile]:-NA}
	local expFilter=F
	local regionsStr=${o[regionsStr]:-promoters,threeUTRsByTranscript,self}
	local regionfle=${o[regionfile]:-NA}
	local seedStart=${o[seedStart]:-2}
	local seedEnd=${o[seedEnd]:-8}
	if [[ $regionsStr =~ threeUTRsByTranscript ]] || [[ $regionsStr =~ cdsBy ]]; then
		miRSeqfle=$PI_HOME/Data/casco/miRBase/mature.fa.gz
	fi
	local genome=${o[genome]:-mm10}

	ln -sf $fle .
	Rscript \
		-e "source('$PI_SCRATCH/Projects/biter_biter_shared/bin/auxslurm.R')" \
		-e "rfle <- ifelse('$regionfle' == 'NA', NA, '$regionfle')" \
		-e "miRSeq.file <- ifelse('$miRSeqfle' == 'NA', NA, '$miRSeqfle')" \
		-e "miRExp.file <- ifelse('$miRExpfle' == 'NA', NA, '$miRExpfle')" \
		-e "gmt <- ifelse('$regionsStr' == 'self', F, T)" \
		-e "caller.getSeqFeatures(expressedTx.file='$fle',mirSeq.file=miRSeq.file,miRExp.file=miRExp.file,exp.filter=$expFilter,regions.str='$regionsStr',rfle=rfle,min.score='$minScore', seed.start=$seedStart, seed.end=$seedEnd, gnm='$genome',mc.cores=${o[nc]}, gmt=gmt, debug=${o[debug]})"
}

PF() {
	echo FUNCTION $FUNCNAME $* $(date):  

	local names=($(bioawk -tc hdr 'NR>1 && $Pair != 2 { print $Name; }' metadata))
	local name=$(bioawk -tc hdr -v i=$((${o[stiv]}+1)) 'NR>1 && $Pair != 2 { k+=1; if (k==i) print $Name; }' metadata)
	local cname=$(bioawk -tc hdr -v i=$((${o[stiv]}+1)) 'NR>1 && $Pair != 2 { k+=1; if (k==i) print $Cname; }' metadata)
	local broad=${o[broad]:-False}

	echo NAME $name

	source $PI_HOME/PythonSandbox/MACS2/bin/activate
	# predictd
	macs2 predictd -g ${o[gsize]:-mm} -i indir/$name.ddsorted.bam --rfile $name\_predictd.R;
	Rscript ${name}_predictd.R
	# callpeak
	# TODO fixme broad does not work; gives error for False True FALSE TRUE
	# TODO run --cutoff-analysis for the optimum --broad-cutoff
	local extsize=$(less ${name}_predictd.R | grep ^altd | sed 's/[c()]//g' | awk '{ print int($3/2); }')
	if [ $cname != NA ]; then
		macs2 callpeak -c indir/$cname.ddsorted.bam -t indir/$name.ddsorted.bam -n $name -g ${o[gsize]:-mm} \
			-q ${o[FDRcut]:-0.01} --nomodel --extsize ${o[extsize]:-$extsize} -B --trackline --SPMR #--broad $broad --broad-cutoff 0.1 
	fi
			          
	macs2 --version
	deactivate
}

# TODO
# atac-seq pipeline: https://docs.google.com/document/d/1f0Cm4vRyDQDu0bMehHD7P7KOMxTOP-HiNoIvL1VcBt8/edit#heading=h.z5vzv9ujdu5r

# Accepts unique (and concordantly paired) mapped reads as input
# TODO fragment coverage for paired-end; currently read coverage 
# https://genome.ucsc.edu/goldenpath/help/trackDb/trackDbDoc.html#wig_-_Signal_Graphing_Track_Settings

# TODO make the rest array job compatible
GB() {
	echo FUNCTION $FUNCNAME $* $(date):  
	local analysisdirbase=$1

	# copy with links archive compress permissions 
	#-----------------------
	local copy2server=${o[copy2server]:-F}
	local webserverdir=${o[webserverdir]:-"biter@cardinal.stanford.edu:/afs/ir/group/rando_lab/WWW/biter"}
	if [[ "$copy2server" == T ]]; then
		rsync -e 'ssh -p 22' -Lazp ${o[analysisdir]} $webserverdir 
		return
	fi

	#-----------------------
	local indir=${o[indir]:-RM}
	local names=($(bioawk -tc hdr 'NR>1 && $Pair != 2 { print $Name; }' metadata))
	local genome=${o[genome]:-mm10}
	local genomel=${o[genomel]:-mouse.mm10}
	local selchr=${o[selchr]:-chr19}
	local bigDataUrl="${o[bigDataUrl]:-"https://web.stanford.edu/group/rando_lab/biter"}/$analysisdirbase"
	local genomefafile=${o[genomefafile]:-"$PI_HOME/Data/casco/HISAT2_indexes/$genome/all/chr.fa"}
	local auxgenome=${o[auxgenome]:-$HOME/aux/genomes/$genomel.genome} 
	local email=${o[email]:-"biter@stanford.edu"}
	local tracksfle=${o[tracksfle]:-tracks.txt}
	local genomesfle=${o[genomesfle]:-genomes.txt}
	local	visibilitybw=${o[visibilitybw]:-full} 
	local	visibilitybam=${o[visibilitybam]:-hide} 
	local	maxHeightPixels=${o[maxHeightPixels]:-50:30:11}
	local viewLimitsmin=${o[viewLimitsmin]:-3}
	local viewLimitsmax=${o[viewLimitsmax]:-50}
	local color=${o[color]:-0,0,0}
	local position=${o[position]:-"$selchr:212000-250000"}
	# hub.txt
	cat > hub.txt <<- END
	hub $analysisdirbase
	shortLabel $analysisdirbase
	longLabel $analysisdirbase
	genomesFile $genomesfle
	email $email
	END

	# genomes.txt
	cat > genomes.txt <<- END
	trackDb $tracksfle
	genome $genome
	END

	# tracks.txt
	echo browser position $position > $tracksfle
	for name in ${names[@]}; do
		cat >> $tracksfle <<- END

		track $name 
		shortLabel $name
		longLabel $name
		type bigWig 
		bigDataUrl $bigDataUrl/$name.bw 
		visibility $visibilitybw
		maxHeightPixels $maxHeightPixels
		viewLimits $viewLimitsmin:$viewLimitsmax
		color $color

		track BAM_$name
		shortLabel BAM_$name
		longLabel BAM_$name
		type bam 
		bigDataUrl $bigDataUrl/${selchr}_$name.sorted.bam 
		visibility $visibilitybam 
		END
	done

	local qcmd=$(cat <<- END
	nc=${o[nc]}
	datadir=${o[datadir]}/$indir
	selchr=$selchr
	auxgenome=$auxgenome
	bigDataUrl=$bigDataUrl
	names=(${names[@]})
	name=\${names[\$${o[stiv]}]}

	bedtools genomecov -bg -ibam \$datadir/\$name.sorted.bam -g \$auxgenome > \$name.bg
	bedGraphToBigWig \$name.bg \$auxgenome \$name.bw

	samtools view -@ \$nc -b -o \${selchr}_\$name.sorted.bam \$datadir/\$name.sorted.bam \$selchr
	samtools index \${selchr}_\$name.sorted.bam

	rm -rf \$name.bg \${selchr}_\$name.bam

	samtools --version
	bedtools --version
	bedGraphToBigWig
	END
	)
	jobnames+=($analysisdirbase)
	qcmds[$analysisdirbase]="$qcmd"
	qarrs[$analysisdirbase]=$(( ${#names[@]} - 1 ))
	qwpts[$analysisdirbase]=F

}


slurmParseOutput() {
	pushd ${o[analysisdir]}
	echo 
	local f
	local tag
	local files
	if [[ ${o[analysis]} == DE ]]; then
		echo Name Down Up DownUp Universe | sed 's/ /\t/g' > summary
		f=$(ls slurm*)
		tag=$(pwd | awk -F "/" '{ print $NF}')
		local d=$(grep "DE at FDR" -A 5 $f | awk -v type=-1 'BEGIN{k=0} $1==type { k=$3} END{ print k; }')
		local u=$(grep "DE at FDR" -A 5 $f | awk -v type=1 'BEGIN{k=0} $1==type { k=$3} END{ print k; }')
		local o=$(grep "DE at FDR" -A 5 $f | awk -v type=0 'BEGIN{k=0} $1==type { k=$3} END{ print k; }')
		echo $d $u $o $tag | awk '{OFS="\t"; print $4,$1,$2,$1+$2,$1+$2+$3}' >> summary
	elif [[ ${o[analysis]} == ER ]]; then
		echo Name ExpressedFeatures ExpressedFeaturesPercent | sed 's/ /\t/g' > summary
		f=$(ls slurm*)
		tag=$(pwd | awk -F "/" '{ print $NF}')
		local ef=($(grep "Clustering table:" $f -A 2 | awk 'NR % 4 == 3 { print $2}'))
		local efp=($(grep "Mixing probabilities:" $f -A 2 | awk 'NR % 4 == 3 { printf("%.1f\n", 100*$2) }'))
		local N=$((${#efp[@]}/2))
		local n=($(grep -m 1 mGMM -A $N $f | awk 'NR>1 { print $2}'))
		for i in `seq 0 $((N-1))`; do
			echo ${tag}_${n[$i]} ${ef[$i]} ${efp[$i]}
		done | sed 's/ /\t/g' >> summary
	elif [[ ${o[analysis]} == FC ]]; then
		echo Name SingleLocusMappingReads FeatureReads FeatureReadsPercent | sed 's/ /\t/g' > summary
		f=$(ls slurm*)
		tag=$(pwd | awk -F "/" '{ print $NF}')
		grep -P "Process BAM|Total fragments|Total reads|Successfully assigned fragments|Successfully assigned reads" $f | \
			awk '{ if (NR % 3 == 0) { ORS="\n" } else { ORS="\t" } print }' | \
			awk '{ OFS="\t"; print $5,$11,$18,$19 }' | sed 's/[(%)]//g' | sed 's/.nsrt.bam...//g' | \
			awk -F "/" -v tag=$tag '{ print tag"_"$2 }' >> summary
	elif [[ ${o[analysis]} == RM ]]; then
		files=($(ls slurm_${o[analysis]}*))
		if [ "$(grep -m 1 "bam_rmdupse_core" ${files[0]} | awk '{ print $4-$2}')" == "" ]; then
			echo Name TotalReads SingleLocusMappingReads SingleLocusMappingReadPercent MultiLocusMappingReadPercent
		else
			echo Name TotalReads SingleLocusMappingReads deduplicatedSingleLocusMappingReads SingleLocusMappingReadPercent deduplicatedSingleLocusMappingReadPercentPercent MultiLocusMappingReadPercent
		fi | sed 's/ /\t/g' > summary
		for f in ${files[@]}; do
			local n=`grep "^NAME" $f | awk '{ print $2}'`; 
			local t=`grep -m 1 "reads; of these:" $f | head | awk '{ print $1}'`; 
			local s=`grep -A 2 "nsrt read mapping statistics:" $f | tail -n 1 | cut -f 1`; 
			local m=`grep -A 2 "nsrt read mapping statistics:" $f | tail -n 1 | cut -f 2`; 
			local d=`grep -m 1 "bam_rmdupse_core" $f | awk '{ print $4-$2}' `; 
			if [ "$d" == "" ]; then 
				echo $n $t $s $m | awk '{ OFS="\t"; printf("%s\t%.0f\t%.0f\t%.1f\t%.1f\n", $1, $2, $3, $3/$2*100, $4/$2*100) }'; 
			else
				echo $n $t $s $d $m | awk '{ OFS="\t"; printf("%s\t%.0f\t%.0f\t%.0f\t%.1f\t%.1f\t%.1f\n", $1, $2, $3, $4, $3/$2*100, $4/$2*100, $5/$2*100) }'; 
			fi
		done | sort -k 1,1 >> summary 
	elif [[ ${o[analysis]} == TR ]]; then
		echo Name TotalReads EffectiveReads EffectiveReadsPercent ReadswAdaptorPercent | sed 's/ /\t/g' > summary
		for f in slurm_${o[analysis]}*; do 
			local n=`grep "^NAME" $f | awk '{ print $2}'`; 
			local t=`grep -P 'Total reads processed:|Total read pairs processed:' $f | awk '{ print $NF}'`
			local rwa1=`grep -P 'Reads with adapters:|Read 1 with adapter:' $f | awk '{ print $NF}'`
			local rwa2=`grep -P 'Read 2 with adapter:' $f | awk '{ print $NF; }'`
			local rwa=$(echo $rwa1 $rwa2 | sed 's/[()%]//g' | awk '{ if (NF ==  1 ) { print $1; } else { print $1"|"$2; } }')
			local es=`grep -P 'Reads written|Pairs written' $f | awk '{ print $(NF-1),$NF}'` 
			echo $n $t $es $rwa | sed 's/,//g' | sed 's/[()%]//g' | sed 's/ /\t/g'
		done | sort -k 1,1 >> summary 
	else
		echo Implement slurmParseOutput for ${o[analysis]}
	fi

	popd
}

slurmSubmitToQueue() {
	mkdir -p ${o[analysisdir]}; pushd ${o[analysisdir]}
	mkdir -p bin

	jobnames=(${o[analysis]})
	[[ ,QA,TR,RM =~ ,${jobnames[0]} ]] && [ ${o[skipfastqc]:-F} == F ] && jobnames+=(QC)

	local ja=0
	local jids=""
	local sr="-c ${o[nc]} -p ${o[queue]} -t ${o[t]}:${o[tmin]}:0";  # schedular resources
	[ ${o[queue]} == gpu ] && sr+=" --gres gpu:${o[nc]}" #node and gpu core count
	for jn in ${jobnames[@]}; do
		local funcs=""
		case $jn in
			QA)
				funcs="fa2fq"
				o[indir]=${o[indir]:-rawData}
				o[QCpattern]=.fastq.gz
				[ ${o[skipfastqc]:-F} == F ] && jobnames+=(QC)
				;;
			TR)
				funcs="fq2trimseqr fq2trimseq5terminalN fq2len"
				o[indir]=${o[indir]:-QA}
				o[QCpattern]=.fastq.gz
				[ ${o[skipfastqc]:-F} == F ] && jobnames+=(QC)
				;;
			RM)
				funcs="fq2rrfq samSplitMapperTypes bssam2merge bssam2merge hisat2header"
				o[indir]=${o[indir]:-QA}
				o[QCpattern]=.sorted.bam
				[ ${o[skipfastqc]:-F} == F ] && jobnames+=(QC)
				;;
			FC)
				funcs="bam2flselbam"
				o[indir]=${o[indir]:-RM}
				;;
			ER)
				o[indir]=${o[indir]:-FC}
				;;
			DE)
				funcs="writeGmt"
				o[indir]=${o[indir]:-ER}
				;;
			PF)
#				funcs="calcsim auxCalcNucSpacing"
				o[indir]=${o[indir]:-RM}
				;;
			DS)
				o[indir]=${o[indir]:-FC}
				;;
			FF)
				o[indir]=${o[indir]:-DE}
				;;
			GS)
				o[indir]=${o[indir]:-DE}
				;;
			*)
				;;
		esac
		local datadirs=($(echo ${o[datadir]} | sed 's/:/ /g'))
		local ja=0
		if [ ${#datadirs[@]} -eq 1 ]; then
			ln -sf ${o[datadir]}/${o[indir]} indir
			ln -sf ${o[datadir]}/${o[indir]}/metadata metadata
			ja=$(bioawk -tc hdr 'NR>1 && $Pair != 2 { k+=1; } END{ print k-1 }' metadata)
		else
			local indirs=($(echo ${o[indir]} | sed 's/:/ /g'))
			for i in ${!datadirs[@]}; do
				ln -sf ${datadirs[$i]}/${indirs[$i]} ${prefix}indir
				ln -sf ${datadirs[$i]}/${indirs[$i]}/metadata ${prefix}metadata
			done
		fi
		[[ ,QA,FC,ER,DE,GS,FF,DS =~ ,$jn ]] || [ ${o[debug]:-F} == T ] && ja=0 # TODO add others
		local srcfle=bin/$jn${o[analysisdir]##*/}.sh
		local qcmd="sbatch -J ${srcfle##*/} $sr -a 0-$ja -e slurm_${jn}_%A_%a.out -o slurm_${jn}_%A_%a.out"
		[[ ! -z $jids ]] && qcmd+=" --dependency=afterOK$jids" 
		qcmd+=" $srcfle"
		echo Doing $qcmd
		cat > $srcfle <<- END
		#!/bin/bash -l 
		# version=${o[version]:-1.0}
		# date=$(date -I)
		# sbatch command: $qcmd
		date
		hostname
		$(declare -p o)
		$(declare -f $funcs $jn)
		$jn
		date
		echo DONE
		END
		if [ ${o[dryrun]:-F} == F ]; then  
			local sbout=$($qcmd)
			jids=$jids:$(echo $sbout | awk '{print $4}')
		fi
	done

	popd
}

# Makes a link of the analysis folder in data folder
slurmPostProc() {
	pushd ${o[datadir]}
	[ ! -e ${o[analysislnk]} ] && ln -sf ${o[analysisdir]} ${o[analysislnk]}
	popd
}

slurmMain() {

	[[ -z "${@}" ]] && auxhelp

	declare -A o=( [scheduler]=slurm [verbose]=F [debug]=F [analysis]=$1 [queue]=normal [nc]=8 [t]=6 [tmin]=0 )
	#declare -A o=( [scheduler]=slurm [verbose]=F [debug]=F [analysis]=$1 [queue]=casco [nc]=8 [t]=6 [tmin]=0 )
	#declare -A o=( [scheduler]=slurm [verbose]=F [debug]=F [analysis]=$1 [queue]=owners [nc]=8 [t]=6 [tmin]=0 )
	declare -a jobnames

	shift

	slurmSetObject "$@"

	[ ${o[verbose]} == T ] && auxPrintObject

	[ ${o[parseOutput]:-F} == T ] && slurmParseOutput || slurmSubmitToQueue

	[ ${o[skippostproc]:=F} == F ] && slurmPostProc

}

