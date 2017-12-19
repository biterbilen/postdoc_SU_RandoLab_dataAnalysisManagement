#!/bin/bash - 
#===============================================================================
#
#          FILE: cmd_PAS_analysis.sh
# 
#         USAGE: ./cmd_PAS_analysis.sh 
# 
#   DESCRIPTION: 
# 
#       OPTIONS: ---
#  REQUIREMENTS: ---
#          BUGS: ---
#         NOTES: ---
#        AUTHOR: Dr. Biter Bilen (bb), biterbilen@yahoo.com
#  ORGANIZATION: Stanford, Palo Alto, USA
#       CREATED: 05/20/2016 17:23
#      REVISION:  ---
#===============================================================================

set -o nounset                              # Treat unset variables as an error
source ~/Projects/biter_biter_shared/bin/auxslurm.sh; # for auxDirNameTune

#1817358
# 1779437
PASdeprecated() {
	# 1.
	files=($( auxDirNameTune ${o[datadir]}/${o[indir]:-RM} | sed 's/,/\/*.sorted.bam /g' | sed 's/$/\/*.sorted.bam/g' ))
	labels=($( echo ${files[@]##*/} | awk -F / '{ print $NF}' | sed 's/.sorted.bam//g' ))
	ln -sf ${files[@]} .
	jobname=genomecov 
	for i in ${!labels[@]}; do
		FT=${labels[$i]}
		qcmd=$(cat <<- END
			parallel -k --delay 1 --plus --tmpdir . -j0 --header : --tag --joblog $jobname.jlog --xapply \
			"samtools view -b -{SST} 0x10 $FT.sorted.bam | \
			bedtools genomecov -bg -ibam stdin -g ~/aux/genomes/mouse.${o[genome]:-mm10}.genome | \
			awk -v OFS='\t' -v dummy='.' -v strand='{ST}' '{ print $1,$2,$3,dummy,$4,strand }' > $FT{ST}.bed" \
		::: ST + - ::: SST F f  
		END
		)
	done
	if [[ ,${o[jobnames]}, =~ ,$jobname, ]]; then qcmds+=( [$jobname]="$qcmd" ); fi

	# 2.
	jobname=merge
	#cat *bg | bedtools sort -i stdin | bedtools merge -i stdin > commonIntervals.bg
}

nucdistr() {
	local fle=$1
	local slopl=${2:--1900}
	local slopr=${3:--900}
	local len=${4:-200}
	local ofletag=${5:-$fle}
	ofletag=${ofletag##*/}
	ofletag=${ofletag%%.*}
	#echo $FUNCNAME $ofletag
	bedtools slop -s -l $slopl -r $slopr -i $fle -g ~/aux/genomes/mouse.mm10.genome | \
		awk -v len=$len '$3-$2==len { print }' | \
		bedtools getfasta -name -s -fi ~/PI_HOME/Data/casco/UCSC_tracks/mm10/assembly/all/chr.fa -bed stdin | \
		if [[ $ofletag == stdout ]]; then 
			less 
		else 
			gzip -c > $ofletag.fa.gz 
		fi
	if [[ $ofletag != stdout ]]; then
		N=$(bioawk -c fastx 'END{ print NR}' $ofletag.fa.gz)
		Rscript -e 'source("~/Projects/biter_biter_shared/bin/auxGenomicFeatures.R"); ' \
						-e 'dta <- caller.nucProbPerPosition("'$ofletag'.fa.gz", start=1, end='$len')' \
						-e 'source("~/Projects/biter_biter_shared/bin/auxggplot2.R"); ggplot <- set.ggplot();' \
						-e 'ggplot(melt(dta, id="position"), aes(x=position, y=value, color=variable, shape=variable)) + geom_point() + geom_line(size=rel(1)) + ' \
						-e 'labs(y="Probability", x="Location", title=paste0("'$ofletag'", " N=", '$N'));' \
						-e 'ggsave(paste0("'$ofletag'","_nucComp.pdf"))' #2> /dev/null
	fi
}

nuc() {
	local fle=$1
	local slopl=${2:--1950}
	local slopr=${3:--1000}
	local len=${4:-200}
	local pattern=${5:-AATAAA}
	local ofletag=${6:-$fle}
	ofletag=${ofletag##*/}
	ofletag=${ofletag%%.*}
	echo $FUNCNAME $ofletag
	bedtools slop -s -l $slopl -r $slopr -i $fle -g ~/aux/genomes/mouse.mm10.genome | \
		awk -v len=$len '$3-$2==len { print }' | \
		bedtools nuc -C -fullHeader -pattern $pattern -s -fi ~/PI_HOME/Data/casco/HISAT2_indexes/mm10/all/chr.fa -bed stdin | \
		sed 's/^#//g' > $ofletag.nuc
	echo TASK $pattern probability in $ofletag.nuc: $(awk 'BEGIN{k=0;} NR>1{ if ($NF>0) k++; } END{ printf("%.2f", k/(NR-1)); }' $ofletag.nuc )
}

SVD4nucComp() {
	local fle5p=$1
	local fle3p=$2
	local ofletag=$3
	local pcuse=${4:-1:2}
	local feaPat=${5:-num_[AT]}
	# TODO rownames is wrong debug
	Rscript -e 'source("~/Projects/biter_biter_shared/bin/auxggplot2.R"); ggplot <- set.ggplot();' \
					-e 'a <- read.table("'$fle5p'", header=T, comment.char="");' \
					-e 'b <- read.table("'$fle3p'", header=T, comment.char="");' \
					-e 'dta <- merge(a, b, by="X4_usercol", suffix=c(".'$ofletag'5p", ".'$ofletag'3p"));' \
					-e 'rownames(dta) <- dta$X4_usercol; ' \
					-e 'dta <- dta[,c(2,3,4,1,5:ncol(dta))];' \
					-e 'N <- nrow(dta); ' \
					-e 'dta[,grep("'$feaPat'", names(dta))] <- dta[,grep("'$feaPat'", names(dta))] / 50; ' \
					-e 'writeData(dta, "'$ofletag'.txt", row.names=F); ' \
					-e 'get.svdplot(dta[grep("'$feaPat'", names(dta))], ks=c(2,4), title=paste0("'$ofletag'", " N=", N), ofle.tag="'$ofle'_svd", verbose=T, pc.use='$pcuse');' 2> /dev/null
}

closest() {
	local fle=$1
	local annotf=$2
	local slopl=${3:--1999}
	local slopr=${4:--1000}
	local ofletag=${5:-$fle}
	ofletag=${ofletag##*/}
	ofletag=${ofletag%%.*}
	local tmp=_$RANDOM$FUNCNAME$ofletag
	echo $FUNCNAME $ofletag
	bedtools slop -s -l $slopl -r $slopr -i $fle -g ~/aux/genomes/mouse.mm10.genome | bedtools sort -i stdin > $tmp.bed
	less $annotf | bedtools sort -i stdin | bedtools closest -d -s -t first -a $tmp.bed -b stdin > $ofletag.closest
	rm $tmp.bed
}

bamsto5pbed() {
	local indir=$1
	local nc=${2:-1}
	local extractSoftClipped=${3:-F}
	local ftag=${4:-sorted}
	local files=($(ls $indir/*$ftag*bam))
	samtools cat ${files[@]} | \
		if [ $extractSoftClipped == T ]; then 
			samtools view -h | awk '$0 ~ "^@" || $6 ~ /S/ { print }' | samtools view -bS 
		else
			samtools view -bS
		fi | \
		bedtools bamtobed -i stdin | awk 'BEGIN{OFS="\t"}{ if ($6=="+"){ $3=$2+1; } else { $2=$3-1; } print $0; }' | \
		sort -k1,1 -k2,2g -k6,6 | bedtools groupby -g 1,2,3,6 -o count -c 1 | \
		awk 'BEGIN{OFS="\t"; k=0; }{ print $1,$2,$3,++k,$5,$4; }'
}

bedtoautodist() {
	local fle=$1
	local d=${2:-50}
	local ofletag=${3:-$fle}
	ofletag=${ofletag##*/}
	ofletag=${ofletag%%.*}
	bedtools slop -s -l $d -r $d -g ~/aux/genomes/mouse.mm10.genome -i $fle | \
		bedtools intersect -s -wao -a stdin -b $fle | \
		awk -v d=$d 'BEGIN{OFS="\t"; l=2*d+1;} $4!=$10 && $3-$2==l { k=$9-$3+d; if (k<0) k=-k; array[k]+=$11;}END{ for(i in array) print i, array[i];} ' | \
		sort -k 1,1g > $ofletag.autodist
	Rscript -e 'source("~/Projects/biter_biter_shared/bin/auxggplot2.R"); ggplot <- set.ggplot();' \
					-e 'f <- "'$ofletag.autodist'"; a <- read.table(f); a <- a[order(a$V1),];' \
					-e 'ggplot(a, aes(x=V1, y=V2)) + geom_point() + geom_line(size=rel(1)) + scale_y_log10() + ' \
					-e 'labs(x="Distance", y="Empirical Distribution", title=paste0("'$ofletag'", ".autodist"));' \
					-e 'ggsave(paste0(f, ".pdf"));' \
					-e 'ggplot(a, aes(x=V1, y=cumsum(V2))) + geom_point() + geom_line(size=rel(1)) + scale_y_log10() + ' \
					-e 'labs(x="Distance", y="Empirical Cumulative Distribution", title=paste0("'$ofletag'", ".autodist"));' \
					-e 'ggsave(paste0(f, "_cs.pdf"));' 2> /dev/null
}

cluster() {
	local fle=$1
	local d=${2:-10}
	local ofletag=${3:-$fle}
	ofletag=${ofletag##*/}
	ofletag=${ofletag%%.*}
	bedtools merge -s -d $d -c 5 -o sum -i $fle | \
		awk 'BEGIN{OFS="\t"; k=0;} { print $1,$2,$3,++k,$5,$4; }' | \
		bedtools intersect -wao -s -a stdin -b $fle | sort -k 1,6 -k11,11gr -T . | \
		bedtools groupby -full -g 1,2,3,4,5,6 -c 11 -o max | \
		awk 'BEGIN{OFS="\t"; } { print $1,$2,$3,$4,$5,$6,$11,0,-1,$8-$2; }' | \
		sort -k 1,1 -k 2,2g -T . > ${ofletag}_cls.narrowPeak
	# reverse complement for DRS
	awk -v c=0 'BEGIN{OFS="\t"}$5>c{ s="+"; if ($6=="+") s="-"; print $1,$2+$NF,$2+$NF+1,$4,$5,s }' ${ofletag}_cls.narrowPeak | \
		bedtools sort -i stdin > ${ofletag}_summit.bed 
#	cluster statistics 
	Rscript -e 'source("~/Projects/biter_biter_shared/bin/auxggplot2.R"); ggplot <- set.ggplot();' \
					-e 'f <- "'${ofletag}_cls.narrowPeak'"; a <- read.table(f);' \
					-e 'ggplot(a, aes(x=V5)) + geom_density(size=rel(1)) + scale_x_log10() + ' \
					-e 'labs(x="Cluster Read Count", y="Empirical Distribution", title=paste0("'$ofletag'", "_cls"));' \
					-e 'ggsave(paste0(f, "_cls_readCount.pdf"));' \
					-e 'ggplot(a, aes(x=V3-V2)) + geom_density(size=rel(1)) + scale_x_log10() + ' \
					-e 'labs(x="Length", y="Empirical Distribution", title=paste0("'$ofletag'", "_cls"));' \
					-e 'ggsave(paste0(f, "_cls_length.pdf"));' \
					-e 'my_bin(a, aes(x=V5, y=(V7/V5)/(V3-V2))) + ' \
					-e 'labs(y="Length Normalized Ratio of Representative to Cluster Read Count", x="Cluster Read Count", title=paste0("'$ofletag'", "_cls"));' \
					-e 'ggsave(paste0(f, "_cls_representative_read_count_density_read_count.pdf"));' \
					-e 'my_bin(a, aes(x=V5, y=V7/V5)) + ' \
					-e 'labs(x="Cluster Read Count", y="Ratio of Representative to Cluster Read Count", title=paste0("'$ofletag'", "_cls"));' \
					-e 'ggsave(paste0(f, "_cls_ratio_readCount.pdf"));' \
					-e 'my_bin(a, aes(x=V3-V2, y=V7/V5)) + ' \
					-e 'labs(x="Length", y="Ratio of Representative to Cluster Read Count", title=paste0("'$ofletag'", "_cls"));' \
					-e 'ggsave(paste0(f, "_cls_length_read_count.pdf"));' 2> /dev/null

}
modelGMMcdens() {
	local fletrain=$1
	local fletest=$2
	Rscript -e 'source("~/Projects/biter_biter_shared/bin/auxggplot2.R"); ggplot <- set.ggplot();' \
	Rscript -e 'source("~/Projects/biter_biter_shared/bin/auxML.R");' \
					-e 'ggsave(paste0(f, "_cls_length_read_count.pdf"));' 2> /dev/null
}

PAS() {
	# FUNCTIONS
	export -f nucdistr nuc SVD4nucComp closest bamsto5pbed bedtoautodist cluster modelGMMcdens	
	
	indir=${o[indir]:-RM}
	ln -sf ${o[datadir]}/$indir .
	files=($(ls $indir/*.sorted.bam))

	hexamer=AATAAA
	hexamerconj=TTTTTT
	repeatMasker=~/PI_HOME/Data/casco/UCSC_tracks/mm10/rmskJoinedBaseline.bed.gz
	[ -e ${repeatMasker##*/} ] || ln -sf $repeatMasker .
	#TESfiles=(~/PI_HOME/Data/casco/UCSC_tracks/mm10/distal/TES.bed.gz ~/PI_HOME/Data/casco/UCSC_tracks/mm10/proximal/TES.bed.gz)
	# TODO prepare
	TESfiles=(~/PI_HOME/Data/casco/UCSC_tracks/mm10/mm10_refGene_distal.gtf.gz ~/PI_HOME/Data/casco/UCSC_tracks/mm10/mm10_refGene_proximal.gtf.gz)

	for TES in ${TESfiles[@]}; do
		echo MAIN TASK 1 features of $TES
		polyAfile=${TES%/*}
		polyAfile=${polyAfile##*/}polyA.bed.gz
		[[ -e $polyAfile ]] && break
		ln -sf $TES $polyAfile
		ofle=${polyAfile##*/}
		ofle=${ofle%%.*}5p
		echo TASK stop codons of of trx from $polyAfile file
		nucdistr $polyAfile -1997 -1000 3 stdout | grep -v "^>" | sort | uniq -c | sort -k 1,1gr | head
		echo TASK nucleotide distribution around PASs from $polyAfile file
		nucdistr $polyAfile -1900 -900 200 $ofle 
		nuc $polyAfile -1950 -1000 50 $hexamer ${ofle}_5p50nc 
		nuc $polyAfile -2000 -950 50 $hexamerconj ${ofle}_3p50nc 
		SVD4nucComp ${ofle}_5p50nc.nuc ${ofle}_3p50nc.nuc ${ofle}_nuc 1:2 ${o[feaPat]:-num_[AT]}
		closest $polyAfile $repeatMasker -1999 -1000 ${ofle}_repeat
	done &> XMTASK1.log

	mergedist=${o[mergedist]:-10} # decided based on the bedtoautodist profile
	for dummy in 1; do
		ofle=pooled5p
		[[ -e ${ofle}.5p.bed ]] && break
		echo MAIN TASK 2 features of $ofle
		# TODO delete???
		filetags=(${files[@]##*/})
		filetags=(${filetags[@]%%.*})
		echo TASK bamsto5pbed
		bamsto5pbed $indir ${o[nc]} ${o[extractSoftClipped]:-NULL} > $ofle.5p.bed
		for f in $indir/*sorted.bam; do fn=${f##*/}; bamsto5pbed $indir ${o[nc]} ${o[extractSoftClipped]:-NULL} ${fn%.*} > ${fn%%.*}.5p.bed; done
		echo TASK autodistance 
		for d in 50 1000; do bedtoautodist $ofle.5p.bed $d ${ofle}_within${d}nc; done
		echo TASK cluster autodistance 
		cluster $ofle.5p.bed $mergedist ${ofle}_dist$mergedist
	done &> XMTASK2.log

	# TODO debug
	ofle=${ofle}_dist${mergedist}_summit
	for dummy in 1; do
		[[ -e ${ofle}_repeat.closest ]] && break
		echo MAIN TASK PAS modeling with GMM
		nucdistr $ofle.bed 100 99 200 $ofle 
		echo TASK nucleotide distribution around PASs from $ofle.bed file
		nuc $ofle.bed 49 0 50 $hexamer ${ofle}_5p50nc 
		nuc $ofle.bed 0 49 50 $hexamerconj ${ofle}_3p50nc 
		SVD4nucComp ${ofle}_5p50nc.nuc ${ofle}_3p50nc.nuc ${ofle}_nuc 2:3 ${o[feaPat]:-num_[AT]}
		closest $ofle.bed $repeatMasker 0 0 ${ofle}_repeat
	done &> XMTASK3.log

	fleTrain=distalpolyA5p_nuc.txt # TODO set?
	fleTest=pooled5p_dist${mergedist}_summit_nuc.txt
	scoreColumn=$(less $fleTest | awk 'NR==1{ print $5;}') 
	echo ${o[feaPat]:-num_[AT]}
	#for readCountLowerLimit in 5 4 3 2 1; do
	for readCountLowerLimit in 5; do
		ofleTag=pooled5p_dist${mergedist}_rcll${readCountLowerLimit}_summit_nuc
		Rscript -e 'source("~/Projects/biter_biter_shared/bin/auxML.R"); ' \
						-e 'GMMcdens("'$fleTrain'", "'$fleTest'", "'$ofleTag'", param.preprocess=list(subset.cond="\''$scoreColumn\''>'$readCountLowerLimit'",col.sel.pat="'${o[feaPat]:-num_[AT]}'",'$readCountLowerLimit'))'
		awk 'BEGIN{OFS="\t"} $NF==2{ print $1,$2,$3,$4,$5,$6 }' ${ofleTag}_GMMcdens.txt > ${ofleTag}_selected.bed
		awk 'BEGIN{OFS="\t"} $NF==1{ print $1,$2,$3,$4,$5,$6 }' ${ofleTag}_GMMcdens.txt > ${ofleTag}_filtered.bed
		nucdistr ${ofleTag}_selected.bed 51 99 200 ${ofleTag}_selected
		nucdistr ${ofleTag}_filtered.bed 51 99 200 ${ofleTag}_filtered
	done &> MTASK4.log


	return

	echo annotate w repeat masker
#	bedtools sort -i $repeatMasker | bedtools closest -d -s -t first -a $ofle.5p.${mergedist}cls.summit.bed -b stdin > $ofle.closestRepeat
#	Rscript	-e 'source("~/Projects/biter_biter_shared/bin/auxggplot2.R"); ggplot <- set.ggplot();' \
#		-e 'a <- read.table("pooled5p.closestRepeat", header=F, comment.char=""); a$type= "pooled5p"; b <- read.table("polyA5p.closestRepeat", header=F, comment.char=""); b$type <- "polyA5p";' \
#		-e 'dta <- rbind(a,b); dta$ps <- 1; ggplot(dta, aes(color=type, x=log2(V13+ps))) + geom_density(size=rel(1)); ggsave("closestRepeat.pdf");'
# repeat distance normalized by the maximum
#	Rscript	-e 'source("~/Projects/biter_biter_shared/bin/auxggplot2.R"); ggplot <- set.ggplot();' \
#		-e 'a <- read.table("pooled5p.closestRepeat", header=F, comment.char=""); a$type= "pooled5p"; b <- read.table("polyA5p.closestRepeat", header=F, comment.char=""); b$type <- "polyA5p";' \
#		-e 'dta <- rbind(a,b); dta$ps <- 1; ggplot(dta, aes(color=type, x=log2((V13+ps)/max(V13)))) + geom_density(size=rel(1)); ggsave("closestRepeatnormByMaxDistance.pdf");'

#	awk -v b=100 -v d=100 'BEGIN{OFS="\t"} $NF>d { print $1,$2-b,$3+b-1,$4,$5,$6; }' $ofle.closestRepeat | \
#		bedtools getfasta -name -s -bed stdin -fi ~/PI_HOME/Data/casco/UCSC_tracks/mm10/assembly/all/chr.fa | \
#		gzip -c > $ofle.5p.${mergedist}cls.summitFarFromRepeat.fa.gz
#	Rscript -e 'source("~/Projects/biter_biter_shared/bin/auxGenomicFeatures.R"); fafile <- "'$ofle.5p.${mergedist}cls.summitFarFromRepeat.fa.gz'"; dta <- caller.nucProbPerPosition(fafile, start=1, end=200)' \
#		-e 'source("~/Projects/biter_biter_shared/bin/auxggplot2.R"); ggplot <- set.ggplot(); dta.m <- melt(dta, id="position"); ggplot(dta.m, aes(x=position, y=value, color=variable)) + geom_point() + geom_line(size=rel(1)); ggsave(paste0(fafile,".pdf"))' > /dev/null

#	echo coverage
#	echo "chr start end name score strand ${filetags[@]}" | sed 's/ /\t/g' > $ofle.multicov
#	awk -v c=1 'BEGIN{OFS="\t"}$5>c{ print $0 }' $ofle.5p.${mergedist}cls.narrowPeak | \
#		bedtools map -s -a stdin -b ${filetags[0]}.5p.bed -c 5 -o sum | \
#		bedtools map -s -a stdin -b ${filetags[1]}.5p.bed -c 5 -o sum | \
#		awk 'BEGIN{OFS="\t"; }{ if ($(NF-1)==".") { $(NF-1)=0; } if ($NF==".") { $NF=0; } print $0; }' >> $ofle.multicov

#	bedtools slop -s -l $lengthPA2hexamer -r 0 -i $ofle.5p.${mergedist}cls.summit.bed -g ~/aux/genomes/mouse.mm10.genome | \
#		bedtools nuc -C -fullHeader -pattern $hexamer -s -fi ~/PI_HOME/Data/casco/HISAT2_indexes/mm10/all/chr.fa -bed stdin | \
#		sed 's/^#//g' > $ofle.nuc
#	bedtools slop -s -l 0 -r $lengthPA2hexamerconj -i $ofle.5p.${mergedist}cls.summit.bed -g ~/aux/genomes/mouse.mm10.genome | \
#		bedtools nuc -C -fullHeader -pattern $hexamerconj -s -fi ~/PI_HOME/Data/casco/HISAT2_indexes/mm10/all/chr.fa -bed stdin | \
#		sed 's/^#//g' > $ofle.3p.nuc
#	Rscript -e 'source("~/Projects/biter_biter_shared/bin/auxggplot2.R"); ggplot <- set.ggplot();' \
#		-e 'a <- read.table("pooled5p.nuc", header=T); a$type= "pooled5p.5p50nc"; b <- read.table("polyA5p.5p50nc.nuc", header=T); b$type <- "polyA5p.5p50nc"; dta <- rbind(a,b); library(ggplot2); ggplot(dta, aes(color=type, x=(X9_num_A)/X15_seq_len)) + geom_density(size=rel(1)); ggsave(paste0("5p.nuc.pdf"))' 
#	Rscript -e 'source("~/Projects/biter_biter_shared/bin/auxggplot2.R"); ggplot <- set.ggplot();' \
#		-e 'a <- read.table("pooled5p.3p.nuc", header=T); a$type= "pooled5p.3p50nc"; b <- read.table("polyA5p.3p50nc.nuc", header=T); b$type <- "polyA5p.3p50nc"; dta <- rbind(a,b); library(ggplot2); ggplot(dta, aes(color=type, x=(X12_num_T)/X15_seq_len)) + geom_density(size=rel(1)); ggsave(paste0("3p.nuc.pdf"))' 

	# TODO filter internal priming sites using the code at PAS_20161118 in auxML.R
	l=$((100-$lengthPA2hexamer))
	for class in 1 2; do
		echo $ofle.5p.${mergedist}cls.summitclass$class.fa.gz
		awk -v class=$class 'BEGIN{OFS="\t"} NR>1 && $NF==class { print $2,$3,$4,$1,$5,$6 }' test.txt | \
			bedtools slop -s -l $l -r 99 -i stdin -g ~/aux/genomes/mouse.mm10.genome | \
			awk '$3-$2==200 { print }' | \
			bedtools getfasta -name -s -bed stdin -fi ~/PI_HOME/Data/casco/UCSC_tracks/mm10/assembly/all/chr.fa | \
			gzip -c > $ofle.5p.${mergedist}cls.summitclass$class.fa.gz
		Rscript -e 'source("~/Projects/biter_biter_shared/bin/auxGenomicFeatures.R"); fafile <- "'$ofle.5p.${mergedist}cls.summitclass$class.fa.gz'"; dta <- caller.nucProbPerPosition(fafile, start=1, end=200)' \
			-e 'source("~/Projects/biter_biter_shared/bin/auxggplot2.R"); ggplot <- set.ggplot(); dta.m <- melt(dta, id="position"); ggplot(dta.m, aes(x=position, y=value, color=variable)) + geom_point() + geom_line(size=rel(1)); ggsave(paste0(fafile,".pdf"))' > /dev/null
	done
	return;


	# TODO continue here
	echo nucleotide distribution
	bedtools slop -s -l 100 -r 99 -i $ofle.5p.${mergedist}cls.summit.bed -g ~/aux/genomes/mouse.mm10.genome | \
		awk '$3-$2==200 { print }' | \
		bedtools getfasta -name -s -bed stdin -fi ~/PI_HOME/Data/casco/UCSC_tracks/mm10/assembly/all/chr.fa | \
		gzip -c > $ofle.5p.${mergedist}cls.summit.fa.gz
	Rscript -e 'source("~/Projects/biter_biter_shared/bin/auxGenomicFeatures.R"); fafile <- "'$ofle.5p.${mergedist}cls.summit.fa.gz'"; dta <- caller.nucProbPerPosition(fafile, start=1, end=200)' \
		-e 'source("~/Projects/biter_biter_shared/bin/auxggplot2.R"); ggplot <- set.ggplot(); dta.m <- melt(dta, id="position"); ggplot(dta.m, aes(x=position, y=value, color=variable)) + geom_point() + geom_line(size=rel(1)); ggsave(paste0(fafile,".pdf"))' > /dev/null


	l=$((100-$lengthPA2hexamer))
#	awk -v ATpct=0.7 'BEGIN{OFS="\t"} NR>1 && $7>ATpct { print $1,$2,$3,$4,$5,$6 }' $ofle.nuc | \
#		bedtools slop -s -l $l -r 99 -i stdin -g ~/aux/genomes/mouse.mm10.genome | \
#		awk '$3-$2==200 { print }' | \
#		bedtools getfasta -name -s -bed stdin -fi ~/PI_HOME/Data/casco/UCSC_tracks/mm10/assembly/all/chr.fa | \
#		gzip -c > $ofle.5p.${mergedist}cls.summitATrich.fa.gz
#	Rscript -e 'source("~/Projects/biter_biter_shared/bin/auxGenomicFeatures.R"); fafile <- "'$ofle.5p.${mergedist}cls.summitATrich.fa.gz'"; dta <- caller.nucProbPerPosition(fafile, start=1, end=200)' \
#		-e 'source("~/Projects/biter_biter_shared/bin/auxggplot2.R"); ggplot <- set.ggplot(); dta.m <- melt(dta, id="position"); ggplot(dta.m, aes(x=position, y=value, color=variable)) + geom_point() + geom_line(size=rel(1)); ggsave(paste0(fafile,".pdf"))' > /dev/null
	awk -v Apct=0.2 'BEGIN{OFS="\t"} NR>1 && $9/$15>Apct { print $1,$2,$3,$4,$5,$6 }' $ofle.nuc | \
		bedtools slop -s -l $l -r 99 -i stdin -g ~/aux/genomes/mouse.mm10.genome | \
		awk '$3-$2==200 { print }' | \
		bedtools getfasta -name -s -bed stdin -fi ~/PI_HOME/Data/casco/UCSC_tracks/mm10/assembly/all/chr.fa | \
		gzip -c > $ofle.5p.${mergedist}cls.summitArich.fa.gz
	Rscript -e 'source("~/Projects/biter_biter_shared/bin/auxGenomicFeatures.R"); fafile <- "'$ofle.5p.${mergedist}cls.summitArich.fa.gz'"; dta <- caller.nucProbPerPosition(fafile, start=1, end=200)' \
		-e 'source("~/Projects/biter_biter_shared/bin/auxggplot2.R"); ggplot <- set.ggplot(); dta.m <- melt(dta, id="position"); ggplot(dta.m, aes(x=position, y=value, color=variable)) + geom_point() + geom_line(size=rel(1)); ggsave(paste0(fafile,".pdf"))' > /dev/null
	return	

	echo annotate w known polyA sites
	bedtools sort -i $polyAfile | bedtools closest -d -s -t first -a $ofle.5p.${mergedist}cls.summit.bed -b stdin > $ofle.closestPolyA
	awk -v b=100 -v d=1 'BEGIN{OFS="\t"} $NF>-1 && $NF<d { print $1,$2-b,$3+b-1,$4,$5,$6; }' $ofle.closestPolyA  | \
		bedtools getfasta -name -s -bed stdin -fi ~/PI_HOME/Data/casco/UCSC_tracks/mm10/assembly/all/chr.fa | \
		gzip -c > $ofle.5p.${mergedist}cls.summitclosest2PA.fa.gz
	Rscript -e 'source("~/Projects/biter_biter_shared/bin/auxGenomicFeatures.R"); fafile <- "'$ofle.5p.${mergedist}cls.summitclosest2PA.fa.gz'"; dta <- caller.nucProbPerPosition(fafile, start=1, end=200)' \
		-e 'source("~/Projects/biter_biter_shared/bin/auxggplot2.R"); ggplot <- set.ggplot(); dta.m <- melt(dta, id="position"); ggplot(dta.m, aes(x=position, y=value, color=variable)) + geom_point() + geom_line(size=rel(1)); ggsave(paste0(fafile,".pdf"))' > /dev/null
	less $ofle.closestPolyA | awk -v d=1 '$NF>-1 && $NF<d { print }' | sort -k 10,10 -k 5,5gr -T . | bedtools groupby -g 10 -full -o max -c 5 | \
		awk -v b=100 'BEGIN{OFS="\t"} { print $1,$2-b,$3+b-1,$4,$5,$6; }' | \
		bedtools getfasta -name -s -bed stdin -fi ~/PI_HOME/Data/casco/UCSC_tracks/mm10/assembly/all/chr.fa | \
		gzip -c > $ofle.5p.${mergedist}cls.repsummitclosest2PA.fa.gz
	Rscript -e 'source("~/Projects/biter_biter_shared/bin/auxGenomicFeatures.R"); fafile <- "'$ofle.5p.${mergedist}cls.repsummitclosest2PA.fa.gz'"; dta <- caller.nucProbPerPosition(fafile, start=1, end=200)' \
		-e 'source("~/Projects/biter_biter_shared/bin/auxggplot2.R"); ggplot <- set.ggplot(); dta.m <- melt(dta, id="position"); ggplot(dta.m, aes(x=position, y=value, color=variable)) + geom_point() + geom_line(size=rel(1)); ggsave(paste0(fafile,".pdf"))' > /dev/null


#	less $ofle.closestPolyA | sort -k 11,11 -T . | bedtools groupby -g 11 -c 13,13,13,13 -o count,mean,sstdev,median

#	samtools view -b ${files[0]} chr10 | bedtools bamtobed -cigar -i stdin | bedtools sort -i stdin | \
#		bedtools closest -d -S -t first -a stdin -b ${polyAfile##*/} > ${files[0]##*/}.closestPolyA

	# check code in auxggplot2 PAS_20161028

	echo DONE
}

#slurmMain PAS --datadir=~/Projects/biter_demorree_PAS/data/deMorree_2012/directRNAseq --t=2 --nc=2 --dryrun=T --skippostproc=T --otag=TAGME \
#	--jobnames=genomecov
#slurmMain PAS23 --datadir=~/Projects/biter_demorree_PAS/data/deMorree_2012/directRNAseq --indir=RM_TR_trim5termT_L06minscoreNSC --t=2 --nc=2 --dryrun=T --skippostproc=T --otag=md15_numAT5p --feaPat="num_[AT].*nuc5p"
#slurmMain PAS23 --datadir=~/Projects/biter_demorree_PAS/data/deMorree_2012/directRNAseq --indir=RM_TR_trim5termT_L06minscoreNSC --t=2 --nc=2 --dryrun=T --skippostproc=T --otag=md15_numAT --feaPat="num_[AT]"
#slurmMain PAS23 --datadir=~/Projects/biter_demorree_PAS/data/deMorree_2012/directRNAseq --indir=RM_TR_trim5termT_L06minscoreNSC --t=2 --nc=2 --dryrun=T --skippostproc=T --otag=md10_numAT --feaPat="num_[AT]" --mergedist=10
#slurmMain PAS23 --datadir=~/Projects/biter_demorree_PAS/data/deMorree_2012/directRNAseq --indir=RM_TR_trim5termT_L06minscoreNSC --t=2 --nc=2 --dryrun=T --skippostproc=T --otag=md15
#slurmMain PAS23 --datadir=~/Projects/biter_demorree_PAS/data/deMorree_2012/directRNAseq --indir=RM_TR_trim5termT_L04minscoreNSC --t=2 --nc=2 --dryrun=T --skippostproc=T --otag=md15
#slurmMain PAS2 --datadir=~/Projects/biter_demorree_PAS/data/deMorree_2012/directRNAseq --indir=RM_TR_trim5termT_L06minscore --t=2 --nc=2 --dryrun=T --skippostproc=T 
slurmMain PAS --datadir=~/Projects/biter_demorree_PAS/data/deMorree_2012/directRNAseq --indir=RM_TR_trim5termT_L06minscore --t=2 --nc=2 --dryrun=T --skippostproc=T --otag=md10_numAT_softClipped --feaPat="num_[AT]" --mergedist=10 --extractSoftClipped=T



