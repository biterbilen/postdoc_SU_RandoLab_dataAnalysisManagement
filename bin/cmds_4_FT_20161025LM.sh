#samtools view -b /home/biter/Projects/biter_biter_stemCellFateRegulation/data/Liu_2015/ATACseq/RM_TR/WT_old_NA_SC_1.sorted.bam chr19 > chr19.bam
#samtools view -b ~/Projects/biter_lingliu_DNAdamage/data/Liu_2014/MNaseseq/RM_TR_ada/WT_young_NA_SC_1.sorted.bam chr19 > chr19_singleEnd.bam
#echo "which pooled data classifier gives the highes number of sites for WT_young_NA_SC_2, the smallest library?" 

# USAGE
#sh cmds_4_FT_20161025LM.sh pairedEndwminus_WT &> log.pairedEndwminus_WT

tag=${1:-pairedEndwminus_WT};
tag=${1:-pairedEndwminus_WT};
odir=PIQ4LM.$tag

mkdir -p $odir; pushd $odir
for type in pdf bed; do
	for f in ../$tag/*.calls*/*$type; do
	#for f in $tag*WT/*.calls_WT_young_60hPI_SC_2/*$type; do
	#for f in $tag*/*.calls/*$type; do
		ftag=${f#*-}
		ftag2=${ftag/.calls/}
		ftag3=${ftag2/\//_}
		[ $type == pdf ] && cp $f $ftag3
		#[ $type == bed ] && awk -F "\t" '$5>700{ print }' $f | sort -k 5,5gr | gzip -c > $odir/${ftag3/-calls.all.bed/.bed.gz}
		[ $type == bed ] && awk -F "\t" '$5>700{ print }' $f | sort -k 1,1 -k 2,2g > ${ftag3/-calls.all.bed/.bed}
	done
done

# ----------------
[ ! -e promoters.bed ] && less /home/biter/PI_HOME/Data/casco/UCSC_tracks/mm10/canonical/promoters.bed.gz | sort -k 1,1 -k 2,2g > promoters.bed
[ ! -e mm10_refGene.seltab_entrez.gz ] && ln -sf /home/biter/PI_HOME/Data/casco/UCSC_tracks/mm10/mm10_refGene.seltab_entrez.gz .

# ----------------
echo TASK1 fisher stats of reproducibility of replicates, not tested or used in LM
sort -k 1,1 -k 2,2g /home/biter/aux/genomes/mouse.mm10.genome > genome.mm10
bedtools fisher -a WT_WT_young_60hPI_SC_1_99-MA00991Fos.bed -b WT_WT_young_60hPI_SC_2_99-MA00991Fos.bed -g genome.mm10
bedtools fisher -a WT_WT_young_NA_SC_1_99-MA00991Fos.bed -b WT_WT_young_NA_SC_2_99-MA00991Fos.bed -g genome.mm10
bedtools fisher -a WT_WT_old_NA_SC_1_99-MA00991Fos.bed -b WT_WT_old_NA_SC_2_99-MA00991Fos.bed -g genome.mm10

# ----------------
echo TASK2 peak count and relative distance to the promoter regions
d=1000
for tag in .RC ""; do
f=WT_WT_old_NA_SC_1_99-MA00991Fos$tag.bed; echo $(wc -l $f) $(bedtools closest -t first -d -a $f -b promoters.bed  | cut -f 1-6,10,13 | awk -v d=$d '$8<d{print }' | wc -l) | awk '{ print $2, $1, 100*$3/$1}'
f=WT_WT_old_NA_SC_2_99-MA00991Fos$tag.bed; echo $(wc -l $f) $(bedtools closest -t first -d -a $f -b promoters.bed  | cut -f 1-6,10,13 | awk -v d=$d '$8<d{print }' | wc -l) | awk '{ print $2, $1, 100*$3/$1}'
f=WT_WT_young_NA_SC_1_99-MA00991Fos$tag.bed; echo $(wc -l $f) $(bedtools closest -t first -d -a $f -b promoters.bed  | cut -f 1-6,10,13 | awk -v d=$d '$8<d{print }' | wc -l) | awk '{ print $2, $1, 100*$3/$1}'
f=WT_WT_young_NA_SC_2_99-MA00991Fos$tag.bed; echo $(wc -l $f) $(bedtools closest -t first -d -a $f -b promoters.bed  | cut -f 1-6,10,13 | awk -v d=$d '$8<d{print }' | wc -l) | awk '{ print $2, $1, 100*$3/$1}'
f=WT_WT_young_60hPI_SC_1_99-MA00991Fos$tag.bed; echo $(wc -l $f) $(bedtools closest -t first -d -a $f -b promoters.bed  | cut -f 1-6,10,13 | awk -v d=$d '$8<d{print }' | wc -l) | awk '{ print $2, $1, 100*$3/$1}'
f=WT_WT_young_60hPI_SC_2_99-MA00991Fos$tag.bed; echo $(wc -l $f) $(bedtools closest -t first -d -a $f -b promoters.bed  | cut -f 1-6,10,13 | awk -v d=$d '$8<d{print }' | wc -l) | awk '{ print $2, $1, 100*$3/$1}'
done | sed 's/ /\t/g'

# ----------------
echo TASK3 distance to the promoter region - ASC target sites gets far away from promoter region!
fle=input4plotDistance2promoterRegion
echo variable value gid | sed 's/ /\t/g' > $fle
for f in *.RC.bed; do
	tag=${f/_99-MA00991Fos.RC.bed}
	tag2=${tag}_99-MA00991Fos
	cat $tag2.RC.bed $tag2.bed | sort -k 1,1 -k 2,2g | bedtools closest -t first -d -a - -b promoters.bed | awk -v tag=$tag 'BEGIN{OFS="\t"}{ print tag,$13, $10}' 
done >> $fle
# TODO check the plotting code in auxggplot2 #0 LM20161025
Rscript -e "source('~/Projects/biter_biter_shared/bin/auxggplot2.R'); fle <- '$fle'; get.bwplot(fle, param.preprocess=list(col.sel.pat=c('value'), ps=1, partial=T), ofle=paste0(fle,'.pdf'),panel='panel.violin', param.plot=list(scales=list(y=list(log=10))));"

# ----------------
echo TASK3 gene annotation w/peak counts

fles="bamreadcounts"
awk 'NR>1{ OFS="\t"; print $1, $3}' ~/Projects/biter_biter_stemCellFateRegulation/data/Liu_2015/ATACseq/RM_TR/summary | sed 's/RM_TR_Liu_2015_ATACseq_//g' > $fles

fle=input4peakCountPerGene
less promoters.bed | cut -f 4 | sort > $fle
tags=(Id)
index=2
for f in *.RC.bed; do
	tag=${f/_99-MA00991Fos.RC.bed}
	tags=(${tags[@]} $tag)
	cat $tag*.bed | sort -k 1,1 -k 2,2g | bedtools closest -t first -d -a - -b promoters.bed | cut -f 10 | sort | uniq -c | awk 'BEGIN{OFS="\t"}{ print $2,$1}' > tmp 
	echo "1-$((index-1)),$((index+1))"
	perl ~/PI_HOME/Applications/bilebi00/Project_EWSR1/scripts/_PAPD5/leftJoin.pl $fle tmp 1 1 "1-$((index-1)),$((index+1))" > tmp2
	mv tmp2 $fle
	index=$((index+1))
done 
mv $fle tmp
echo ${tags[@]} | sed "s/ /\t/g" > $fle
cat tmp >> $fle
awk '$3+$4+$5+$6+$7+$8>0 || NR==1{ print }' $fle | cut -f 1,3- > tmp 
mv tmp $fle
head $fle
Rscript -e 'source("~/Projects/biter_biter_shared/bin/auxggplot2.R"); fle <- "'$fle'"; fles <- "bamreadcounts"; N <- t(readBigTable(fles,row.names=1)); dta <- readBigTable(fle, header=T, row.names="Id");  dta1 <- rbind(colSums(dta),N); rownames(dta1) <- c("TFBS", "uniqueMapperFileSizeinGB"); xt <- chisq.test(dta1); xt <- chisq.test(dta1); print(list(statistic=xt$statistic, pvalue=xt$p.value, data=xt$observed) ); dtat <- apply(dta, 1, function(x) { xt <- chisq.test(rbind(x,N)); c(statistic=xt$statistic, pvalue=xt$p.value)} ); dtaa <- cbind(dta,t(dtat)); dtaa$FDR <- p.adjust(dtaa$pvalue, method = "fdr", n = length(dtaa$pvalue)); writeData(dtaa,paste0(fle, ".", fles, "normalized.chisqDE.gz"), row.names=T)'  &> yASCenrichmentStat.log
less $fle.${fles}normalized.chisqDE.gz | sort -k 10,10g | perl ~/PI_HOME/Applications/bilebi00/Project_EWSR1/scripts/_PAPD5/leftJoin_stdin.pl mm10_refGene.seltab_entrez.gz 1 1 '1-10,12' | awk -v fdr=0.01 'BEGIN{OFS="\t"} $10<fdr || NR==1{ if (NR==1) $11="gid"; print $0}' > $fle.significant 

Rscript -e 'source("~/Projects/biter_biter_shared/bin/auxggplot2.R"); fle <- "'$fle'.significant"; dta1 <- readBigTable(fle, header=T); dta1 <- dta1[,c(grep("gid", names(dta1)), grep("WT", names(dta1)))]; dta1$rsums <- rowSums(dta1[,-1]); dta1$gid <- factor(dta1$gid, levels=(dta1$gid)[order(dta1$rsums)]); dta1$rsums <- NULL; ggplot(melt(dta1), aes(y=gid, x=variable, fill=value)) + geom_tile(size=5) + scale_fill_gradientn(colours=c("gray30", "white", "gold")) + theme(axis.text=element_text(size=rel(1),lineheight=2),axis.text.x=element_text(angle=90),axis.title=element_blank()); ggsave(paste0(fle,".pdf"))' &>> yASCenrichmentStat.log

popd



