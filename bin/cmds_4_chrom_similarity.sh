#!/bin/bash -l 
## sbatch command: sbatch -J FC_Cheung_2013_RNAseq -p normal -t 6:0:0 -c 16 			-e slurm_FC_Cheung_2013_RNAseq.out -o slurm_FC_Cheung_2013_RNAseq.out bin/sbatch_FC_Cheung_2013_RNAseq.sh
# sbatch command: sbatch -J featureCounts_genome -p normal -t 6:0:0 -c 8 -e slurm_featureCounts_genome.out -o slurm_featureCounts_genome.out cmds_4_chrom_similarity.sh  

indir=../RM_TR_ada_Liu_2015_ATACseq # for 20160113_4TomRandoGrant in 2015-12-09_chromatin
cutoff=10
tag=_FC_wlongFragments
tag=_FC_wlongFragments_cutoff20
tag=_FC_wlongFragments_cutoff5
tag=_FCaging_wlongFragments_cutoff$cutoff
mkdir -p $tag; pushd $tag


less ~/aux/genomes/mouse.mm10.genome.txt  | grep -w chr10 | awk '{ print $1"\t"$2}' > chr10.txt
bedtools makewindows -g chr10.txt -w 300 -s 150 | \
	awk 'BEGIN{k=0; OFS="\t"; } { k=k+1; print $1,"custom","exon",($2+1),$3,1000,".",".","gene_id \""k"\"; transcript_id \""k"\""; }' | \
	less > mm10.gtf 

## paired-end sequencing overlaps not allowed
##featureCounts -T 16 -t exon -g transcript_id -p -B -C 		-s 0 -a mm10.gtf -o featureCounts.txt RM/02m_05hdPI_SC_1.sorted.bam RM/02m_05hdPI_SC_2.sorted.bam RM/02m_NA_SC_1.sorted.bam RM/02m_NA_SC_2.sorted.bam;
## paired-end sequencing overlaps allowed
##featureCounts -T 16 -t exon -g transcript_id -p -B -C -s 0 -O -a mm10.gtf -o featureCounts.txt RM_Liu_2015_ATACseq/WT_old_NA_SC_1.sorted.bam  RM_Liu_2015_ATACseq/WT_young_60hPI_SC_1.sorted.bam  RM_Liu_2015_ATACseq/WT_young_NA_SC_1.sorted.bam RM_Liu_2015_ATACseq/WT_old_NA_SC_2.sorted.bam  RM_Liu_2015_ATACseq/WT_young_60hPI_SC_2.sorted.bam  RM_Liu_2015_ATACseq/WT_young_NA_SC_2.sorted.bam
##featureCounts -T 16 -t exon -g transcript_id -p -B -C -D 2000 -s 0 -O -a mm10.gtf -o featureCounts.txt ../RM_Liu_2015_ATACseq/WT_old_NA_SC_1.sorted.bam  ../RM_Liu_2015_ATACseq/WT_young_60hPI_SC_1.sorted.bam  ../RM_Liu_2015_ATACseq/WT_young_NA_SC_1.sorted.bam ../RM_Liu_2015_ATACseq/WT_old_NA_SC_2.sorted.bam  ../RM_Liu_2015_ATACseq/WT_young_60hPI_SC_2.sorted.bam  ../RM_Liu_2015_ATACseq/WT_young_NA_SC_2.sorted.bam
ln -sf $indir input
featureCounts -T 16 -t exon -g transcript_id -p -B -C -D 2000 -s 0 -O -a mm10.gtf -o featureCounts.txt input/WT_old_NA_SC_1.srt.bam  input/WT_old_NA_SC_2.srt.bam input/WT_young_NA_SC_1.srt.bam input/WT_young_NA_SC_2.srt.bam 
## single-end sequencing overlaps allowed
##featureCounts -T 16 -t exon -g transcript_id         -s 0 -O -a mm10.gtf -o featureCounts.txt RM_Liu_2015_ATACseq/WT_old_NA_SC_1.sorted.bam  RM_Liu_2015_ATACseq/WT_young_60hPI_SC_1.sorted.bam  RM_Liu_2015_ATACseq/WT_young_NA_SC_1.sorted.bam RM_Liu_2015_ATACseq/WT_old_NA_SC_2.sorted.bam  RM_Liu_2015_ATACseq/WT_young_60hPI_SC_2.sorted.bam  RM_Liu_2015_ATACseq/WT_young_NA_SC_2.sorted.bam

featureCounts -v;
gzip featureCounts.txt;
rm -rf mm10.gtf;

echo DONE

Rscript 			-e "source('~/Projects/biter_biter_shared/bin/auxML.R')" 			-e "param.preprocess <- list(col.sel.pat='.srt.bam$', colname.split.rank=2, colname.split.pat='.')" 			-e "spike <- list(spikeId='ERCC', removeSpike=F)" 			-e "param.process <- list(method='filterRowMeanMin', cutoff=$cutoff)" -e "dta.exp <- caller.getExpressed('featureCounts.txt.gz', param.preprocess=param.preprocess, ofle.tag='featureCounts', nc=8, spike=spike, row.names='Geneid', param.process=param.process)"
#Rscript 			-e "source('~/Projects/biter_biter_shared/bin/auxML.R')" 			-e "param.preprocess <- list(col.sel.pat='.srt.bam$', colname.split.rank=2, colname.split.pat='.')" 			-e "spike <- list(spikeId='ERCC', removeSpike=F)" 			-e "param.process <- list(method='filterRowMeanMin', cutoff=$cutoff)" -e "dta.exp <- caller.getExpressed('../_FC_wlongFragments_cutoff5/featureCounts.txt.gz', param.preprocess=param.preprocess, ofle.tag='featureCounts', nc=8, spike=spike, row.names='Geneid', param.process=param.process)"

echo DONE2
popd
