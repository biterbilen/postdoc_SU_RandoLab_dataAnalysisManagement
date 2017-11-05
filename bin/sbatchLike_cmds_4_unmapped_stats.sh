#!/bin/bash -l
nc=8; t=8

s1=ATCGGAAGAGC # truseq adaptor
s2=CGGAAGAGCACACGTCT # truseq adaptor
s3=GCACACGTCTGAACTCCA # http://www.ncbi.nlm.nih.gov/nucleotide/685047575?report=genbank&log$=nuclalign&blast_rank=1&RID=G89EZTTP01R
s4=ACGTCTGAACTCCAGTCACCTT # http://www.ncbi.nlm.nih.gov/nucleotide/690029403?report=genbank&log$=nuclalign&blast_rank=1&RID=G89Z0F3001R
s5=ACGTCTGAACTCCA
ss1=TTAGGG; rss1=CCCTAA
ss2=TTTTCA; rss2=TGCCCC
#ds=($(ls -d RM_Liu_2015_ChIPseq RM_Liu_2016_ChIPseq RM_Liu_2016_ChIPseq2 RM_TR_ada_Liu_2015_ChIPseq2 RM_*Zhang*ChIPseq RM_Liu_2016_RNAseq RM_TR_ada_vanVelthoven_2016_PullDownTrx))
ds=($(ls -d RM_de* RM_wLS_deMorree_2012_directRNAseq* RM_wLS2_deMorree_2012_directRNAseq*))


# Stats excluding adaptors
sbatch -J unmappedStatswoTSAdaptors -p normal -t $t:0:0 -c $nc <<SBATCH
#!/bin/bash -l
	for d in ${ds[@]}; do
		for f in \$d/*_un.fastq.gz \$d/*_un_[12].fastq.gz; do 
			echo \$f
			less \$f | awk 'NR % 4 == 2{ print }' | \
				grep -v $s1 | grep -v $s2 | \
				sort  -T . | uniq -c | sort -k 1,1gr -T . | head -n 100
		done > \$d.unmapped_stats.txt;
	done
SBATCH

