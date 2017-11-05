source("~/Projects/biter_biter_shared/bin/bigData.R")
library(GenomicRanges)
library(rtracklayer)
library(dplyr)

test.selectTrxAddGeneSymbol2gtf <- 
	function() { 
	source("~/Projects/biter_biter_shared/bin/auxRefGene.R")
	# prep canonical's gtf with gene symbols 
	gtfFile <- "mm10_refGene.gtf.gz"
	seltabFile <- "mm10_refGene.seltab.gz"
	type <- "canonical"; version <- ".v2"
	outFile <- paste0("mm10_refGene_",type,version,".gtf.gz")
	selectTrxAddGeneSymbol2gtf(gtfFile, seltabFile, outFile, type, verbose=F)

	# for PAS analysis of Antoine
	seltabFile <- "mm10_refGene.seltab.gz"
	type <- "canonical"; version <- ".v2"
	outFile <- paste0("mm10_refGene_",type,version,".gtf.gz")
	selectTrxAddGeneSymbol2gtf(gtfFile, seltabFile, outFile, type, verbose=F)
}

# seltab file should have the following ordered information for type=canonical 
#	"refSeqId","geneId","description","kcClusterId","kcTranscriptId","ensemblId","status","bioType"     
# seltab file should have the following ordered information for type=proximalPAS
#	"refSeqId","geneId"
# canonical transcripts are already mRNAs and long non-coding RNAs and mRNAs!
selectTrxAddGeneSymbol2gtf <- 
	function(gtfFile="mm10_refGene.gtf.gz", seltabFile="mm10_refGene.seltab.gz",
					 outFile="mm10_refGene_canonical.gtf.gz", type="canonical", verbose=F, ...) { 
	# 20161129 duplicates removed with select and distinct		
	selected <- readBigTable(seltabFile,sep="\t", quote="",
		col.names=c("transcript_id", "geneId", "description", "kcClusterId", "kcTranscriptId", "ensemblId", "status", "bioType"), verbose=verbose) %>% dplyr::select(transcript_id, geneId, kcTranscriptId) %>% arrange %>% distinct
	dim(selected)

	# select by alphabetic order the first if there are >1 canonical transcripts per gene
	selected <- selected %>% filter(kcTranscriptId!="n/a") %>% 
		arrange(geneId, transcript_id) %>% group_by(geneId) %>% slice(1)

	# merge with selected by transcript_id
  selected.gtf <- import.gff(gtfFile) %>% as.data.frame %>% inner_join(selected) 
	selected.gtf %>% count

	# select first chr if a transcript_id maps to multiple chromosomes
	selected.gtf <- dplyr::select(selected.gtf, geneId,transcript_id,seqnames) %>%
		arrange %>% distinct %>% group_by(geneId,transcript_id) %>% slice(1) %>%
		ungroup %>% inner_join(selected.gtf)
	selected.gtf %>% count
	selected.gtf %>% head

	# TODO 
	trxPAS <- filter(selected.gtf, type=="exon") %>% group_by(geneId,transcript_id) %>%
		mutate(PAS=ifelse(strand=="+",end,start)) %>% arrange(geneId,transcript_id,PAS) %>%
		slice(ifelse(strand=="+",n(),1)) %>% ungroup %>% distinct %>% group_by(geneId)
	trxPAS %>% dim

	proximal <- slice(trxPAS, ifelse(strand=="+",1,n())) %>% ungroup %>% distinct
	distal   <- slice(trxPAS, ifelse(strand=="+",n(),1)) %>% ungroup %>% distinct
	proximaldistal <- inner_join(proximal, distal, by="geneId", suffix=c(".x", ".y")) %>% filter(PAS.x==PAS.y)
	proximaldistal %>% count

	# revise by set difference 
	proximal <- anti_join(proximal, proximaldistal, by=c("transcript_id"="transcript_id.x")) %>% mutate(PAStype="proximal")
	distal <- anti_join(distal, proximaldistal, by=c("transcript_id"="transcript_id.y")) %>% mutate(PAStype="distal") 
	proximaldistal <- dplyr::select(proximaldistal, -matches(".y$")) %>% mutate(PAStype="proximaldistal")
	names(proximaldistal) <- names(proximal) # get rid of .x
	proximal %>% count
	distal %>% count
	proximaldistal %>% count

	selected.gtf <- inner_join(selected.gtf, 
														 select(get(type), PAStype,transcript_id), 
														 by="transcript_id")
	# add attributes 
	selected.gtf <- mutate(selected.gtf, 
												 attribute=paste0('gene_id "',geneId,'"; transcript_id "',transcript_id, '";'))

	if (type == "proximaldistal")
		selected.gtf <- mutate(selected.gtf, attribute=paste0(attribute, ' PAStype "', PAStype, '";'))
	 
	#write as file
	writeData(dplyr::select(selected.gtf, seqnames,source,type,start,end,score,strand,phase,attribute),
						ofle=outFile, row.names=F, col.names=F, na=".")

}

