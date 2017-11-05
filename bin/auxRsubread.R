source("~/Projects/biter_biter_shared/bin/bigData.R")

# TODO write me
caller.annotations <-
	function(gr.a, gr.b) {


}

test.caller.featureCounts <- 
	function() {
	source("~/Projects/biter_lingliu_DNAdamage/bin/auxRsubread.R")
	bamdir <- "~/Projects/biter_lingliu_DNAdamage/results/2015-05-18/HISAT_mapOut_Hr_cKO/"
	# param list for unstranded paired RNAseq
	param.featureCounts <- list(annot.inbuilt="mm10", annot.ext="mm10_refGene_canonical.gtf",  isGTFAnnotationFile=T, GTF.featureType="exon", GTF.attrType="transcript_id", isPairedEnd=T, requireBothEndsMapped=T, checkFragLength=T, maxFragLength=50000, strandSpecific=0, read2pos=NULL)
	ofleTag <- ""
	mc.cores <- 1
	fc <- caller.featureCounts(bamdir, param.featureCounts=param.featureCounts, ofleTag=ofleTag, mc.cores=mc.cores)

}

# TODO write me
# use it for 3' end sequencing bias; plot ecdf of TSS and TES read density in different RNAseq libraries 
# use it for antibody testing in ChIP-seq
# atac-seq peaks
caller.featureCounts <- 
	function(bamdir, pattern=".sorted.bam$", param.featureCounts=NULL, ofleTag=NULL, mc.cores=1) {
	library(Rsubread)
	files <- dir(bamdir, pattern=pattern, full.names=T)
	files
# TODO use these to get statistics; propmapped is slightly higher than HISAT reported mapping rates
# atcgContent
# qualityScores
#lapply(1:length(files), function(i) { print(propmapped(output.files[i]))})

	fc.SR <- do.call("featureCounts", c(list(files=files, nthreads=mc.cores), param.featureCounts))
	dta <- cbind(fc.SR$counts, fc.SR$annotation) 
	dta <- dta[c("GeneID", setdiff(names(dta), "GeneID"))]

	if (!is.null(ofleTag)) {
		commentLine <- paste0(lapply(attributes(param.featureCounts), function(x) paste(x, param.featureCounts[x], sep="=")), collapse=" ")
		
		# write stats
		ofle <- paste(ofleTag, "featureCounts", ".stat", sep="")
		writeData(fc.SR$stat, ofle=ofle, commentLine=commentLine, row.names=F)

		# write count data
		ofle <- paste(ofleTag, ".featureCounts.gz", sep="")
		writeData(dta, ofle=ofle, commentLine=commentLine, row.names=F)
	}
	dta
}
