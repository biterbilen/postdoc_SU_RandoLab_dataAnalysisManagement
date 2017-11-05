source("~/Projects/biter_biter_shared/bin/bigData.R")
library("ENCODExplorer")

test.getENCODETreatControl <- function() {
	source("~/Projects/biter_biter_shared/bin/auxENCODExplorer.R")
	caller.getENCODETreatControl(target="H3K9me2", file_format="fastq")
}

caller.getENCODETreatControl <- function(oflet="metadata", nc=1, ...) {

	# treatment
	dtat <- getENCODEmetadata(...)

	# controls
	ff <- mergeSetwDefault(list(...), list(file_format="fastq"))$file_format
	res <- mclapply(unique(dtat$controls), function(x)   
		getENCODEmetadata(set_accession=strsplit(x, "/")[[1]][3], file_format=ff), 
		mc.cores=nc)
	dtac <- do.call(rbind, res)

	ofle <- paste(c(oflet, paste.namedList(list(...), sp=":", sp2="_")), collapse=".");
	dta <- rbind(dtat, dtac)
	writeData(dta, ofle, row.names=F, col.names=T)
	sessionInfo()
}

test.getENCODEmetadata <- function() {
	source("~/Projects/biter_biter_shared/bin/auxENCODExplorer.R")
	getENCODEmetadata(set_accession="ENCSR569IIN", file_format="fastq")
	dta <- getENCODEmetadata(target="H3K9me", file_format="fastq")
	# ghost files
	dta <- getENCODEmetadata(set_accession="ENCSR471VHW", file_format="fastq")
	# all file types
	dta <- getENCODEmetadata(set_accession="ENCSR471VHW")["biosample_name"]
	dta <- getENCODEmetadata(set_accession="ENCSR000COJ", file_format="bed")
	(dta)
#dta <- queryEncode(set_accession="ENCSR471VHW", file_format="fastq")$experiment
	dta <- getENCODEmetadata(assay="ChIP-seq")
	dta <- dta[!is.na(dta$controls), c("accession", "biosample_name", "target", "assay")] 
	writeData(dta, "dta.txt", row.names=T)
	head(dta)
	dim(dta)

	dta1 <- dta[order(dta[,2], dta[,3]),]
	dim(dta1)
	writeData(dta1, "dta1.txt", row.names=T)
	head(dta1)

	nrow(dta1)
#sessionInfo()
	# all mouse accessibility
	dta <- getENCODEmetadata(organism="Mus musculus", assay="DNase-seq+ATAC-seq")
	dta <- dta[, c("accession", "biosample_name", "target", "assay")] 
	dta
	writeData(dta, "dta_DNAaccessibility.txt", row.names=T)

	# C2C12 all mouse DNA methylation 
	dta <- getENCODEmetadata(file_format="bed", biosample="C2C12")
	dta <- dta[, c("accession", "biosample_name", "target", "assay")] 
	writeData(dta, "dta_DNAaccessibility.txt", row.names=T)

	# GM12878LCL 
	dta <- getENCODEmetadata(file_format="bed", biosample="GM12878", assay="ChIP-seq")
	writeData(dta, "GM12878", row.names=T)

}

getENCODEmetadata <- function(...) {
	#dta <- queryEncode(fixed=F, set_accession="ENCSR623GZY", file_format="fastq")$experiment
	#dta <- queryEncode(fixed=F, biosample="C2C12", file_format="bed")$experiment
	#head(dta)
	print(list(...))
	dta <- queryEncode(fixed=F, ...)$experiment
	dta <- dta[,c("platform", "controls", "href", "accession", "biosample_name", "biological_replicate_number", "technical_replicate_number", "paired_end", "paired_with", "target","biosample_type", "organism", "date_released", "lab", "assay")]
	dta
}


