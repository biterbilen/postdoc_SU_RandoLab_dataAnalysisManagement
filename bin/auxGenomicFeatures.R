source("~/Projects/biter_biter_shared/bin/bigData.R")

library(Biostrings)
library(GenomicFeatures)
library(parallel)
library(rtracklayer)

findUpdatedMirName <- 
	function(x, y, colName="name") {
		x.updated <- 
			y[grepl("MIMAT", y$V1) & 
				grepl(paste(x,";",sep=""), y$V2, fixed=T) &
				(grepl("*", x, fixed=T) | (!grepl(paste(x,"*;",sep=""), y$V2, fixed=T))), colName] 
		if (length(x.updated) == 0) 
			x.updated <- as.character(x)
		x.updated[1]
	}

caller.updateMirIds <- 
	function(mirExp.file="originalTaqManRodentmiRNA_QSC_ASC.txt.gz",
					 ofle="TaqManRodentmiRNA_QSC_ASC.txt.gz", 
					 mirAliases.file="../miRBase/aliases.txt.gz",
					 mc.cores=1) 
	{

		dta <- readBigTable(mirExp.file, header=T)

		# read aliases file to update the mirNames in expression file 
		aliases <- readBigTable(mirAliases.file)
		# last entry in V2 should be the most updated name of miR
		aliases$name <- gsub("(([^;]+);)+", "\\2", aliases$V2, perl=T)

		dta[,1] <- unlist(lapply(dta[,1], 
														 function(x) findUpdatedMirName(x,aliases,"name")))
		print(dim(dta))

		if (! is.null(ofle))
			writeData(dta, ofle, row.names=F)

		dta
	}

test.get.mirSeedRevComp <-
	function() {

		source("~/Projects/biter_biter_shared/bin/auxGenomicFeatures.R")
		mirSeq.file <- "/home/biter/PI_HOME/Data/casco/miRBase/mature.fa.gz"
		org <- "mmu"
		seed.start=2;
		seed.end=8
		dta <- get.mirSeedRevComp(mirSeq.file, org, seed.start, seed.end, mirExp.name="Name")

	}


get.mirSeedRevComp <- 
	function(mirSeq.file, org="mmu", seed.start=2, seed.end=8, mirExp.name="Name", select=NULL, frmt="data.frame") {
		# read sequence file
		seqs <- readRNAStringSet(mirSeq.file)
		names(seqs) <- gsub("(.+) MIMAT.+", "\\1", names(seqs))
		if (!is.null(org))
			seqs <- subset(seqs, grepl(org, names(seqs)))
		if (!is.null(select)) 
			seqs <- subset(seqs, grepl(select, names(seqs)))
		# extract seeds
		seedRevComps <- DNAStringSet(reverseComplement(Biostrings::subseq(seqs, start=seed.start, end=seed.end)))

		if (frmt == "data.frame") {
			seedRevComps <- as.data.frame(seedRevComps) 
			seedRevComps[,mirExp.name] <- row.names(seedRevComps)
			names(seedRevComps)[1] <- "seedRevComp"
		}

		seedRevComps
	}

test.caller.getkmerFreq <-
	function() {
		source("~/Projects/biter_biter_shared/bin/auxGenomicFeatures.R")

		mirSeq.file <- "/home/biter/PI_HOME/Data/casco/miRBase/mature.fa.gz"

		seedRevComps <- get.mirSeedRevComp(mirSeq.file, org="mmu", frmt="DNAStringSet")

		pdict0 <- PDict(seedRevComps)
		genseq <- getBSgenome("mm10")
		genseq[[3]]
		chrY <- genseq$chrY
		#mi0 <- matchPDict(pdict0, chrY)     # Search...
		names
		mic0 <- vcountPDict(pdict0, DNAStringSet(list(genseq[["chrM"]], genseq[["chrY"]])))     # Search...
		rownames(mic0) <- names(pdict0)
		colnames(mic0) <- seqnames(genseq)

	}

caller.getkmerFreq <- 
	function(fafile, fafilepat, type="DNA", typepat="DNA", start=1, end=200, ...) {
		f <- ifelse (type == "DNA", get("readDNAStringSet"), get("readRNAStringSet")) 
		seqs <- f(fafile)

		fpat <- ifelse (typepat == "DNA", get("readDNAStringSet"), get("readRNAStringSet")) 
		seqspat <- fpat(fafilepat)

		mirseqs <- subset(mirseqs, grepl(org, names(seqs)))
		#alphabetFrequency(seqs)
		seq.tab <- data.frame(t(consensusMatrix(subseq(seqs, start=start, end=end), as.prob=T)[1:4,]))
		seq.tab$position <- 1:nrow(seq.tab) - (end-start+1)/2
		seq.tab
	}

test.caller.nucProbPerPosition <-
	function() {
		source("~/Projects/biter_biter_shared/bin/auxGenomicFeatures.R")
		fafile <- "./tmpself.seqs.gz"
		dta <- caller.nucProbPerPosition(fafile, start=1, end=200)
	}

# case insensitive
caller.nucProbPerPosition <- 
	function(fafile, type="DNA", start=1, end=200, ...) {
		if (type == "DNA") { 
			f <- get("readDNAStringSet")
		} else {
			f <- get("readRNAStringSet")
		}
		seqs <- f(fafile)
		#alphabetFrequency(seqs)
		seq.tab <- data.frame(t(consensusMatrix(subseq(seqs, start=start, end=end), as.prob=T)[1:4,]))
		seq.tab$position <- 1:nrow(seq.tab) - (end-start+1)/2
		seq.tab
	}


getBSgenome <- 
	function(gnm) {
		if (gnm == "mm10") {
			library("BSgenome.Mmusculus.UCSC.mm10")
			genSeq <- Mmusculus
		}
		if (gnm == "mm9") {
			library("BSgenome.Mmusculus.UCSC.mm9")
			genSeq <- Mmusculus
		}
		if (gnm == "hg19") {
			library("BSgenome.Hsapiens.UCSC.hg19")
			genSeq <- Hsapiens
		}
		if (gnm == "hg18") {
			library("BSgenome.Hsapiens.UCSC.hg18")
			genSeq <- Hsapiens
		}
		genSeq
	}


test.get.revCompSeedsOfExpressedMirFamilies <-
	function() {
		source("auxGenomicFeatures.R")
		# w/ expression filter, TODO
		mirSeq.file <- "~/PI_HOME/Data/casco/miRBase/mature.fa.gz"
		mirExp.file <- "~/Projects/biter_jbrett_epigenetics/data/Cheung_2014_Nature/TaqManRodentmiRNA_QSC_ASC.txt.gz"
		ofle="expressedwExpressionFilter.txt.gz"

		# wo/ expression filter, OD
		mirSeq.file <- "~/PI_HOME/Data/casco/miRBase/mature.fa.gz"
		mirExp.file <- NA
		ofle <- "expressedwoExpressionFilter.txt.gz"
		dta <- get.revCompSeedsOfExpressedMirFamilies(mirSeq.file, mirExp.file=mirExp.file, ofle=ofle)

		mirSeq.file <- "~/PI_HOME/Data/casco/miRBase/mature.fa.gz"
		mirExp.file <- "~/Projects/biter_wosczyna_miRNAsInMuscleCellFate/data/Wosczyna_2013/selected.txt.gz"
		ofle <- "expressedwoFilter.txt.gz"
		dta <- get.revCompSeedsOfExpressedMirFamilies(mirSeq.file, mirExp.file=mirExp.file, ofle=ofle)

 }	 

get.revCompSeedsOfExpressedMirFamilies <- 
	function(mirSeq.file="mature.fa.gz", mirExp.file=NA,
					 gnm="mm10", mirExp.diffLogRatio=5, exp.filter=F,
					 #mirExp.background="MammU6", mirExp.cols=c("AvgctQSC","Avgct3.5D"), mirExp.name="Geneid", mirExp.max.CT=30, # TODO
					 mirExp.background="MammU6", mirExp.cols=c("AvgctQSC","Avgct3.5D"), mirExp.name="Geneid", 
					 seed.start=2, seed.end=8, mc.cores=1, ofle=NULL)
	{

		if (gnm == "mm10" || gnm == "mm9") {
			org <- "mmu"
		} else if (gnm == "hg19" || gnm == "hg18") { 
			org <- "hsa"
		} else {
			message("miR seed reverse complement extraction is not supported for ", gnm)
			quit();
		}

		seedRevComps <- get.mirSeedRevComp(mirSeq.file, org, seed.start, seed.end, mirExp.name)

		if (!is.na(mirExp.file)) {
			expr <- readBigTable(mirExp.file, header=T)
			print(dim(expr))
			# read expression file and filter for expressed and normalize expression
			if (exp.filter) {
				# filter for expressed and normalize expression
				expr <- expr[,c(mirExp.name,mirExp.cols)]
				expr <- subset(expr, apply(expr, 1, function(x,y=mirExp.max.CT) any(x < y))) 
				# normalize expression wrt mirExp.background
				normFactor <- unlist(expr[expr[mirExp.name]==mirExp.background, -1])
				normFactor <- matrix(rep(normFactor,  dim(expr)[1]), nrow=dim(expr)[1], byrow=T)
				relExprInLog2 <- normFactor - expr[,-1]
				names(relExprInLog2)[1:2] <- paste(names(relExprInLog2)[1:2], ".wrt", mirExp.background, sep="")
				relExprInLog2$CTdiff <- relExprInLog2[,2] - relExprInLog2[,1] # relative to the first assuming first one is the control
				relExprInLog2[mirExp.name] <- expr[mirExp.name]
				print(dim(relExprInLog2))
				# get org specific and differentially expressed by putting a cut-off 
				# TODO add more advanced options to select differentially expressed; also working for >2 groups
				diffExpressedOrg <- subset(relExprInLog2, 
																	 abs(CTdiff) > mirExp.diffLogRatio & grepl(org, relExprInLog2[,mirExp.name]))
				print(dim(diffExpressedOrg))

			} else {
				expr %<>% select_(mirExp.name)
			}
			diffExpMirSeedRevComps <- merge(expr, seedRevComps)
		} else {
			diffExpMirSeedRevComps <- seedRevComps
		}
		rownames(diffExpMirSeedRevComps) <- diffExpMirSeedRevComps$Geneid

		# patch for aggregate
		diffExpMirSeedRevComps[,mirExp.name] <- as.character(diffExpMirSeedRevComps[,mirExp.name])
		diffExpMirSeedRevComps[,"seedRevComp"] <- as.character(diffExpMirSeedRevComps[,"seedRevComp"])
		print(c("diffExpMirSeedRevComps", dim(diffExpMirSeedRevComps)[1], dim(diffExpMirSeedRevComps)[2]))

		# finally get the families aggregating the seeds
		dta <- aggregate(. ~ seedRevComp, diffExpMirSeedRevComps, paste, collapse="|")
		print(c("aggdiffExpMirSeedRevComps", dim(dta)[1], dim(dta)[2]))
		if (! is.null(ofle))
			writeData(dta, ofle, row.names=F)

		dta
	}

test.caller.expressedTxRegionsOrSeqs <- 
	function() {
		#caller.expressedTxRegionsOrSeqs(expressedTx.file="~/Projects/biter_jbrett_epigenetics/results/2014-11-28/back.bs20_edgeRnorm_Q18_Q4/STAR_Q4_Q18_DE.txt", regions=c("promoters", "threeUTRsByTranscript", "intronsByTranscript", "fiveUTRsByTranscript"))

		#expressedTx.file="~/Projects/biter_jbrett_epigenetics/results/2014-11-28/back.bs20_edgeRnorm_Q18_Q4/STAR_Q4_Q18_DE.txt"
		#expressedTx.file="./head_DEseq_WT_KO_DE.txt"
		#regions=c("promoters")
		#expressedTx.colInd=1
		#excludeTx.file=NULL 
		#regions.only=F
		#ofleTag=NULL
		#gnm="mm10"
		#remove.ambiguous=T
		#mc.cores=1
		# Returns regions list of DNAStringSet
		source("~/Projects/biter_biter_shared/bin/auxGenomicFeatures.R")

		# extract promoters using refGene table of UCSC
		expressedTx.file <- "~/PI_HOME/Data/casco/UCSC_tracks/mm10/mm10_refGene_canonical.gtf.gz"
		regions <- "promoters"
		regions.only <- T
		cord <- caller.expressedTxRegionsOrSeqs(expressedTx.file, regions=regions, 
																						regions.only=T,ofleTag="promoters")[[1]]
		# extract promoters using refGene table of UCSC
		expressedTx.file <- "~/PI_HOME/Data/casco/UCSC_tracks/mm10/mm10_refGene_4test.gtf.gz"
		regions <- "transcripts"
		regions.only <- F
		cord <- caller.expressedTxRegionsOrSeqs(expressedTx.file, regions=regions, 
																						regions.only=T,ofleTag="tmp")[[1]]

		# extract sequences
		regions.only <- F
		regions <- "self"
		#expressedTx.file <- cord # GRanges object
		expressedTx.file <- "~/PI_HOME/Data/casco/UCSC_tracks/mm10/mm10_refGene_canonical_promoters.bed.gz"

		expressedTx.colInd <- 1 #mcols of GRanges
		dss <- caller.expressedTxRegionsOrSeqs(expressedTx.file, expressedTx.colInd=expressedTx.colInd, regions=regions, regions.only=regions.only, ofleTag="promoterseqs")

		dss <- caller.expressedTxRegionsOrSeqs(expressedTx.file, expressedTx.colInd=expressedTx.colInd, regions=regions, regions.only=regions.only, ofleTag="promoterseqs")

		fle <- "~/PI_HOME/Data/casco/ENCODE/hg19/Synapse_syn6131484/labels/MYC.train.labels.tsv.gz" 
		caller.expressedTxRegionsOrSeqs(fle, expressedTx.frmt="bed", nrows=10, header=T, verbose=T)[[1]] 
		readBigTable(fle, col.names=c("chr","start","end","name","score","strand"), frmt="GRanges", nrows=10, verbose=T, header=T)
		readBigTable(fle, col.names=c("chr","start","end","name","score","strand"), nrows=10, verbose=T, header=T)
		readBigTable(fle, col.names=c("chr","start","end"), nrows=10, verbose=T, header=T)

		# extract polyA sequences using refGene table of UCSC
		source("~/Projects/biter_biter_shared/bin/auxGenomicFeatures.R")
		#expressedTx.file <- "~/PI_HOME/Data/casco/UCSC_tracks/mm10/canonical/polyA.bed.gz"
		expressedTx.file <- "~/PI_HOME/Data/casco/UCSC_tracks/mm10/canonical/a.bed.gz"
		regions <- "self"
		seqs <- caller.expressedTxRegionsOrSeqs(expressedTx.file, regions=regions, 
																						regions.only=F,ofleTag="tmp")[[1]]

		# Download different regions
		source("~/Projects/biter_biter_shared/bin/auxGenomicFeatures.R")
		fle <- "~/PI_HOME/Data/casco/SRA/Wu_2017_BBA/bmsc_d00_15_histone_rx2_ctcf_smc1a_all_time_mm10_dense.bed.gz"
		regionStr <- "self"
		ofleTag <- "dummy"
		gnm <- "mm10"
		col.names=c("chr","start","end","name","score","strand", "start2", "end2", "color")
		col.names=c("chr","start","end","name","score","strand")

		k <- readBigTable(fle, header=T, skip=1)
		atxr <- readBigTable(fle, col.names=col.names, frmt="GRanges", skip=1, filterStr="name==7")
		head(atxr)
		caller.expressedTxRegionsOrSeqs(expressedTx.file=fle, regions=regionStr, gnm=gnm, ofleTag=ofleTag, debug=T, col.names=col.names, skip=1)

		source('~/Projects/biter_biter_shared/bin/auxGenomicFeatures.R')
		fle <- "/share/PI/casco/Data/casco/SRA/Wu_2017_BBA/bmsc_d00_15_histone_rx2_ctcf_smc1a_all_time_mm10_dense.bed.gz"
		regions <- 'self'
		ofleTag <- regions
		genome <- 'mm10'
		linesSkipped <- 1
		col.names=c("chr","start","end","class"); filterStr="class==7"
		verbose <- T
		atxr <- readBigTable(fle, col.names=col.names, filterStr=filterStr, skip=linesSkipped,frmt="GRanges", posStrand=T, verbose=verbose)
	  regionSeqs <- caller.expressedTxRegionsOrSeqs(expressedTx.file=fle, regions=regions, gnm=genome, skip=linesSkipped, ofleTag=ofleTag, debug=T, filterStr="class==7", col.names=col.names, verbose=T, posStrand=T)


		# CDS
		fle <- "/share/PI/casco/Data/casco/UCSC_tracks/mm10/mm10_refGene_canonical.v2.gtf.gz"
		regions <- "cdsBy"
	  regionSeqs <- caller.expressedTxRegionsOrSeqs(expressedTx.file=fle, regions=regions, gnm=genome, skip=linesSkipped, ofleTag=ofleTag, debug=T, filterStr="class==7", col.names=col.names, verbose=T, posStrand=T)
	}

caller.expressedTxRegionsOrSeqs <- 
	function(expressedTx.file, expressedTx.colInd=1, expressedTx.frmt=NULL, excludeTx.file=NULL, excludeTx.colInd=1, 
					 col.names=c("chr","start","end","name","score","strand"),
					 regions=c("self", "transcripts", "promoters", "threeUTRsByTranscript", "intronsByTranscript", "fiveUTRsByTranscript", "cdsBy"),
					 regions.only=F, ofleTag=NA, promoter.up = 2000, promoter.down = 1000,
					 gnm="mm10", remove.ambiguous=T, mc.cores=1, debug=F, verbose=F, ...)
	{
		#print(expressedTx.file)
		# Expressed transcriptIds
		transcriptIds <- NULL
		atxr <- NULL
		if ( is.character(expressedTx.file)) {
			if (any(grep(".gtf", expressedTx.file, fixed=T)) || ((!is.null(expressedTx.frmt)) && expressedTx.frmt == "gtf")) {
				atxr <- import.gff(expressedTx.file)
				#transcriptIds <- gsub(".*transcript_id \\\"(.*?)\\\";.*", "\\1", atxr$group)
				transcriptIds <- unique(as.character(atxr$transcript_id))
			} else if (any(grep(".bed", expressedTx.file, fixed=T)) || ((!is.null(expressedTx.frmt)) && expressedTx.frmt == "bed")) {
				#atxr <- readBigTable(expressedTx.file, col.names=c("chr","start","end","name","score","strand"), frmt="GRanges")
				atxr <- readBigTable(expressedTx.file, col.names=col.names, frmt="GRanges", verbose=verbose, ...)
				#atxr <- import.bed(expressedTx.file) # this doesn't work for chr4_L456216_random
				transcriptIds <- as.character(atxr$name)
			} else {
				#transcriptIds <- as.character(readBigTable(expressedTx.file, header=T)[, expressedTx.colInd])
				transcriptIds <- as.character(readBigTable(expressedTx.file, header=T, verbose=verbose, ...)[, expressedTx.colInd])
			}
		} else {
			atxr <- expressedTx.file
			transcriptIds <- mcols(atxr)[,expressedTx.colInd]
		}

		if (debug) {
			transcriptIds <- transcriptIds[1:10]
			atxr <- atxr[1:10,]
		}

		print(length(transcriptIds))

		# Remove excludeTx 
		if ( !is.null(excludeTx.file)) {
			ex.transcriptIds <- as.character(readBigTable(excludeTx.file, header=T)[, excludeTx.colInd])
			transcriptIds <- setdiff(transcriptIds, ex.transcriptIds)
		}
		print(length(transcriptIds))

		# Get expressed transcript structures if self is not the only structure
		tx <- NULL 
		if (!any(grep("self", regions))) 
			tx <- makeTxDbFromUCSC(gnm, tablename="refGene", transcript_ids=transcriptIds)

		# extract sequences
		dta <- NULL 
		genSeq <- getBSgenome(gnm)
		r <- "self"
		r <- "promoters";
		r <- "threeUTRsByTranscript"
		r <- "transcripts"
		r <- "cdsBy"
		for (r in regions) {
			print(r)
			f <- NULL
			if (r != "self") 
				f <- get(r)
			if (r == "self") {
				txr <- atxr
				txr <- mclapply(1:length(txr), function(i) { txi <- txr[i]; } , mc.cores=mc.cores) 
				names(txr) <- transcriptIds
				txr <- GRangesList(txr)
			} else if (r == "promoters") {
				# Make GRangesList from Granges object for extractTranscriptSeqs function call
				txr <- f(tx, upstream=promoter.up, downstream=promoter.down)
				nms <- values(txr)[,2]
				txr <- mclapply(1:length(txr), function(i) { txi <- txr[i]; } , mc.cores=mc.cores) 
				names(txr) <- nms
				txr <- GRangesList(txr)
			} else if (r == "transcripts") {
				txr <- f(tx, columns="TXNAME") #RefSeq Accession
				names(mcols(txr)) <- "name"
				names(txr) <- as.character(txr$name)
			} else {
				txr <- f(tx, use.names=T)
			}

			# some refseq IDs map to multiple positions in the genome
			# remove ambiguous entries 
			if (remove.ambiguous) 
				txr <- txr[!duplicated(names(txr))]

			print(length(txr))
			#txr <- txr[1:3,]
			if (regions.only) {
				if (! is.na(ofleTag))
					writeData(unlist(txr), paste(ofleTag, r, ".bedplus.gz", sep=""), frmt="bedplus", col.names=F, row.names=F)
				dta <- c(dta, unlist(txr))
				rm(txr)
			} else { 
				txrs <- mclapply(1:length(txr), 
												 function(i, x=genSeq) {
														tx1 <- txr[i]
														# TODO makes sure this works
														#save(list=ls(), file=paste(ofleTag, ".RData",sep=""), envir=environment())
														extractTranscriptSeqs(x, tx1)
												}, mc.cores=mc.cores)
				dss <- DNAStringSet(mclapply(txrs, function(x)  x[[1]] ))
				names(dss) <- mclapply(txrs, function(x)  names(x) )
				if (! is.na(ofleTag))
					writeXStringSet(dss, paste(ofleTag, r, ".fa.gz", sep=""), compress=T) 
				dta <- c(dta, dss)
				rm(txrs)
				rm(dss)
			}
		}
		names(dta) <- regions

		dta
	}


test.getDataFromUCSC <- 
	function() {
		source("~/Projects/biter_biter_shared/bin/auxGenomicFeatures.R") 
		range <- GRanges(seqnames=Rle(c('chr1', 'chr2', 'chr3'), c(1, 1, 1)),IRanges(1:3, width=5000000))
		TSSregion <- getDataFromUCSC(track="refGene",table="refGene",process=list(name="TSS",slopu=2000,slopd=1000),frmt="GRanges")
		TSS <- getDataFromUCSC(track="refGene",table="refGene",process=list(name="TSS",slopu=0,slopd=0),frmt="GRanges")
		head(TSS)
		TRXasBed <- getDataFromUCSC(track="refGene",table="refGene",process=list(name="TRX",slopu=0,slopd=0), ofle="refGene.gz")
		head(TRXasBed)
		TSSasBed <- getDataFromUCSC(track="refGene",table="refGene",process=list(name="TSS",slopu=0,slopd=0))
		head(TSSasBed)
		cons <- getDataFromUCSC(track="cons60way", table="multiz60way", range=range)
		head(cons)
		refgene <- getDataFromUCSC(track="refGene",table="refGene",ofle="refGene.gz")
		head(refgene)
		cytoBand <- getDataFromUCSC(track="cytoBand",table="cytoBand",ofle="cytoBand.gz")
		head(cytoBand)
	}

# TODO make a local database if downloaded for the first time
getDataFromUCSC <- 
	function(track="refGene",table="refGene",range=NULL,gnm="mm10",frmt=NULL,process=NULL,printAllTracks=T,printAllTables=T, ofle=NULL) {

		if (is.null(range))
			range <- gnm

		session <- browserSession()
		genome(session) <- gnm
		if (printAllTracks) 
			print(trackNames(session))

		query <- ucscTableQuery(session, track, range=range)
		if (printAllTables)
			print(tableNames(query))

		# Extract table
		tableName(query) <- table;
		dta <- getTable(query)
		print(dim(dta))
		print(head(dta,3))

		# process
		if (!is.null(process)) {
			if (process$name == "TSS") { # hard coded for table refGene fields 
				startEnd <- t(apply(dta[,c("strand","txStart","txEnd")], 1, 
														function(x) { 
															if(x[1]=="+") c(as.numeric(x[2]) - process$slopu,
																							as.numeric(x[2]) + process$slopd + 1) 
														else c(as.numeric(x[3]) - process$slopd - 1,
																	 as.numeric(x[3]) + process$slopu) } ))
				dta <- data.frame(dta$chrom, startEnd, dta[,"name"], 1000, dta[,"strand"])
				names(dta) <- c("chr", "start", "end", "name", "score", "strand")
			} else if (process$name == "TES") { # hard coded for table refGene fields 
				startEnd <- t(apply(dta[,c("strand","txStart","txEnd")], 1, 
														function(x) { 
															if(x[1]=="-") c(as.numeric(x[2]) - process$slopu,
																							as.numeric(x[2]) + process$slopd + 1)
														else c(as.numeric(x[3]) - process$slopd - 1,
																	 as.numeric(x[3]) + process$slopu) } ))
				dta <- data.frame(dta$chrom, startEnd, dta[,"name"], 1000, dta[,"strand"])
				names(dta) <- c("chr", "start", "end", "name", "score", "strand")
			} else if (process$name == "TRX") { # hard coded for table refGene fields 
				startEnd <- t(apply(dta[,c("strand","txStart","txEnd")], 1, 
														function(x) { 
															if(x[1]=="+") c(as.numeric(x[2]) - process$slopu,
																							as.numeric(x[3]) + process$slopd) 
														else c(as.numeric(x[2]) - process$slopd,
																	 as.numeric(x[3]) + process$slopu) } ))
				dta <- data.frame(dta$chrom, startEnd, dta[,"name"], 1000, dta[,"strand"])
				names(dta) <- c("chr", "start", "end", "name", "score", "strand")
			} else if (process$name == "5HEXONS") { # hard coded for table refGene fields 
				# TODO continue
				startEnd <- t(apply(dta[4:6,c("strand","exonStarts","exonEnds")], 1, 
														function(x) { 
															exonStarts <- as.numeric(strsplit(x[2],",")[[1]]);
															exonEnds <-   as.numeric(strsplit(x[3],",")[[1]]);
															lens <- exonEnds - exonStarts
															csumlens <- cumsum(lens)
															halflen <- floor(sum(lens) / 2)
															index <- grep(T, csumlens > halflen)[1]
															nExonStarts <- NULL
															nExonEnds   <- NULL
															breach <- halflen
															if (index > 1) {
																nExonStarts <- exonStarts[1:(index-1)]
																nExonEnds   <- exonEnds[1:(index-1)]
																breach <- breach - cumsum(lens[1:(index-1)]) 
															}
															nExonStarts <- paste(c(nExonStarts,exonStarts[index]), collapse=",")
															nExonEnds <- paste(c(nExonEnds, exonStarts[index]+breach), collapse=",")
															print(c(lens,0,halflen,"100",nExonStarts,"000",nExonEnds))
															print(nExonEnds)
														} ))
				startEnd
				dta[4:6,1:10]
				dta <- data.frame(dta$chrom, startEnd, dta[,"name"], 1000, dta[,"strand"])
				names(dta) <- c("chr", "start", "end", "name", "score", "strand")
			} else {
				message("ERROR:", process$name, " process is not implemented")
			}
		}

		if (!is.null(frmt))
			dta <- convertDataFormat(dta, frmt=frmt)

		if (!is.null(ofle))
			writeData(dta, ofle, row.names=F, col.names=F)

		dta
	}	


aux.rtracklayer.dummy <- 
	function() {
		gnm <- "mm10"
		library(rtracklayer)
		session <- browserSession()
		genome(session) <- gnm
		#	tm <- trackNames(session)
		#	names(tm)
		tm["RefSeq Genes"]
		#tm["miRNA"]
		tm["Conservation"]
		tm

		query <- ucscTableQuery(session, "refGene")
		tableNames(query)
		tableName(query) <- "refGene"
		table.dta <- getTable(query)
		TSS <- t(apply(table.dta, 1, function(x, slop=1000) { if(x[4]=="+") c(TSSstart=as.numeric(x[5])-slop,TSSend=as.numeric(x[5])+slop) else c(TSSstart=as.numeric(x[6])-slop,TSSend=as.numeric(x[6])+slop) } ))

		TSSwcoord <- cbind(table.dta$chrom, TSS, table.dta[,c("name","exonCount","strand")])
		TSSwcoord
		#names(TSSwcoord) <- c("seqname", "start", "end", "name", "score", "strand")
		names(TSSwcoord) <- c("chr", "start", "end", "name", "score", "strand")
		TSS.gr <- convertDataFormat(TSSwcoord, frmt="GRanges")
		head(TSS.gr)

		head(table.dta,2)

		length(names(tm))
		query <- ucscTableQuery(session, "Conservation")

		Conservation.gr <- track(query)
		query <- ucscTableQuery(session, "Conservation", GRangesForUCSCGenome(gnm, "chr12",IRanges(57795963, 57815592)))
		tableName(query) <- "multiz30waySummary"
		?tableName
		names(tm[91:123])
		getTable(query)

		tableNames(query)
		ucscSchema(query)                                                     
		getTable(query)
		query <- ucscTableQuery(session, "cpgIslandExt") 
		CpGislands.gr <- track(query)
		CpGislands.gr

		data(cpneTrack)
		plot(start(cpneTrack), score(cpneTrack))

		annoF <- "./mmu.gff3.gz"
		annoGR <- import(annoF, format="gff3", asRangedData=FALSE)
		head(annoGR)

		canonical <- import.gff("~/Projects/biter_jbrett_epigenetics/data/UCSC_tracks/mm10/mm10_refGene_canonical.gtf.gz")
		head(canonical)
	}

fun.rphast.dummy <- 
	function() {

		library(rphast)
		mf <- read.feat("./mmu.gff3")
		mf.pri <- mf[mf$feature=="miRNA_primary_transcript",]
		head(mf)
		f <- read.feat("a.gtf", pointer.only=F)
		f <- add.UTRs.feat(f)
		f <- add.introns.feat(f)
		fu <- unique(f)
		dim(f)
		dim(fu)

		f3 <- f[f$feature=="3'UTR", ]
		#table(f$feature)
		composition.feat(mf.pri, f)
		f.ov <- overlap.feat(mf.pri, f)
		plot.gene(f[1:9,])
		plot.feat(f[1:9,])

		exampleArchive <- system.file("extdata", "examples.zip", package="rphast")
		filename <- "rev.mod"
		unzip(exampleArchive, filename)
		tm <- read.tm(filename)
		plot(tm)
		plot(tm, show.eq.freq=FALSE)
		plot(tm, max.cex=20, eq.freq.max.cex=1,
				 col=matrix(1:16, nrow=4),
				 eq.freq.col=c("red", "green"),
				 filled=TRUE, add=TRUE)

		m <- msa(seqs=c("ACGTAT", "AGGTAA", "AGGTAG"),
						 names=c("human", "mouse", "rat"))
		print(m)
		print(m, format="FASTA")
		print(m, format="PHYLIP", pretty.print=TRUE)

	}

