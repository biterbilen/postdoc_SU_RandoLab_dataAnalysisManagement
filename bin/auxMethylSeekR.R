library(MethylSeekR)
source("/home/biter/Projects/biter_biter_shared/bin/auxggplot2.R")
source("/home/biter/Projects/biter_biter_shared/bin/bigData.R")
source("/home/biter/Projects/biter_biter_shared/bin/auxGenomicFeatures.R")

reduceCpGCounts2PlusStrand <- function(dta, genSeq, num.cores=1) {
	# GRanges object to extract sequences
	all.gr <- GRanges(seqnames=unlist(dta[1]), ranges=IRanges(unlist(dta[2]+1), unlist(dta[3])))
	values(all.gr) <- data.frame(T=unlist(dta[4]+dta[5]), M=unlist(dta[4]))

	# only segments chromosomes with at least nCpG.smoothing covered CpGs
	nCGsPerChr <- table(as.character(seqnames(all.gr)))
	chrs <- names(nCGsPerChr)

	res <- mclapply(chrs, function(chr){
		sel <- which(as.character(seqnames(all.gr))==chr)
		
		gr  <- all.gr[sel,]

		sq <- getSeq(genSeq, gr)
		if((sq=="C"||sq=="G") == F)
			message("Error in input: positions corresponds to non-C bases")

		# Prep for reduction by extending +(C) downstream and -(G) upstream
		end(gr)[sq=="C"] <- end(gr)[sq=="C"]+1
		start(gr)[sq=="G"] <- start(gr)[sq=="G"]-1

		# Reduce
		combined.gr <- reduce(gr,min.gapwidth=0)
		ov <- findOverlaps(gr, combined.gr)
		T <- tapply(values(gr[queryHits(ov)])[, 1], subjectHits(ov), sum)
		M <- tapply(values(gr[queryHits(ov)])[, 2], subjectHits(ov), sum)
		values(combined.gr) <- data.frame(T=T, M=M)

		#Revert to single nucleotide 
		end(combined.gr) <- end(combined.gr) - 1

		# remove 0 coverage positions
		combined.gr <- combined.gr[values(combined.gr)$T>0]
		combined.gr
	}, mc.cores=num.cores)
	all.combined.gr <- do.call(c, unname(res))
}

extractData4CombinedSegments <- 
	function(rdsfiles, all.segs, cols=c("median.meth"), pdfsplomFilename=NULL,
					 pdfdendrogramFilename=NULL, TableFilename=NULL, verbose=F) {

		library(GenomicRanges)
		dta <- NULL
		nms <- NULL
		for (f in rdsfiles) {
			nm <- gsub(".rds", "", basename(f))
			nms <- c(nms,
							 unlist(lapply(cols, function(cl) { paste(nm, cl, sep=".") })))
			x <- readRDS(f)
#			writeData(x, paste(f,".bedplus.gz", sep=""), frmt="bedplus", row.names=F, verbose=T)
			ov <- findOverlaps(all.segs,x)
			cur.col <- data.frame(matrix(rep(NA,length(all.segs)*length(cols)),nrow=length(all.segs),ncol=length(cols)))
			cur.col[queryHits(ov),] <- unlist(values(x[subjectHits(ov)])[cols])
			if (is.null(dta)) {
				dta <- data.frame(cur.col)
			} else {
				dta <- data.frame(dta, cur.col)
			}
		}
		names(dta) <- nms
		values(all.segs) <- dta

		if (!is.null(pdfsplomFilename)) 
			get.splom(dta, outfileName=pdfsplomFilename);

		if (!is.null(pdfdendrogramFilename)) 
			clusterHier(dta, outfileName=pdfdendrogramFilename,verbose=T);

		if (!is.null(TableFilename)) 
			writeData(all.segs, TableFilename, frmt="bedplus", row.names=F, verbose=verbose)
		all.segs		
	}

getCombinedSegments <-
	function(files, regionType=NULL, fileFormat="rds", combinationType="union", verbose=F) {
		segments.gr <- NULL
		cur.segments.gr <- NULL
		for (f in files) {
			if (fileFormat == "rds")
				cur.segments.gr <- readRDS(f)
			#TODO add other file formats
			if (is.null(segments.gr)) {
				segments.gr <- cur.segments.gr;
			} else {
				if (! is.null(regionType)) 
					cur.segments.gr <- cur.segments.gr[values(cur.segments.gr)$type==regionType]
				if (combinationType=="union")
					segments.gr <- union(segments.gr, cur.segments.gr)
				if (combinationType=="intersection")
					segments.gr <- intersect(segments.gr, cur.segments.gr)
			}
		}
		if (verbose == T) 
			message(paste("Combined segments N=", length(segments.gr), sep=""))
		segments.gr
	}

# TODO set DNAmethyl=T for Jamie's data analysis
getCommonSegments <-
	function(seg1, seg2, GRangesFilename=NULL, TableFilename=NULL, pdfFilename=NULL, 
					 xlabtag=NULL, ylabtag=NULL, DNAmethyl=F) {

		# get common segments; should map one to one
		ov <- findOverlaps(seg1, seg2)
		message("Common segments N=", length(ov))
		segs <- subsetByOverlaps(seg1[queryHits(ov)], seg1)

		if (DNAmethyl) {
			values(segs) <- 
				DataFrame(nCG.segmentation=values(segs)$nCG.segmentation, 
									median.meth.1=values(segs)$median.meth,
									median.meth.2=values(seg2[subjectHits(ov)])$median.meth)
		} 

		# save GRanges object
		if(!is.null(GRangesFilename))
			saveRDS(segs, GRangesFilename)

		# save a tab-delimited table
		if(!is.null(TableFilename)){
			#D=data.frame(chr=as.character(seqnames(segs)), start=start(segs), end=end(segs), type=values(segs)$type, nCG.segmentation=values(segs)$nCG.segmentation, nCG.seq=values(segs)$nCG, mean.meth=100*values(segs)$pmeth, median.meth=100*values(segs)$median.meth)
			D <- data.frame(chr=as.character(seqnames(segs)), start=start(segs), end=end(segs), 
											nCG.segmentation=values(segs)$nCG.segmentation, 
											median.meth.1=100*values(segs)$median.meth.1,
											median.meth.2=100*values(segs)$median.meth.2)
			write.table(D, file=TableFilename, quote=FALSE, sep="\t", row.names=FALSE)
		}
		# scatterplot
		if (!is.null(pdfFilename)) {
			pdfFilename <- paste(tag,"_",type,".pdf", sep="")
			jet.colors <- colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan","#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))
			pdf(pdfFilename, height=5, width=5)
			smoothScatter(100*values(seg.gr1[queryHits(ov),])$median.meth, 
										100*values(seg.gr2[subjectHits(ov),])$median.meth, 
										colramp=jet.colors, 
										xlab=paste("median methylation (%) -",xlabtag, sep=""), 
										ylab=paste("median methylation (%) -",ylabtag, sep="")) 
			dev.off()
		}

		#return
		segs
	}

plotScatter <- function(segments.gr, x="nCG", y="median.meth", xlab="number of CpGs in segment (log2)",
												ylab="median methylation (%)", main="", pdfFilename = NULL) {
	jet.colors <- colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan","#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))
	if(!is.null(pdfFilename)){
		pdf(pdfFilename, height=5, width=5)
	}
	#smoothScatter(log2(values(segments.gr)[x]), 100*values(segments.gr)[y], colramp=jet.colors, xlab=xlab, ylab=ylab);
	smoothScatter(log2(values(segments.gr)[,x]), values(segments.gr)[,y], colramp=jet.colors, xlab=xlab, ylab=ylab, main=main);
	#abline(v=log2(nCG.classification), lty=5)
	if(!is.null(pdfFilename))
		dev.off()
}


getMethylation <- 
	function(methFile, segments.gr, gnm, main="",
					 nCpG.smoothing = 3, minCover = 5, nCpG.cutoff = 3, PMDs = NULL, pdfFilename = NULL, num.cores = 1) {

		genSeq <- getBSgenome(gnm)
		seqLengths <- seqlengths(genSeq)

		# ready methylome	
		m <- readMethylome(FileName=methFile, seqLengths=seqLengths)
		if (! is.null(PMDs))
			m <- setdiff(m, PMDs)

		# select CpGs with minimal coverage
	  nCpGall <- length(m)
		m <- m[values(m)[, 1]>=minCover]
		nCpG <- length(m)
		CpGwminCover <- nCpG/nCpGall*100
		message(paste("Above minCover(",minCover,") CpGs N=",nCpG, "(",CpGwminCover,")",sep=""))

		ov <- findOverlaps(m, segments.gr)

		# update segments.gr in case some segments are notmapped by any reads
		segments.gr <- segments.gr[unique(subjectHits(ov)),]

		if (length(segments.gr) == 0)
			return(NULL)

		ov <- findOverlaps(m, segments.gr)

		# Fields
		#nCG <- vcountPattern("CG", getSeq(genSeq, segments.gr))
		nCG <- vcountPattern("CG", getSeq(genSeq, resize(segments.gr, width(segments.gr), fix="start"), as.character=FALSE))#

		T <- tapply(values(m[queryHits(ov)])[, 1], subjectHits(ov), sum)
		M <- tapply(values(m[queryHits(ov)])[, 2], subjectHits(ov), sum)
		nCG.segmentation <- tapply(values(m[queryHits(ov)])[, 1], subjectHits(ov), length)

		median.meth <- tapply(as.vector(runmean(Rle(values(m[queryHits(ov)])[, 2]/values(m[queryHits(ov)])[, 1]), nCpG.smoothing, endrule="constant")), subjectHits(ov), median)

		median.meth <- pmax(0, median.meth) # conversion to Rle causes rounding errors that sometimes give very small negative numbers

		if(!isTRUE(all.equal(as.numeric(names(T)), 1:length(segments.gr)))) {
			message("error in calculating methylation levels")
		}

		# Patch for C(G) map at the end of the segment
		indx <- nCG.segmentation > nCG
		nCG[indx] <- nCG.segmentation[indx]

		values(segments.gr) <- DataFrame(nCG.segmentation, nCG, T, M, pmeth=M/T, median.meth=median.meth, coveredCG.ratio=nCG.segmentation/nCG)

		# remove segments represented with less than nCpG.cutoff sufficiently covered CpGs
		if (! is.null(nCpG.cutoff))
			segments.gr <- segments.gr[values(segments.gr)$nCG.segmentation >= nCpG.cutoff,]		

		plotScatter(segments.gr, pdfFilename=pdfFilename, main=main)

		segments.gr
	}

# TODO rewrite using combineRegions function
getFMRs <- 
	function(UMRsLMRs, PMDs=NULL, num.cores=1, fine.tune=F) {
		library(GenomicRanges)
		library(parallel)

		# segments in between UMRsLMRs
		FMRs <- gaps(UMRsLMRs)
		FMRs <- FMRs[strand(FMRs)=="*"]

		# restrict FMRs on chromosomes of UMRsLMRs
		chrs <- as.character(runValue(seqnames(UMRsLMRs)))
		res <- mclapply(chrs, 
										function(chr) { FMRs[seqnames(FMRs)==chr] }, 
										mc.cores=num.cores)
		FMRs <- do.call(c, res)

		# remove segments overlapping with PMDs
		if (!is.null(PMDs)) {
			notPMDs <- gaps(PMDs)
			notPMDs <- notPMDs[strand(notPMDs)=="*"]
			FMRs <- subsetByOverlaps(FMRs, notPMDs)
		}

		# remove too short segments; had coded here as 5
		if (fine.tune == T) {
			FMRs <- FMRs[width(FMRs)>5,]
			print("TODO fine tune based on CG count etc")	
		}
		FMRs
	}

# TODO write
finetune <- function( meth.cutoff = 0.5, nCpG.smoothing = 3) {
	print("TODO: organize")
	return;

	# only segments chromosomes with at least nCpG.smoothing covered CpGs
	nCGsPerChr <- table(as.character(seqnames(m)))
	chrs <- names(nCGsPerChr)[nCGsPerChr>=nCpG.smoothing]

	res <- mclapply(chrs, function(chr) {
									#select chr
									sel <- which(as.character(seqnames(m))==chr)
									# methylation rate of CpGs
									rle <- Rle(values(m)[sel, 2]/values(m)[sel, 1])
									# smoothed methylation rate
									mean.meth <- runmean(Rle(values(m)[sel, 2]/values(m)[sel, 1]), k=nCpG.smoothing, endrule="constant")
									# most likely hypomethylated CpGs 
									indx <- mean.meth < meth.cutoff

									# remove segments that are too short
									# first hypomethylated regions
									runValue(indx)[runLength(indx)<nCpG.cutoff & runValue(indx)==TRUE]=FALSE
									# then methylated regions
									runValue(indx)[runLength(indx)<nCpG.cutoff & runValue(indx)==FALSE]=TRUE

									# ids 
									tmp.ids <- rep(1:length(runLength(indx)), runLength(indx))
									tmp <- split(1:length(sel), tmp.ids)
									tmp <- tmp[runValue(indx)==TRUE]
									if(length(tmp)>0) {
										coords <- cbind(sapply(tmp, min), sapply(tmp, max))
										starts <- round((start(m)[sel[pmax(1, coords[, 1]-1)]]+start(m)[sel[coords[, 1]]])/2) # take middle between two CpGs (if at the beginning or end of chromosome, don't extend
										ends <- round((start(m)[sel[coords[, 2]]]+start(m)[sel[pmin(length(sel), coords[, 2]+1)]])/2)

										hmr.gr <- GRanges(seqnames=unique(seqnames(m[sel])), strand="*", ranges=IRanges(starts, ends), seqlengths=seqLengths)
									}
									else{
										hmr.gr <- GRanges(, seqlengths=seqLengths)
									}
									hmr.gr
										}, mc.cores=num.cores)

	segments.gr <- do.call(c, unname(res))
}

combineColumnsFromCommonIntervals <- 
	function(all.gr, subset=list(), select=NULL, RDSfiles=NULL, nms=NULL, outfle.tag=NULL, debug=F, minSegmentCount=100) {

	if (debug) 
		save(list=ls(), file=paste("combineColumnsFromCommonIntervals", outfle.tag, ".RData", sep=""), envir=environment())	

	if (! is.null(RDSfiles))  
		all.gr <- lapply(RDSfiles, function(x) readRDS(x))

	if (! is.null(nms))
		names(all.gr) <- nms

	# one subset criteria is supported
	x <- all.gr[[1]]
	all.gr <- lapply(all.gr, 
									 function(x) { 
										 tmp <- subset(x, get(names(subset)) > subset[[1]], select=select)
										 print(paste("subset", length(tmp)))
										 tmp
									 })
	print(length(all.gr))

	# common in between all
	#intersect.gr <- Reduce(subsetByOverlaps, all.gr)
	# exact match
	intersect.gr <- Reduce(function(x,y) subsetByOverlaps(x, y, type="equal"), all.gr)
	print(paste("common", length(intersect.gr)))

	if (debug) 
		save(list=ls(), file=paste("combineColumnsFromCommonIntervals2", outfle.tag, ".RData", sep=""), envir=environment())	

#load("combineColumnsFromCommonIntervals2pmeth_OtherClass.RData"); library(GenomicRanges)
	# extract data columns
	dta <- as.data.frame(lapply(all.gr,
															function(x) {
																#ov <- findOverlaps(intersect.gr, x) 
																# exact match
																ov <- findOverlaps(intersect.gr, x, type="equal")
																print(ov)
																tmp <- values(x[subjectHits(ov)])
																print(paste("common", dim(tmp)[1]))
																tmp
															}))
	names(dta) <- paste(rep(select, length(all.gr)), rep(names(all.gr), each=length(select)), sep=".")

	# Graphics output
	if (length(select) == 1 && (!is.null(outfle.tag)) && dim(dta)[1] > minSegmentCount) {
		print("Graphics output")
		get.exploratoryAnalysisPlots(dta, ofle.tag=outfle.tag, title=outfle.tag)
	}

	# Text output
	print("Text output")
	nms <- c("seqnames", "start", "end", names(dta))
	ov <- findOverlaps(intersect.gr, all.gr[[1]], type="equal")
	dta <- cbind(as.character(seqnames(intersect.gr[queryHits(ov)])),
								start(intersect.gr[queryHits(ov)]) - 1,
								end(intersect.gr[queryHits(ov)]),
								dta)
	names(dta) <- nms

	print(head(dta))
	writeData(dta, paste(outfle.tag, ".tab.gz", sep=""), row.names=F)	

	return(dta)
}

caller.estimate.combine.callDMR <- 
	function(regionfile, msrFiles, msrFiles.names=NULL, regionfile.wheader=F, gnm="mm10", 
					 nCpG.cutoff=NULL, min.CpGcover=5, 
					 coveredCG.ratio=0, cols2save="pmeth,T:M", #cols2save=list(c("pmeth"), c("T","M")),
					 DMRmodel=NULL, samplename.pat="Q4,A4,R5", bs=1, nc=1, 
					 outfle.tag=NULL, save.image=F, debug=F) {
	msrFiles <- strsplit(msrFiles, ",")[[1]]
	samplename.pat <- strsplit(samplename.pat, ",")[[1]]

	print(date())

	# read region file
	seg <- readBigTable(regionfile, header=regionfile.wheader) 	
	names(seg) <- c("chr", "start", "end", "name", "score", "strand"); #patch for conversion
	seg.gr <- convertDataFormat(seg, frmt="GRanges")

	print(seg.gr)
	if (length(seg.gr) < 100) {
		message("Not processing region with ", length(seg.gr), " segments")
		return()
	}

	if (is.null(msrFiles.names)) {
		msrFiles.names <- 
		 unlist(lapply(msrFiles, function(x) sub(".*group_(.*).msr.gz", "\\1", x,perl=TRUE)))
	} else {
		msrFiles.names <- strsplit(msrFiles.names, ",")[[1]]
	}

	print("msrFiles.names")	
	print(msrFiles.names)	
	i <- 1
	all.gr <- lapply(1:length(msrFiles), 
									 function(i) {
										 x <- msrFiles[i]
										 print(i)
										 tag <- paste(msrFiles.names[i], strsplit(basename(regionfile), ".", fixed=T)[[1]][1], sep="_")
										 getMethylation(x, seg.gr, gnm, num.cores=nc,
																	 nCpG.cutoff=nCpG.cutoff, minCover=min.CpGcover, 
																	 main=tag)

									 })
	names(all.gr) <- msrFiles.names

	if (debug) 
		save(list=ls(), file=paste("caller.estimate.combine.callDMR", outfle.tag, ".RData", sep=""), envir=environment())	

	cols2save <- strsplit(cols2save, ",")[[1]]
	x <- cols2save[1]
	dta <- lapply(cols2save,
								function(x) {
									tag <- paste(paste(x, collapse=""), strsplit(basename(regionfile), ".", fixed=T)[[1]][1], sep="_")
									print("tag")
									print(tag)
									combineColumnsFromCommonIntervals(all.gr, subset=list(coveredCG.ratio=coveredCG.ratio), select=strsplit(x,":")[[1]], outfle.tag=tag, debug=debug)
								})
	print("ahahha")
	print(length(dta))


	print("Starting DMR model with this data")
	print(head(dta[[1]]))
	print(head(dta[[2]]))
	
	saveRDS(dta, paste("dta", outfle.tag, ".RDS", sep=""))

	if (! is.null(DMRmodel)) { 
		fle <- paste(outfle.tag, ".tab.gz", sep="")
		
		if (DMRmodel == "pairedBaySeq") {
			pairCD <- baySeqMDdensity(fle=fle, bs=bs, nc=nc, pats=samplename.pat) 
		} else {
			message("DMRmodel", DMRmodel, "is not defined");
			quit()
		}
	}

	if (save.image)
		save.image(paste(outfle.tag,".RData",sep=""))
	#	saveRDS(dta, paste(outfle.tag,".RDS", sep=""))

	sessionInfo()
	print(date())
	dta
}
