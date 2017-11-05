source("~/Projects/biter_biter_shared/bin/auxTFBSTools.R"); 
source("~/Projects/biter_biter_shared/bin/bigData.R")
#source("~/Projects/biter_biter_shared/bin/auxggplot2.R") plot functions should not be here
library(parallel)
library(dplyr)
library(tidyr)
library(mclust)

test.all <-
	function() {
		# TODO run
		source("auxML.R")
		test.GMMcdens()
		test.mGMM()
		test.GMM()
		test.caller.LM()
		test.caller.getOverlap()
		test.caller.getMEnrichment.wFET4SegmentCountPerMatrix()
		test.getDNAshapeVectorInRegion()
		test.getMotifCountsInRegion()
		test.caller.getSeqFeatures()
		test.caller.getExpressed()
		test.caller.DESeq()
		test.caller.DE.edgeR()
		test.GSEA()
		test.caller.normalize()
		test.caller.normalize()
		test.caller.profileNormalization.DE()
	}

test.GMMcdens <- 
	function() {
		source("auxML.R")
		fle.train <- "polyA5p_nuc.txt"
		fle.test <- "pooled5p_dist10_summit_nuc.txt"
		ofleTag <- "pooled5p_dist10_summit_nuc"
		param.Mclust=list(G=c(1))
		param.preprocess=list(col.sel.pat="num_[ACT].*nuc5p")
		param.preprocess=list(col.sel.pat="num_[AT]")
		GMMcdens("polyA5p_nuc.txt", "pooled5p_dist10_summit_nuc.txt", "pooled5p_dist10_summit_nuc")

		# PAS for huge test file
		fle.train <- "distalpolyA5p_nuc.txt"
		fle.test <- "pooled5p_dist10_summit_nuc.txt"
		train.param.preprocess <- list(col.sel.pat='num_[AT]')
		test.param.preprocess <- c(subset.cond='X5_usercol.pooled5p_dist10_summit_nuc5p>5',train.param.preprocess)
		test.param.preprocess <- c(subset.cond='X5_usercol.pooled5p_dist10_summit_nuc5p>5')
		ofleTag <- NULL;
		GMMcdens(fleTrain, fleTest, ofleTag, param.preprocess=param.preprocess)
}			      

# GMM for density estimation
GMMcdens <- 
	function(fle.train, fle.test=NULL, ofleTag=NULL, 
					 train.param.preprocess=list(col.sel.pat="num_[ACT]"), 
					 test.param.preprocess=list(subset.cond='X5_usercol.pooled5p_dist10_summit_nuc5p>5', col.sel.pat="num_[ACT]"), 
					 param.Mclust=list(G=1), scoreLimit=0, verbose=T, ...) {
		print(ofleTag)

		# PAS_20161118   
		dta.train <- preprocess(dta=readBigTable(fle.train, header=T), param=train.param.preprocess)
		dta.test.all <- readBigTable(fle.test, header=T)

		# TODO fix the memory problem in function calling! is there a less IP-dependent version?
		#source("auxML.R")
		#dta.test.all <- preprocess(dta=dta.test.all, param=test.param.preprocess, verbose=T) %>% head
		dta.test.all = preprocess(dta=dta.test.all, param=test.param.preprocess, verbose=T) %>% head

		# TODO
		dta.test.1 <- dta.test.all %>% filter_("X5_usercol.pooled5p_dist10_summit_nuc5p>0")
		dta.test.2 <- dta.test.all %>% filter(X5_usercol.pooled5p_dist10_summit_nuc5p>0)
		dta.test.3 <- dta.test.all %>% filter_("X5_usercol.pooled5p_dist10_summit_nuc5p>0")
		dta.test.4 <- dta.test.all %>% filter(X5_usercol.pooled5p_dist10_summit_nuc5p>5)
		dta.test %>% dplyr::select(X5_usercol.pooled5p_dist10_summit_nuc5p) %>% head(10)
		# plot 
		source("~/Projects/biter_biter_shared/bin/auxggplot2.R")
		if (! is.null(ofleTag)) {
			plot.new(); pdf(paste(ofleTag,"_GMMcdens.pdf", sep=""),h=10,w=10) 
		}
		my_corrheat_seg(dta.train, size=rel(4), title="train")
		get.splom(dta.train)
		# model
		mdl <- do.call(Mclust, c(list(data=dta.train), param.Mclust))
 		print(summary(mdl, parameters=T))
		plot(mdl, "density")
		#plot(mdl, "classification")
		plot(mdl, "BIC")
		#plot(mdl, "uncertainty") #does not work

		# tweak from https://github.com/fisproject/anomaly-detection/blob/master/ch03/mclust.R
		mdl.d <- cdens(modelName=mdl$modelName, dta.test, parameters=mdl$parameters)

		print("mdl.d")
		print(summary(mdl.d))
		ps <- 1
		dta.test.all$PASdensity <- log2(as.matrix(mdl.d+ps))
		min4shift <- min(dta.test.all$PASdensity)
		min4shift <- ifelse (min4shift < 0, -min4shift, 0)
		print(summary(dta.test.all$PASdensity))
		# pi is mixture prob but we use G=1 for our case so this is not necessary
		# anomary <- -log(as.matrix(mdl.d) %*% as.matrix(pi)) 
		dta.test <- cbind(dta.test, PASdensity.PASdensity=dta.test.all$PASdensity)
		get.splom(dta.test)
		my_corrheat_seg(dta.test, size=rel(4), title="test", noheatmap=T)

		cor(dta.test)

		dta.test.all$class <- GMM(dta.test.all[,"PASdensity", drop=F], title="PASdensity", selected=F, verbose=T)

		if (! is.null(ofleTag)) {
			dev.off()
			writeData(dta.test.all, paste0(ofleTag,"_GMMcdens.txt"), row.names=F)
		}
	}

test.mGMM <-
	function() {
		source("~/Projects/biter_biter_shared/bin/auxML.R")
		# multivariate
		data <- iris; columns=1:4
		mGMM(data, columns=columns) %>% head

		data <- iris; columns=1:4
		mGMM(data[1:4]) %>% head

		fle <- "~/Projects/biter_wosczyna_miRNAsInMuscleCellFate/results/2017-04-05/FC_Wosczyna_2015_PullDown/featureCounts.txt.gz"
		fle <- "~/Projects/biter_wosczyna_miRNAsInMuscleCellFate/results/2017-04-05/FC_Wosczyna_2015_RNAseq/featureCounts.txt.gz"
		ps <- 1 
		data <- readBigTable(fle, header=T, row.names='Geneid') %>% sample_frac(1)
		columns <- grep(".nsrt.bam", names(data))
		#--------------
		#1.
		#param.Mclust <- list(G=2, modelNames="V")
		#data %<>% mutate_at(columns, funs(log((.+ps)/Length)))
		#2. 
		param.Mclust <- list(G=2)
		data %<>% mutate_at(columns, funs(log(.+ps)))
			
		k <- mGMM(data, columns=columns, param.Mclust=param.Mclust)
		k <- mGMM(data[,columns], param.Mclust=param.Mclust)

	}

mGMM <- 
	function(data, columns=1:ncol(data), param.Mclust=list(G=2, modelNames="V")) {
		GMMa <- function(dta) { print(colnames(dta)); GMM(dta, param.Mclust=param.Mclust, selected=F, factored=F, plot=F) }
		print(data.frame(mGMM=colnames(data)[columns]))
		datac <- 
			data %>% 
			mutate_at(columns, funs(GMMa)) %>%
			mutate(NotClass1=as.factor(ceiling(rowMeans(.[columns]))))
		datac[,columns] <- data[,columns]
		datac
	}

test.GMM <- 
	function() {
		source("~/Projects/biter_biter_shared/bin/auxML.R")
		dta1 <- data.frame(value=c(15, rnorm(100, 6, 1), rnorm(9, 15, 1)))
		hist(dta1$value)
		head(dta1)
		ofleTag <- "GMM.test"
		param.Mclust <- list(G=2, modelNames="V")
		selected <- GMM(dta1,ofleTag=ofleTag, param.Mclust=param.Mclust, title="TITLE")
		# The first component has a larger variance
		dta <- data.frame(value=c(15, rnorm(90, 6, 1), rnorm(1000, 9, 1)))
		hist(dta$value)
		param.Mclust <- list(modelNames="V")
		param.Mclust <- list(G=2, modelNames="V")
		selected <- GMM(dta, param.Mclust=param.Mclust, title="TITLE")

		dta[,"value"] <- dta.test.all[,"PASdensity", drop=F] + 1
		head(dta)
		hist(dta$value)
		param.Mclust <- list(G=2, modelNames="V")
		param.Mclust <- list(G=1:2)
		dta$group <- GMM(dta, param.Mclust=param.Mclust, title="TITLE", selected=F)

	}

# Gaussian Mixture Model for classification
# input is a data.frame with row.names; otherwise the function has an intrinsic unknown order
# TODO generalize for G>2
# integrate to mGMM and remove this
GMM <- 
	function(dta, ofleTag=NULL, param.Mclust=list(G=2, modelNames="V"), selected=T, factored=T,
					 verbose=T, plot=F, ...) {

		dta <- unname(dta)

		mdl <- do.call(Mclust, c(list(data=dta), param.Mclust))

		#1. This discrimination filtered little
		#border <- mdl$parameters$mean[2] - (1 * sqrt(mdl$parameters$variance$sigmasq[2]))
		#mdl$classification <- ifelse (mdl$data > border , 2, 1)
		#2. Correct for very highly expressed genes
		ind <- mdl$data > mdl$parameters$mean[2] & mdl$classification == 1
		mdl$classification[ind] <- 2
		#3. Correct for very low expressed genes
		ind <- mdl$data < mdl$parameters$mean[1] & mdl$classification == 2
		mdl$classification[ind] <- 1

		#plot(mdl, "density")
		#plot(mdl, "classification")
		#plot(mdl, "BIC")
		#plot(mdl, "uncertainty")

		# select bigger GMM component 
		if (verbose) 
			print(summary(mdl, parameters=T))

		# plot whole data and class densities
		if (plot) {
			dta <- rbind(data.frame(value=unname(mdl$data), group=0),
									 data.frame(value=unname(mdl$data), group=mdl$classification))

			dta$group <- as.factor(dta$group)

			my_dens(dta, 
							mapping=aes(x=value, color=group), 
							columns=1, 
							long=F, 
							size=rel(2), ...) %>% print
		}

		if (selected) 
			rownames(dta)[which(mdl$classification==2)]
		else
			if (factored) 
				as.factor(mdl$classification)
			else
				mdl$classification
	}

# TODO continue
test.caller.LM <- 
	function() {
		source("~/Projects/biter_biter_shared/bin/auxML.R")
		regions.file="promoters.bed.gz"
		pmeth.file="canonical_Q4Q10Q16Q22_promoters_pmeth_promoters.tab.gz"
		TFBS.file="./promotersJASPAR2014.txt.gz"
		TFBS.file="JASPAR2014_promoters.txt.gz"
		g1 <- "pmeth.Q22"
		g2 <- "pmeth.Q4"
		ps <- 0.0001
		TFBS.name <- "Homer"
		ofleTag <- paste(TFBS.name, "_", g1, "_o_", g2, sep="")
		caller.LM(regions.file=regions.file, pmeth.file=pmeth.file, TFBS.file=TFBS.file, g1=g1, g2=g2, ofleTag=ofleTag, ps=ps, TFBS.name=TFBS.name)
}

# TODO generalize
caller.LM <- 
	function(regions.file="promoters.bed.gz", 
					 pmeth.file="./canonical_A4Q4R5_pmeth_promoters.tab.gz", 
					 TFBS.file="JASPAR2014_promoters.txt.gz",
					 g1="pmeth.A4", g2="pmeth.Q4", ofleTag=paste(g1,g2,sep="_o_"), ps=0.0001, TFBS.name=".MA") {

		# TODO call from bigData.R joinFiles function
		regions <- readBigTable(regions.file, col.names=c("seqnames","start","end","name","score","strand"))
		pmeth <- readBigTable(pmeth.file, header=T)
		TFBS <- readBigTable(TFBS.file, header=T); names(TFBS)[1] <- "name"
		print(head(TFBS,2))
		adta <- merge(merge(regions[,1:4], pmeth), TFBS) 
		print("names adta")
		print(names(adta))

		g1s <- strsplit(g1, "|", fixed=T)[[1]]
		g2s <- strsplit(g2, "|", fixed=T)[[1]]
		ig1 <- grep(g1, names(adta), perl=T, value=T)
		ig2 <- grep(g2, names(adta), perl=T, value=T)
		iTFBS <- grep(TFBS.name, names(adta), fixed=T,value=T)
		print(iTFBS)
		
		# dta structure prep
		library(ggplot2)
		library(reshape2)
		plot.new(); pdf(paste(ofleTag, ".pdf",sep=""),h=8,w=8)
		summary(adta[c(ig1,ig2)])
		if (length(g1s) == 1 && length(g2s) == 1) {
			dta.pmeth <- data.frame(pmeth.g1=rowMeans(adta[ig1]), pmeth.g2=rowMeans(adta[ig2]))
		} else if (length(g1s) == 2 && length(g2s) == 2) {
			dta.pmeth <- data.frame(pmeth.g1=rowMeans(adta[grep(g1s[1], ig1, value=T)]) / rowMeans(adta[grep(g1s[2], ig1, value=T)]),
															pmeth.g2=rowMeans(adta[grep(g2s[1], ig2, value=T)]) / rowMeans(adta[grep(g2s[2], ig2, value=T)]))
		} else {
			stop("Operation not defined for", g1, " or ", g2)
		}
		dta.pmeth$pmeth.g1og2 <- (dta.pmeth[,1]+ps)/(dta.pmeth[,2]+ps)
		dta.pmeth$log2.pmeth.g1og2 <- log2(dta.pmeth$pmeth.g1og2)
		dta.pmeth$pmeth.g1dg2 <- dta.pmeth[,1] - dta.pmeth[,2]
		summary(dta.pmeth)
		ggplot(melt(dta.pmeth[,1:2]), aes(x=value, color=variable)) + geom_density()
		ggplot(melt(dta.pmeth[,3]), aes(x=value)) + geom_density() + labs(x="pmeth.g1og2")
		ggplot(melt(dta.pmeth[,4]), aes(x=value)) + geom_density() + labs(x="log2.pmeth.g1og2")
		ggplot(melt(dta.pmeth[,5]), aes(x=value)) + geom_density() + labs(x="pmeth.g1dg2")
		ggplot(dta.pmeth, aes(x=pmeth.g1dg2, y=log2.pmeth.g1og2)) + geom_bin2d() + labs(x="pmeth.g1dg2", y="log2.pmeth.g1og2")

		dev.off()
		print(head(dta.pmeth))
		print(head(adta[,iTFBS]))

		# Gaussian LM
		i <- 5
		for (i in 4:5) { 
			dta <- data.frame(y=dta.pmeth[i], adta[,iTFBS])
			print(names(dta))
			frml <- paste(names(dta)[1], paste0(names(dta)[-1], collapse="+"), sep="~")
			#frml <- paste(names(dta)[1], paste0(names(dta)[-1], collapse="*"), sep="~")
			print(frml)
			mdl <- lm(as.formula(frml), data=dta)
			k <- options("max.print")[[1]]
			options(max.print=100000)
			print(summary(mdl))
			options(max.print=k)
		}


		# TODO check the model parameter plots
		#plot(mdl)
	}

test.caller.getOverlap <- 
	function() {
	de.fle <- "~/Projects/biter_lingliu_DNAdamage/results/2015-10-12/DE_none_Liu_2015_RNAseq/norm_DEG_scnone_RSRUVs_DE_edgeR.txt.gz"
	peak.fle <- "~/Projects/biter_lingliu_DNAdamage/results/2015-05-19/JAMM_peaks_2011_yQSC_H3K9me3_vs_yQSC_DNA/outdir/peaks/filtered.peaks.narrowPeak"
	FDRcut <- 0.02
	FDRcut <- 0.1
	scorecut <- 100
	source("~/Projects/biter_biter_shared/bin/auxML.R")
	#range <- GRanges(seqnames=Rle(c('chr1'), c(1)),IRanges(1:1, width=500000))

	# gene coordinates
	TRX.gr <- getDataFromUCSC(track="refGene",table="refGene",process=list(name="TRX",slopu=0,slopd=0), frmt="GRanges")

	#can.gr <- import.gff("~/Projects/biter_jbrett_epigenetics/data/UCSC_tracks/mm10/mm10_refGene_canonical.gtf.gz")
	#ids <- gsub(".* \"(.+)\";", '\\1', unique(mcols(can.gr)$group), perl=T)
	exp <- readBigTable(de.fle, header=T)
	exp.gr <- NULL
	exp.gr$sig <- TRX.gr[which(elementMetadata(TRX.gr)$name %in% exp[exp$FDR < FDRcut,"Id"])]
	exp.gr$insig <- TRX.gr[which(elementMetadata(TRX.gr)$name %in% exp[exp$FDR >= FDRcut,"Id"])]

	# cytoband coordinates
	cytoBand <- getDataFromUCSC(track="cytoBand",table="cytoBand")
	gieStains <- as.character(unique(cytoBand$gieStain))

	# chip coordinates
	peak <- readBigTable(peak.fle)
	peak.gr <- NULL
	peak.gr$H3K9me3.foreground <- convertDataFormat(peak[peak$V5>scorecut,1:6], frmt="GRanges", nms=c("chr", "start", "end", "name", "score", "strand"))
	peak.gr$H3K9me3.background <- convertDataFormat(peak[peak$V5<=scorecut,1:6], frmt="GRanges", nms=c("chr", "start", "end", "name", "score", "strand"))

	# subtelomeric and interstitial
	subint <- markSubIntTelomeres(gnm="mm10", len=150000)
	types <- as.character(unique(subint$type))

	# gieStains
	m <- data.frame(matrix(0, ncol=length(types), nrow=length(names(exp.gr))), row.names=c("sig", "insig")) 
	colnames(m) <- types
	gs <- types[1]
	for (gs in types) {
		cb <- subint[subint$type == gs,1:3];
		cb.gr <- convertDataFormat(cb, frmt="GRanges", nms=c("chr", "start", "end"))
		s <- names(exp.gr)[2]
		for (s in names(exp.gr)) {
			# One to one mapping (???)
			m[s,gs] <- length(getCommonSegments(exp.gr[s][[1]], cb.gr))
		}
	}
	print(m)
	print(chisq.test(m, simulate.p.value = T))
	print(m[1,] / m[2,])

	# si 
	m <- data.frame(matrix(0, ncol=length(gieStains), nrow=length(names(exp.gr))), row.names=c("sig", "insig")) 
	colnames(m) <- gieStains
	gs <- gieStains[1]
	for (gs in gieStains) {
		cb <- cytoBand[cytoBand$gieStain == gs,1:3];
		cb.gr <- convertDataFormat(cb, frmt="GRanges", nms=c("chr", "start", "end"))
		print(cb)
		s <- names(exp.gr)[1]
		for (s in names(exp.gr)) {
			# One to one mapping (???)
			m[s,gs] <- length(getCommonSegments(exp.gr[s][[1]], cb.gr))
		}
	}
	print(m)
	print(chisq.test(m, simulate.p.value = T))
	print(m[1,] / m[2,])


	# chip
	m <- data.frame(matrix(0, ncol=length(names(peak.gr)), nrow=length(names(exp.gr))), row.names=c("sig", "insig")) 
	colnames(m) <- names(peak.gr)
	gs <- names(peak.gr)[1]
	for (gs in names(peak.gr)) {
		s <- names(exp.gr)[1]
		for (s in names(exp.gr)) {
			# One to one mapping (???)
			m[s,gs] <- length(getCommonSegments(exp.gr[s][[1]], peak.gr[gs][[1]]))
		}
	}
	print(m)
	print(chisq.test(m, simulate.p.value = T))
	print(m[1,] / m[2,])

}

caller.getOverlap <- 
	function() {


}


markSubIntTelomeres <- 
	function(gnm="mm10", mc.cores=1, len=15000) {
	chrs <- paste("chr", c(1:19, "X", "Y"), sep="")
	lens <- seqlengths(getBSgenome(gnm))
	si <- mclapply(chrs, function(x) { rbind(c(x, 0, len, "subtelomeric"), 
																	 c(x, len, lens[x][[1]]-len, "interstitial"),
																	 c(x, lens[x][[1]]-len, lens[x][[1]], "subtelomeric")) }, mc.cores=mc.cores )
	
	res <- data.frame(do.call("rbind", si))
	names(res) <- c("chr", "start", "end", "type")
	res$type <- as.character(res$type)
	res$chr <- as.character(res$chr)
	res$start <- as.integer(as.character(res$start))
	res$end <- as.integer(as.character(res$end))
	res

}

test.caller.getMEnrichment.wFET4SegmentCountPerMatrix <- 
	function() {
		source("~/Projects/biter_biter_shared/bin/auxTFBSTools.R")
		bedfile <- "./DMRwDSSwPresetCoverageCutoffs_Q4Q22_FDR05/Q4Q22_FDR05.txt.gz"
		bedwHeader=F
		slop=12
		caller.getMEnrichment.wFET4SegmentCountPerMatrix(bedfile, bedwHeader=bedwHeader, slop=slop)
	}

# TODO convert list to PWMatrixList object for PFM shuffle method, OK
# TODO read JASPAR Weirauch_2014_Cell PSM and PWM and use
# TODO sequence shuffle method
# TFBS enrichment other than Fisher's Exact test
caller.getMEnrichment.wFET4SegmentCountPerMatrix <- 
	function(bedfile, matrixSource="JASPAR2014", background="shuffledMatrix",
					 Nbackground=1, gnm="mm10", bedwHeader=T,
					 min.score="95%",
					 bed2=F, slop=NULL, topN=NULL, select=NULL, mc.cores=1, 
					 decreasing=T, ofletag=NULL) { 
		print(date())

		# Read coordinate file 
		dta <- readBigTable(bedfile, header=bedwHeader)
		if (! is.null(select))
			dta <- dta[order(dta[select], decreasing=decreasing),]; #name
		dim(dta)

		if (is.null(topN) || dim(dta)[1] < topN) {
			Nsegments <- dim(dta)[1]
		} else { 
			Nsegments <- topN
		}

		# Get bed3
		if (bed2) # zero-based bed files
			dta[,3] = dta[,2] + 1;
		dta <- dta[,1:3]

		# Convert to Granges for sequence search
		names(dta) <- c("chr", "start", "end") # patch for Granges conversion
		dta.gr <- convertDataFormat(dta, frmt="GRanges")

		print(length(dta.gr))

		#seqLogo(toICM(PFMatrixList[[1]]))
		#seqLogo(toICM(PFMatrixList.b[[1]]))

		# Count segments with Matrix
		# Foreground
		PFMatrixList <- getMatrixList(src=matrixSource, gnm=gnm, matrixtype="PFM", all_versions=T)

		nrows <- 1:Nsegments
		seqs <- getSeqsFromGRanges(dta.gr[nrows,], gnm, wnames=T, slop=slop)
		sites.counts <- getSiteCountswMatrixAndDNAStringSet(PFMatrixList, seqs, min.score=min.score, add.seqname=F, mc.cores=mc.cores)

		Mpresence <- as.data.frame(colSums(sites.counts), optional=T)
		names(Mpresence) <- paste("segmentCounts","topNsegments","originalMatrixList", sep=".")
		cbind(Mpresence, Mpresence)

		rm(sites.counts)

		n <- NULL
		if (background == "shuffledMatrix") {
			print(background)
			for (i in 1:Nbackground) {
				print(paste(background, i))
				PFMatrixList.b <- permuteMatrix(PFMatrixList)
				sites.counts <- getSiteCountswMatrixAndDNAStringSet(PFMatrixList.b, seqs, min.score=min.score, add.seqname=F, mc.cores=mc.cores)
				Mpresence <- cbind(Mpresence, as.data.frame(colSums(sites.counts)))
				names(Mpresence)[i+1] <- paste("segmentCounts","topNsegments","shuffledMatrixList", i, sep=".")

				rm(PFMatrixList.b)
				rm(sites.counts)
			}
			n <- paste("segmentCounts", "topNsegments","shuffledMatrixList", "median", sep=".") 
		}
		rm(seqs)

		if (background == "nextTopNsegments") {
			print(background)
			for (i in 1:Nbackground) {
				nrows <- (i*Nsegments+1):(i*Nsegments+Nsegments)
				print(head(nrows))
				if (dim(dta)[1] < nrows[length(nrows)]) {
					Nbackground <- i - 1
					break
				}
				seqs.b <- getSeqsFromGRanges(dta.gr[nrows,], gnm, wnames=T, slop=slop)


				sites.counts <- getSiteCountswMatrixAndDNAStringSet(PFMatrixList, seqs.b, min.score=min.score, add.seqname=F, mc.cores=mc.cores)

				Mpresence <- cbind(Mpresence, as.data.frame(colSums(sites.counts)))
				names(Mpresence)[i+1] <- paste("segmentCounts", "nextTopNsegments", i, "originalMatrixList", sep=".") 

				rm(seqs.b)
				rm(sites.counts)
			}
			n <- paste("segmentCounts", "nextTopNsegments", "median", "originalMatrixList", sep=".") 
		} else if (background == "mutuallyExclusiveRandomNsegments") {
			print(background)
			for (i in 1:Nbackground) {
				nrows <- sample(1:dim(dta)[1], Nsegments, 
												prob=c(0*(1:Nsegments), 0*(1:(dim(dta)[1]-Nsegments))+1))
				print(head(nrows))
				seqs.b <- getSeqsFromGRanges(dta.gr[nrows,], gnm, wnames=T, slop=slop)

				sites.counts <- getSiteCountswMatrixAndDNAStringSet(PFMatrixList, seqs.b, min.score=min.score, add.seqname=F, mc.cores=mc.cores)
				Mpresence <- cbind(Mpresence, as.data.frame(colSums(sites.counts)))
				names(Mpresence)[i+1] <- paste("segmentCounts", "mutuallyExclusiveRandomNsegments", i, "originalMatrixList", sep=".") 

				rm(seqs.b)
				rm(sites.counts)
			}
			n <- paste("segmentCounts", "mutuallyExclusiveRandomNsegments", "median", "originalMatrixList", sep=".") 
		}
		rm(PFMatrixList)

		# Add median counts as the last segment count for the background model
		Mpresence <- cbind(Mpresence, apply(Mpresence[,-1], 1, median))
		names(Mpresence)[Nbackground+2] <- n
		Mpresence

		# Counts Fisher's exact test statistics for contingency table of segment counts w/wo Matrix in foreground and background
		ncols <- c(1,Nbackground+2)
		cs <- cbind(            Mpresence[,ncols], 
								Nsegments - Mpresence[,ncols])
		Mpresence <- cbind(Mpresence, getEnrichment(cs, alternative="g")) 

		# order based on p-value
		pvaluei <- grep("p.value", names(Mpresence))
		Mpresence <- Mpresence[order(Mpresence[, pvaluei]),]
		Mpresence

		# write
		msStr <- strsplit(basename(matrixSource), ".", fixed=T)[[1]][1]
		if (! is.null(ofletag)) {
			ofle <- paste(ofletag, "_matrixSource", msStr, "_minscore", min.score, "_background", background, "_top",topN,"Segments", "_slop", slop, sep="")
			ofle <- paste(ofle,".txt", sep="")
			writeData(Mpresence, ofle, row.names=T, col.names=NA)

			# plot segment counts with motif
			if (background == "nextTopNsegments") {
				pvali <- grep("p.value", names(Mpresence))
				Mpresence.sel <- Mpresence[Mpresence[pvali]<0.05, grep("names|segmentCounts", names(Mpresence))]
				if (dim(Mpresence.sel)[1]>0) { 
					names(Mpresence.sel) <- gsub("segmentCounts.","", names(Mpresence.sel))

					library(reshape2)
					Mpresence.sel.t <- melt(Mpresence.sel)

					library(latticeExtra)
					p <- barchart(variable~value|names, data=Mpresence.sel.t)

					ofle <- paste(ofle,".pdf",sep="")
					plot.new(); pdf(ofle, h=8, w=8)
					print(p)
					dev.off()
				}
			}
		}

		print(sessionInfo())

		print(date())
		#	save.image(paste(basename(bedplusf),".Rdata",sep=""))
		#load("./MT_coveredCG.Ratio0.5_bs20_DE.txt.gz.Rdata")
	}

test.getDNAshapeVectorInRegion <- 
	function() {

		source("~/Projects/biter_biter_shared/bin/auxTFBSTools.R")

		# fasta
		fle <- "tmp.fa"
		caller.getDNAshapeVectorInRegion(fle, ofleTag="haha")
		fle <- "FOXA2.train.labels.tsv.gz"
		readBigTable(fle, nrows=10, verbose=T, header=T, col.names=c("chr","start","end","d1","d2","d3"))
		caller.expressedTxRegionsOrSeqs(fle, nrows=10, verbose=T, header=T, col.names=c("chr","start","end"))[[1]] 
#		caller.getDNAshapeVectorInRegion(fle, ofleTag="haha",  gnm="hg19", regions=c("self"), nrows=10, verbose=T, header=T)
		# TODO ? how about 
		# test.caller.getSeqFeatures

}
	

# TODO
caller.getDNAshapeVectorInRegion <-
	function(fle, dss=NULL, as.fasta=F,
					 featureType=c("1-mer", "1-shape"),
					 ofleTag=NULL, mc.cores=1, ...) {
		if (is.null(dss)) {
			if (as.fasta) {
				dss <- readDNAStringSet(fle)
			} else {
				print(list(...))
				dss <- caller.expressedTxRegionsOrSeqs(fle, ...)[[1]]
			}
		}
		print("head(dss)")
		print(head(dss))

#		cnts.tfs <- getSiteCountswMatrixAndDNAStringSet(, dss, ...)
		library(DNAshapeR);
		pred <- getShape(fle, "All", parse=T)
		pred <- getShape(fle)
		featureVector.normT <- encodeSeqShape(fle, pred, featureType, normalize=T)
		featureVector.normF <- encodeSeqShape(fle, pred, featureType, normalize=F)

		if (! is.null(ofleTag)) {
			writeData(featureVector.normF, paste(ofleTag, ".normF.txt.gz", sep=""), row.names=T)
			writeData(featureVector.normT, paste(ofleTag, ".normT.txt.gz", sep=""), row.names=T)
		}

		print(sessionInfo())

# regression TODO
#function() <- regressor(fle, dss=NULL, featureType=c("1-mer", "1-shape"), otag="DNAshapeR", genome="mm10", ...) {
	#plotShape(pred$MGW)
	#heatShape(pred$ProT, base <- length)
	# This is for linear regression and its correlation with the ground truth
	#fn4 <- system.file("extdata", "SELEXsample_short.s", package = "DNAshapeR")
	#experimentalData <- read.table(fn4)
	#df <- data.frame(affinity=experimentalData$V1, featureVector)
	#
	#library(caret)
	#trainControl <- trainControl(method = "cv", number = 3, savePredictions = TRUE) 
	#model <- train(affinity~ ., data = df, trControl=trainControl, method="lm", preProcess=NULL) 
	#model
}

test.getMotifCountsInRegion <- 
	function() {
		source("~/Projects/biter_biter_shared/bin/auxTFBSTools.R")
		# test sequence extraction
		#fle <- "~/PI_HOME/Data/casco/UCSC_tracks/mm10/canonical/promoters.bed.gz"
		fle <- "~/PI_HOME/Data/casco/UCSC_tracks/mm10/a.bed.gz"; region="self"; regions.only=F
		matrixSource <- "JASPAR2014"; expressedTx.colInd <- 1 #mcols of GRanges
		caller.getMotifCountsInRegion(fle, expressedTx.colInd=1, region=region, regions.only=regions.only,
																	src=matrixSource, all_versions=T, matrixtype="PWM",
																	ofleTag="dummy_") 
		caller.getMotifCountsInRegion(fle, expressedTx.colInd=1, region=region, regions.only=regions.only,
																	src=matrixSource, all_versions=T, matrixtype="PWM", gnm="hg19",
																	ofleTag="dummy_hg19_") 

		# test getting sequences from fasta 
		fle <- "~/PI_HOME/Data/casco/UCSC_tracks/mm10/promoterseqsself.seqs.gz"
		caller.getMotifCountsInRegion(fle, as.fasta=T, 
																	ofleTag="dummy_default_") 

		# test getting Matrix from file
		matrixSource <- "~/PI_HOME/Data/casco/Homer/custom.motifs"
		caller.getMotifCountsInRegion(fle, as.fasta=T, 
																	src=matrixSource, byrow=T, matrixtype="PWM",
																	ofleTag="dummy_mFile_", mc.cores=1) 

		fle <- "~/Projects/biter_biter_shared/results/2016-02-21/LM_ATACseq/TFBS/MTAG/JASPAR2014/RTAG/Acc4yAyQoQwSPMRBroadwSHIFT0/callpeak_peaks.broadPeak"
		matrixSource <- "JASPAR2014"; expressedTx.colInd <- 1; region="self"; regions.only=F;  #mcols of GRanges
		caller.getMotifCountsInRegion(fle, expressedTx.colInd=1, region=region, regions.only=regions.only,
																	src=matrixSource, all_versions=T, matrixtype="PWM",
																	ofleTag="dummy_") 
	}

caller.getMotifCountsInRegion <-
	function(fle, dss=NULL, as.fasta=F, 
					 ofleTag=NULL, mc.cores=1, ...) { 
		if (is.null(dss)) {
			if (as.fasta) {
				dss <- readDNAStringSet(fle)
			} else {
				dss <- caller.expressedTxRegionsOrSeqs(fle, ...)[[1]]
			}
		}
		print("head(dss)")
		print(head(dss))

		MatrixList <- getMatrixList(...)
		print("MatrixList[[1]]")
		print(MatrixList[[1]])
		print("MatrixList")
		print(MatrixList)

		cnts.tfs <- getSiteCountswMatrixAndDNAStringSet(MatrixList, dss, ...)

		if (! is.null(ofleTag)) 
			writeData(cnts.tfs, paste(ofleTag, ".txt.gz", sep=""), row.names=T)

		print(sessionInfo())

	}

test.caller.getExpressed <- 
	function() {
	source("./auxML.R")
	nc <- 1

	fle <- "~/Projects/biter_lingliu_DNAdamage/results/2015-10-12/featureCounts.txt.gz"
	ofle.tag="Liu_Cheung"
	spike <- list(id="ERCC", rm=F)
	param.preprocess <-  list(col.sel.pat=".sorted.bam$", colname.split.rank=2, colname.split.pat=".")
	dta.exp <- caller.getExpressed(fle, param.preprocess=param.preprocess, ofle.tag=ofle.tag, nc=nc, row.names="Geneid")

	fle <- "~/Projects/biter_jbrett_DNAmethylation/results/2015-09-14/FC_Brett_2015_RNAseq/featureCounts.txt.gz"
	ofle.tag <- "Brett_2015_RNAseq"
	spike <- list(id="ERCC", rm=F)
	param.preprocess <-  list(col.sel.pat=".sorted.bam$", colname.split.rank=2, colname.split.pat=".")
	dta.exp <- caller.getExpressed(fle, param.preprocess=param.preprocess, ofle.tag=ofle.tag, nc=nc, spike=spike, row.names="Geneid")

	fle <- "./_FCaging_wlongFragments_cutoff5/featureCounts.txt.gz"
	param.preprocess <- list(col.sel.pat='.sorted.bam$', colname.split.rank=1, colname.split.pat='.')
	spike <- list(id='ERCC', rm=F)
	param.process <- list(method="filterRowMeanMin", cutoff=10) 
	dta.exp <- caller.getExpressed(fle, param.preprocess=param.preprocess, ofle.tag='featureCounts', nc=8, spike=spike, row.names='Geneid', param.process=param.process)

	fle <- "./featureCounts.txt.gz"
	param.preprocess <- list(col.sel.pat='.srt.bam$', colname.split.rank=1, colname.split.pat='.', scale=T)
	spike <- list(id='ERCC', rm=T)
	param.process <- list(method='filterRowMeanMin', cutoff=0)
	dta.exp <- caller.getExpressed(fle, param.preprocess=param.preprocess, ofle.tag='featureCounts_scaled', nc=nc, spike=spike, row.names='Geneid', param.process=param.process)

	fle <- "/home/biter/Projects/biter_biter_shared/results/2016-05-04/PF_Soleimani_2012_DevCell_ChIPseq/calcsim/genome.chr19.bedcov.gz"
	param.preprocess <- list(col.sel.pat='.sorted.bam$', colname.split.rank=1, colname.split.pat='.')
	spike <- list(id='ERCC', rm=T)
	param.process <- list(method='filterRowMeanMin', cutoff=0)
	dta.exp <- caller.getExpressed(fle, dta=dta, param.preprocess=param.preprocess, ofle.tag='featureCounts', nc=nc, spike=spike, param.process=param.process, density.norm.colname=NULL)

	param.preprocess <- list(col.sel.pat='.sorted.bam$', colname.split.rank=2, colname.split.pat='.')
	spike <- list(spikeId='ERCC', rm=T)
	param.process <- list(method='filterRowMeanMin', cutoff=10) 
	dta.exp <- caller.getExpressed(fle='genome.chr19.bedcov.gz', param.preprocess=param.preprocess, ofle.tag='similarityW1000S1000MRC10', nc=1, spike=spike, density.norm.colname=NULL, param.process=param.process)

	#20170306 Brett_2015_RNAseq
	param.preprocess <- list(col.sel.pat='.nsrt.bam', colname.split.rank=2, colname.split.pat='.', mk.uniq.names=T, rm.col.sel.pat='.', sort.names=T)
	spike <- list(spikeId='ERCC', rm=T)
	dta.exp <- caller.getExpressed('2featureCounts.txt.gz', param.preprocess=param.preprocess, ofle.tag='2featureCounts', nc=1, spike=spike, row.names='Geneid')

	#20170503 Liu_2017_ATACseq
	param.preprocess <- list(col.sel.pat='.nsrt.bam', colname.split.rank=2, colname.split.pat='.', mk.uniq.names=T, sort.names=T)
	fle <- 'featureCounts.txt.gz'
	ofle.tag <- "features"
	spikeId <- NA
	dta.exp <- caller.getExpressed(fle, ofle.tag=NULL, nc=1, spikeId=spikeId, density.norm.colname=NULL, param.process=param.process)
}

# returns count data of expressed transcripts
# GMM fit to log(counts) does not always work; fit GMM to read DENSITY! (i.e.count/length)
caller.getExpressed <- 
	function(fle, dta=NULL, param.preprocess=NULL, nc=1, spikeId=NA, rep.sub='_\\d+$',
					 ofle.tag="featureCounts", density.norm.factor=NULL, density.norm.colname="Length", 
					 param.process=list(method="GMM"), ...) {
		ggplot <- set.ggplot() 

		if (is.null(dta))
			dta <- readBigTable(fle, header=T, ...)
		#dta <- readBigTable(fle, header=T, row.names="Geneid")

		if (is.null(density.norm.factor)) 
			if (is.null(density.norm.colname)) {
				density.norm.factor <- data.frame(Length=dta$End - dta$Start)
			} else {
				density.norm.factor <- dta[,density.norm.colname, drop=F]
			}

		if (!is.null(param.preprocess)) 
			dta <- preprocess(dta=dta, param=param.preprocess, verbose=T)

		# plot spike-in-rate
		if (! is.na(spikeId)) {
			# plot spikeIn Ratios if not they are removed
			spikei <- grepl(spikeId, rownames(dta))
			if (any(spikei)) {
				spikeInp <- data.frame(spikeInp=t(t(colSums(dta[spikei,]) / colSums(dta) * 100)))
				spikeInp$group <- sub(rep.sub, '', rownames(spikeInp), perl=T)
				spikeInp$Id <- sub(paste("(",rep.sub,")",sep=""), '\\1', rownames(spikeInp), perl=T)
				spikeInp$Id <- rownames(spikeInp)
				print(spikeInp)
				title <- paste(ofle.tag, ";", "N=", nrow(dta), sep="")
				ggplot(spikeInp, aes(x=Id, y=spikeInp, color=group)) + geom_bar(stat="identity", fill="white", size=2) +
				labs(x="", y="Spike-in reads wrt total(%)", title=title)
				ggsave(paste(ofle.tag,"_spikep.pdf",sep=""))
			}
		}

		# plot raw data
		param.EDApreprocess <- list(sort.names=T, ps=1, logf='log')
		get.exploratoryAnalysisPlots(dta=dta, ofle.tag=paste(ofle.tag, "raw",sep=""),
																 rep.sub=rep.sub,
																 param.preprocess=param.EDApreprocess)
																 #param.preprocess=param.EDApreprocess, gcol=group)

		if (! is.null(ofle.tag)) {
			if (param.process$method == "GMM") {
				plot.new(); pdf(paste(ofle.tag,"_",param.process$method,".pdf",sep=""),h=8,w=8) 
			} else if (param.process$method == "filterRowMeanMin") {
				plot.new(); pdf(paste(ofle.tag, "_",param.process$method,"_",param.process$cutoff,".pdf",sep=""),h=8,w=8) 
			} else {
				stop(param.process$method, " is not implemented")
			}
		}

		# classify expressed and unexpressed using GMM2
		if (param.process$method == "GMM") {
			cl <- 1
			res <- mclapply(1:ncol(dta), 
											function(cl) {
												message("Doing ", cl, " ", colnames(dta)[cl])
												# TODO FIXME isn't it log2(X)/length(X) ??? 
												rw <- dta[,cl] > 0 # non-zero 
												dta.exp <- log2((1+dta[rw,cl, drop=F]) / density.norm.factor[rw,])
												GMM(dta.exp, title=names(dta)[cl], verbose=T)
											},
											mc.cores=1) #Does not do parallel plotting FIXME
#print("res:")
#print(head(res))
			expTxIds <- Reduce(union, unlist(res))
		} else if (param.process$method == "filterRowMeanMin") {
			expTxIds <- rownames(dta)[rowMeans(dta) > param.process$cutoff]
		} else {
			stop(param.process$method, " is not implemented")
		}
#print("expTxIds:")
#print(head(expTxIds))
#print("dta:")
#print(head(dta))
		dta.exp <- dta[as.numeric(expTxIds),]

		message("Number of expressed genes=", length(expTxIds))
		print("Expressed-gene-count-summary:")
		print(summary(dta.exp)) 

		if (! is.null(ofle.tag)) {
			dev.off()
			if (param.process$method == "GMM") {
				ofle <- paste(ofle.tag,"_",param.process$method,"2expressed.txt.gz",sep="")
			} else if (param.process$method == "filterRowMeanMin") {
				ofle <- paste(ofle.tag, "_",param.process$method,"_",param.process$cutoff,"expressed.txt.gz",sep="")
			} else {
				stop(param.process$method, " is not implemented")
			}
			writeData(dta.exp, ofle, row.names=T)
		}

		group <- dta.exp[,1, drop=F]
		group[,1] <- NULL
		if ( ! is.na(spikeId) ) {
			group <- data.frame(metadata=ifelse(grepl(spike$id, rownames(dta.exp)), "spikeIn", "Trx"))
		}

		# TODO add metadata to plot spikeIn info
		if (param.process$method == "GMM") {
			ofle.tag <- paste(ofle.tag,"_",param.process$method,"2expressed_EDA",sep="")
		} else if (param.process$method == "filterRowMeanMin") {
			ofle.tag <- paste(ofle.tag, "_",param.process$method,"_",param.process$cutoff,"expressed_EDA",sep="")
		} else {
			stop(param.process$method, " is not implemented")
		}

		get.exploratoryAnalysisPlots(dta=dta.exp, ofle.tag=ofle.tag,
																 rep.sub=rep.sub,
																 param.preprocess=param.EDApreprocess, gcol=group)

		get.exploratoryAnalysisPlots(dta=dta.exp, ofle.tag=paste(ofle.tag, "scaled",sep=""),
																 rep.sub=rep.sub,
																 param.preprocess=c(scale=T, param.EDApreprocess), gcol=group)
		print(sessionInfo())
		dta.exp
	}

preprocess.4debug <-
	function(fle, dta=NULL, param=NULL, verbose=F) {
		print("in auxML.R")

		if (is.null(dta))
			dta <- readBigTable(fle, header=T)

    if (verbose) print(head(dta,2))

		dta %>% filter_(dta, param$subset.cond)

	}


test.caller.DESeq <- 
	function() {
		source("~/Projects/biter_biter_shared/bin/auxRNAseq.R")
		padj=0.01
		condA <- "WT";
		design <- rbind(data.frame(Sample=paste(condA,1:5,sep="_"), Group=condA),
										data.frame(Sample=paste("KO",1:5,sep="_"), Group="KO"))
		dta.exp <- data.frame(matrix(as.integer(rnorm(100000,100,5)),ncol=10))
		names(dta.exp) <- design$Sample
		a <- cbind(data.frame(id=1:10000, len=as.integer(rnorm(10000,1000,100))), dta.exp)
		a <- dta.exp
		head(a)
		gid <- NULL
		ofleTag <- NULL 
		caller.DESeq(dta.exp, mdta, gid=gid, condA=condA, padj=padj, ofleTag=ofleTag, debug=debug)

	}

#Uses DEseq for calling differentially expressed genes
#Expression values are normalized using geometric mean of the genes as a reference sample
#New variance stabilizing transformation feature is used
#Variance in designs wo/replicates is fit by ignoring the sample labels and outliers are not
caller.DESeq <- 
	function(a, design, gid=NULL, condA="siUntreated", padj=0.2, fitType="parametric", ofleTag=NULL, debug=F) 
	{
		library(DESeq)
		library("RColorBrewer")
		library("gplots")

		#filter out not expressed genes
		all <- a[rowSums(a)>0,];

		#dimensions of data
		ngenes   <- dim(all)[1]
		nsamples <- dim(all)[2]

		print(dim(all))
		print(head(all,2))
		print(design)

		print(length(unique(design$Group)))
		print(length(design$Group));

		if (debug) 
			save(list=ls(), file="debugme.RData", envir=environment())

		cds <- newCountDataSet(all,design$Group)
		cds <- estimateSizeFactors(cds)

		print(cds)

		if (length(unique(design$Group)) == length(design$Group)) { #without replicates
			if (fitType == "local")
				#local dispersion fit as in the paper
				cds <- estimateDispersions(cds, method="blind", sharingMode="fit-only", fitType="local") 
			else
				#parametric dispersion fit: new method
				cds <- estimateDispersions(cds, method="blind", sharingMode="fit-only")
		} else {
			if (fitType == "local")
				#local dispersion fit as in the paper
				cds <- estimateDispersions(cds,fitType="local")
			else
				#parametric dispersion fit: new method
				cds <- estimateDispersions(cds)
		}

		print(str(fitInfo(cds)))

		if (! is.null(ofleTag)) {
			plot.new(); pdf(paste(ofleTag,"_",condA,".pdf", sep=""),h=10,w=10) 
		}
		plotDispEsts(cds, main="Dispersion estimates")
		#sizeFactors(cds)
		for (condB in unique(design$Group)) {
			if (condB == condA) {
				all_norm <- counts(cds, normalized=TRUE) 
				hmcol = colorRampPalette(brewer.pal(9, "GnBu"))(100)
				#		heatmap( as.matrix( dist(t(all)) ), symm=TRUE, main="raw" )
				heatmap.2( as.matrix( dist(t(all_norm)) ), main="normalized", col = hmcol, trace="none", margin=c(10,6))
			} else {
				#vsd <- getVarianceStabilizedData( cds )
				res <- nbinomTest (cds, condA, condB)
				print(plotMA(res, main=paste(condB,"vs",condA)))
				print(hist(res$pval, breaks=100, col="skyblue", border="slateblue", main=paste(condB,"vs",condA)))
				print(hist(res$padj, breaks=100, col="skyblue", border="slateblue", main=paste(condB,"vs",condA)))

				if (is.null(gid)) {
					deseq_DE <- res;
				} else {
					deseq_DE <- merge(res, gid, by="id")
				}
				# order based on padj
				deseq_DE <- deseq_DE[order(deseq_DE$padj),]
				ofle <- paste(ofleTag,"_", condB,"vs", condA,"_DE.txt",sep="")
				writeData(deseq_DE, ofle=ofle, row.names=F,
									commentLine=paste("A=", condA, " B=", condB, " foldChange=B/A", sep=""))
			}
		}
		if (! is.null(ofleTag))
			dev.off()

		deseq_DE

	}


# TODO integrate the rest
test.caller.DE.edgeR <- 
	function() {
		source("./auxggplot2.R")

		# Mike's data
		#		fle <- "Wosczyna_2015_PullDown_wospikeIn_none_count.txt.gz"
		#		tag <- sub('(.*)_count.txt.gz', '\\1', fle)
		#		rep.sub <- '_\\d+$'
		#		param <- list(out="top", method="none", contrast=c(1,-1, -1, 1), FDRcut=0.05, species="Mm")
		#		param <- list(out="top", method="TMM", contrast=c(1,-1, -1, 1), FDRcut=0.01, species="Mm")
		#		plot2file <- T
		#		top <- caller.DE.edgeR(fle=fle, rep.sub=rep.sub, param=param, tag=tag, plot2file = plot2file)

		fle <- "~/Projects/biter_lingliu_DNAdamage/results/2015-10-12/EG_Liu_2016_RNAseq/featureCounts_GMM2expressed.txt.gz"
		tag <- sub('/EG_(.*)/', '\\1', fle)
		tag <- "Liu_2015_RNAseq"
		rep.sub <- '_(\\d+)$'
		param <- list(out="top", scale="none", contrast=c(1,-1), FDRcut=0.05, species="Mm", tag="dummy")
		plot2file <- F
		top <- caller.DE.edgeR(fle=fle, rep.sub=rep.sub, param=param, tag=tag, plot2file = plot2file)

	}

#CURR EDA plot call for normalized counts
caller.DE.edgeR <-
	function(set, group, fle=NULL, rep.sub='\\d+$', param=NULL, tag=NULL, plot2file=F) {

		print(paste("Doing caller.DE.edgeR param:", paste(param, collapse=" "), sep=""))

		if (plot2file) {
			pdf(paste(tag, "_DE_edgeR.pdf", sep=""), h=8, w=8)
		}

		ggplot <- set.ggplot()

		if (!is.null(fle)) {
			dta <- readBigTable(fle, header=T, row.names="Id")
			column.labels <- colnames(dta) 
			group <- sub(rep.sub, '', column.labels, perl=T)
			set <- newSeqExpressionSet(as.matrix(dta),
																 phenoData = data.frame(group, row.names=column.labels))
		}

		ret <- NULL
		# If multiple groups are present and topNinvert is set, get negative controls by correlation to total library size
		if (param$out == "topNinvert" && length(unique(group)) > 2) {
			print("Correlation based negative control gene selection")
			# TODO add scaling based colSums
			dta <- counts(set)
			if (!any(is.na(normCounts(set))))
				dta <- normCounts(set)
			cs <- colSums(dta)

			# negative controls are at the end; correlations are sorted in acsending order
			top <- data.frame(cor=sort(apply(dta, 1, function(x) { cor(cs, x) })))
			hist(top[,1], n=20)

			top$negCon <- c(rep(F, param$topN), rep(T, nrow(dta) - param$topN ))

			if (any(grep("species", names(param))) && !is.na(param$species))
				top <- annotate(top, param$species)

			writeData(top, ofle=paste(tag, "_RUVg_negCon.txt.gz", sep=""), 
								row.names=T) 

			# reorder
			top <- top[rownames(dta),] 
			return(top[,2])
		}

		group.factor <- as.factor(group)

		N.factors <- length(varLabels(set@phenoData))

		# TODO generalize; the following does not work
		#frml <- paste("~0", paste(group.factor, collapse="+"), sep = "+")
		#frml <- paste(frml, paste0(varLabels(set@phenoData)[-1], collapse="+"), sep="+")
		#print(frml)
		#design <- model.matrix(as.formula(frml), data=pData(set))
		if (N.factors == 1) {
			design <- model.matrix(~0+group.factor, data=pData(set))
		} else if (N.factors == 2) {
			design <- model.matrix(~0+group.factor+W_1, data=pData(set))
		} else if (N.factors == 3) {
			design <- model.matrix(~0+group.factor+W_1+W_2, data=pData(set))
		} else {
			message("ERROR: N.factors >3")
			design <- model.matrix(~0+group.factor+W_1+W_2+W_3, data=pData(set))
		}
		colnames(design) <- gsub("group.factor", "", colnames(design))

		print(design)

		y <- DGEList(counts=counts(set), group=group)
		y <- calcNormFactors(y, method=param$scale)
		y <- estimateGLMCommonDisp(y, design)
		#y <- estimateGLMTrendedDisp(y, design) # added to vignette code
		y <- estimateGLMTagwiseDisp(y, design)

		fit <- glmFit(y, design) # neg binomial fit w EdgeR
		#fit <- glmQLFit(y, design, robust=T) # replaced with vignette code 

		if (param$out == "res") { #residuals
			ret <- residuals(fit, type="deviance")
		} else if (param$out == "topNinvert" || param$out == "top") {
			if (any(grep("contrast", names(param)))) {
				contrast <- c(param$contrast)
			} else if (length(unique(group)) == 2) {
				contrast <- c(1,-1) #when RUVg for two groups
			} else {
				stop("ERROR: param$contrast is not set for >2 groups")
			}

			if (N.factors > 1) {
				contrast <- c(param$contrast, 0*(1:(N.factors-1)))
			} else if (N.factors > 3) {
				contrast <- c(param$contrast, 0*(1:3))
			}
			print(contrast)

			lrt <- glmLRT(fit, contrast=contrast) # Likelihood ratio test
			top <- topTags(lrt, n=nrow(set))$table
			if (param$out == "topNinvert") {
				ret <- !(rownames(set) %in% rownames(top)[1:param$topN])
			} else if (param$out == "top") {
				#plotQLDisp(fit)
				print(plotBCV(y))		
				plotMDS(y) # noprint
				#plotMDS(y,col=as.numeric(targets$Genotype)) # noprint
				title <- paste(tag, paste(param,collapse=" "), 
											 paste(colnames(lrt$design), collapse="\n"), sep="\n")

				# TODO plot as layers
				#				http://stackoverflow.com/questions/15706281/controlling-order-of-points-in-ggplot2-in-r

				top$FDRcut <- param$FDRcut # patch for ggplot
				p <- ggplot(top, aes(x=logCPM, y=logFC, size=-log10(FDR), color = FDR < FDRcut)) + 
					geom_hline(yintercept=0, colour="gray", linetype=2) + 
					geom_point(size=1, alpha=0.5) + 
					geom_smooth() + labs(title=title)
				top$FDRcut <- NULL
				print(p)

				p <- ggplot(top, aes(x=logFC, y=-log10(FDR), color=logCPM>8)) + # && logFC > 5 && logFC < -5)) + 
					geom_hline(yintercept=2) + 
					geom_vline(xintercept=c(-2,2)) + 
					geom_point(size=1, alpha=0.5) + 
					labs(title=title)
				print(p)

				print(paste("DE at FDR<", param$FDRcut, sep=""))
				print(summary(dt <- decideTestsDGE(lrt, p.value=param$FDRcut)))

				if (any(grep("species", names(param))) && !is.na(param$species)) 
					top <- annotate(top, param$species)

				commentLine <- gsub("\n",";", title)

				writeData(top, ofle=paste(tag, "_DE_edgeR.txt.gz", sep=""), 
									commentLine=commentLine, row.names=T)

				if (any(grep("species", names(param))) && !is.na(param$species)) {
					GSEA(top, FDRcut=param$FDRcut, keytype="ENTREZID", type="goana", species=param$species, ofle.tag=tag) 
					GSEA(top, FDRcut=param$FDRcut, keytype="ENTREZID", type="kegga", species=param$species, ofle.tag=tag) 
				}

				ret <- top
			}
		} else {
			stop(param$out, " is not defined in function caller.DE.edgeR")
		}

		if (plot2file) 
			dev.off()

		ret
	}


# TODO write caller.GSEA for annotatation of coordinate files files and GSEA check
test.GSEA <- function() {
	source("./auxggplot2.R")

	d <- "~/Projects/biter_lingliu_DNAdamage/results/2015-10-12/EG_Wosczyna_2015_RNAseq/" 
	fle <- paste0(d, "/", "featureCounts.txt.gz")
	efle <- paste0(d, "/", "")
	param <- list(out="top", scale="upper", contrast=c(1,0,-1,0), FDRcut=0.05, species="Mm")
	plot2file <- T
	top <- caller.DE.edgeR(fle=fle, efle=efle, param=param, tag=tag, plot2file = plot2file)
	type <- "goana"
	type <- "kegga"
	species <- "Mm"

	# test for species not set
	GSEA(top, FDRcut=0.05, keytype="ENTREZID", type=type) 
	# test for wrong species
	GSEA(top, FDRcut=0.05, keytype="ENTREZID", type=type, species="Hs") 
	# test for correct species
	GSEA(top, FDRcut=0.05, keytype="ENTREZID", type=type, species=species, ofle.tag=tag) 

	topGSEAwpcut <- GSEA(top, FDRcut=0.05, keytype="ENTREZID", type=type, species=species) 

}

# TODO test for top=GLMfit object
# universe is restricted when FDRcut != NULL and universe is not set
GSEA <- 
	function(top, FDRcut=NULL, universe=NULL, keytype="ENTREZID",  
					 type="goana", species=NA, pcutoff=0.05, ofle.tag=NULL, valuetype="SYMBOL",
					 title="") {
		#	library(KEGGREST)
		#	listDatabases()
		#	org <- keggList("organism")
		#	hsa <- keggList("hsa")
		#	query <- keggGet(c("hsa:10458"))
		if (is.na(species)) {
			message("species is required for GSEA")
			return()
		}

		ggplot <- set.ggplot()

		if (! is.null(FDRcut)) {
			unambig <- !duplicated(top[,keytype])
			if (is.null(universe))
				universe <- top[unambig, keytype]
			# reduce top to a list
			topl <- top[unambig & top$FDR < FDRcut, keytype]
		}

		# set functions
		f <- get(type) #kegga or goana
		if (type == "goana") {

			topGSEAf <- "topGO"
			lib <- paste('org.',species,'Mm.eg.db')
			library(lib)
			GO2ALLEGS <- paste("org", species, "egGO2ALLEGS", sep = ".")
			EG <- AnnotationDbi::toTable(get(GO2ALLEGS))
		} else if (type == "kegga") {
			topGSEAf <- "topKEGG"
			species.KEGG <- switch(species, "Hs"="hsa", "Dm"="dme", "Mm"="mmu", "Rn"="rno", "Pt"="ptr")
			EG <- getGeneKEGGLinks(species.KEGG) 
		}
		tf <- get(topGSEAf)

		gsig <- f(topl, universe=universe, species = species)
		topGSEA <- tf(gsig, n=nrow(gsig)) # order categories

		topGSEAwpcut <- topGSEA[topGSEA$P.DE<pcutoff,]
		head(topGSEAwpcut)

		names(EG)[1:2] <- c("gid", "GSid")
		EG.topl <- EG[EG$gid %in% topl & EG$GSid %in% rownames(topGSEAwpcut), ]
		head(EG.topl)
		topGSEAwpcut$DEgenes <- as.character(lapply(rownames(topGSEAwpcut), 
																								function(GSid) { 
																									paste(sort(top[top[,keytype] %in% unique(EG.topl[EG.topl$GSid %in% GSid,"gid"]),valuetype]), collapse="|") 
																								} ))

		if (is.null(ofle.tag)) {
			n <- min(nrow(topGSEAwpcut), 30)
			print(topGSEAwpcut[1:n,])
		} else {
			writeData(topGSEAwpcut, ofle=paste(ofle.tag, "_", topGSEAf, ".txt.gz", sep=""), row.names=T)
		}

		if (type == "kegga") 
			topGSEA %<>% rename(Term=Pathway)
		print("topGSEA:")
		print(dim(topGSEA))
		print(nrow(topGSEA))
		print("topGSEAwpcut:")
		print(dim(topGSEAwpcut))
		print(nrow(topGSEAwpcut))
		nP.DE <- nrow(topGSEAwpcut)
		topGSEA$pcutoff <- pcutoff
		topGSEA$rank <- 1:nrow(topGSEA)
		topGSEA$label <- paste0(topGSEA$rank,":",topGSEA$Term)
		title <- paste0(title,"\nGSEA:",topGSEAf,";species:",species,";pcutoff=",pcutoff,";N=", nrow(topGSEA),";Nsig=",nP.DE)
#		p <- my_bin(topGSEA, aes(x=N, y=DE/N, color=P.DE < pcutoff, label=label),
#								textRepelFilter=paste0("rank < ", min(nP.DE,20)), smooth=F,
#								long=F, slope=0, intercept=-1, size=rel(2), bins=nrow(topGSEA)) + labs(title=title)
#		print(p)
		p <- my_bin(topGSEA, aes(x=log(N), y=DE/N, color=P.DE < pcutoff, label=label),
								textRepelFilter=paste0("rank < ", min(nP.DE,20)), smooth=F,
								long=F, slope=0, intercept=-1, size=rel(2), bins=nrow(topGSEA)) + labs(title=title)
		print(p)

		top$pcutoff <- NULL
		top$rank <- NULL
		top$label <- NULL

		topGSEAwpcut

	}

test.caller.normalize <- 
	function() {
		source("~/Projects/biter_biter_shared/bin/auxggplot2.R");

		# Liu_2015_RNAseq
		print("---------------------------------------------------------------")
		spike <- NA
		rep.sub <- '_\\d+$'
		fle <- "/home/biter/Projects/biter_lingliu_DNAdamage/results/2015-10-12/FC_Liu_2015_RNAseq/featureCounts.txt.gz"
		#tag <- sub("FC_", "", grep("FC_", strsplit(fle, "/")[[1]], value=T)) 
		tag <- "Liu_2014_RNAseq_Cheung_2013_RNAseq"
		print(paste("Doing ", tag, sep=""))
		caller.normalize(tag,fle=fle,rep.sub=rep.sub,spike=spike, species="Mm", contrast=c(1,-1))
		print("---------------------------------------------------------------")

		# Liu_2014_RNAseq_Cheung_2013_RNAseq
		print("---------------------------------------------------------------")
		spike <- NA
		rep.sub <- '_\\d+$'
		fle <- "/home/biter/Projects/biter_lingliu_DNAdamage/results/2015-10-12/featureCounts.txt.gz"
		#tag <- sub("FC_", "", grep("FC_", strsplit(fle, "/")[[1]], value=T)) 
		# add DE gene names to the GSEA table
		GO2ALLEGS <- paste("org", species, "egGO2ALLEGS", sep = ".")
		EG.GO <- AnnotationDbi::toTable(get(GO2ALLEGS))
		head(EG.GO)
		EG.GO.topl <- EG.GO[EG.GO$gene_id %in% topl & EG.GO$go_id %in% rownames(topGSEAwpcut), ]
		topGSEAwpcut$DEgenes <- as.character(lapply(rownames(topGSEAwpcut), 
																								function(go_id) { 
																									paste(sort(top[top[,keytype] %in% unique(EG.GO.topl[EG.GO.topl$go_id %in% go_id,"gene_id"]),valuetype]), collapse="|") 
																								} ))

		if (is.null(ofle.tag)) {
			n <- min(nrow(topGSEAwpcut), 30)
			print(topGSEAwpcut[1:n,])
		} else {
			writeData(topGSEAwpcut, ofle=paste(ofle.tag, "_", topGSEAf, ".txt.gz", sep=""), row.names=T)
		}

		# volcano plot to see if there's any geneset size bias
		topGSEA$pcutoff <- pcutoff
		p <- ggplot(topGSEA, aes(x=log2(N), y=DE/N, color=P.DE < pcutoff)) + # && logFC > 5 && logFC < -5)) + 
			geom_point(size=2, alpha=0.5) + labs(title=paste(topGSEAf, "analysis: pcutoff=",pcutoff, sep=""))
		print(p)
		top$pcutoff <- NULL

		topGSEAwpcut

	}

test.caller.normalize <- 
	function() {
		source("~/Projects/biter_biter_shared/bin/auxggplot2.R");

		# Liu_2015_RNAseq
		print("---------------------------------------------------------------")
		spike <- NA
		rep.sub <- '_\\d+$'
		fle <- "/home/biter/Projects/biter_lingliu_DNAdamage/results/2015-10-12/FC_Liu_2015_RNAseq/featureCounts.txt.gz"
		#tag <- sub("FC_", "", grep("FC_", strsplit(fle, "/")[[1]], value=T)) 
		tag <- "Liu_2014_RNAseq_Cheung_2013_RNAseq"
		print(paste("Doing ", tag, sep=""))
		caller.normalize(tag,fle=fle,rep.sub=rep.sub,spike=spike, species="Mm", contrast=c(1,-1))
		print("---------------------------------------------------------------")

		# Liu_2014_RNAseq_Cheung_2013_RNAseq
		print("---------------------------------------------------------------")
		spike <- NA
		rep.sub <- '_\\d+$'
		fle <- "/home/biter/Projects/biter_lingliu_DNAdamage/results/2015-10-12/featureCounts.txt.gz"
		#tag <- sub("FC_", "", grep("FC_", strsplit(fle, "/")[[1]], value=T)) 
		tag <- "Liu_2014_RNAseq_Cheung_2013_RNAseq"
		print(paste("Doing ", tag, sep=""))
		caller.normalize(tag,fle=fle,rep.sub=rep.sub,spike=spike, species="Mm", contrast=c(1,-1,0,0,0))
		print("---------------------------------------------------------------")

		# Liu_2014_RNAseq
		print("---------------------------------------------------------------")
		spike <- NA
		rep.sub <- '_\\d+$'
		fle <- "/home/biter/Projects/biter_lingliu_DNAdamage/results/2015-10-12/FC_Liu_2014_RNAseq/featureCounts.txt.gz"
		tag <- sub("FC_", "", grep("FC_", strsplit(fle, "/")[[1]], value=T)) 
		print(paste("Doing ", tag, sep=""))
		caller.normalize(tag,fle=fle,rep.sub=rep.sub,spike=spike, species="Mm", contrast=c(0,-1,0,1))
		print("---------------------------------------------------------------")

		# Cheung_2013_RNAseq
		print("---------------------------------------------------------------")
		spike <- NA
		rep.sub <- '_\\d+$'
		fle <- "/home/biter/Projects/biter_lingliu_DNAdamage/results/2015-10-12/FC_Cheung_2013_RNAseq/featureCounts.txt.gz"
		tag <- sub("FC_", "", grep("FC_", strsplit(fle, "/")[[1]], value=T)) 
		print(paste("Doing ", tag, sep=""))
		caller.normalize(tag,fle=fle,rep.sub=rep.sub,spike=spike, species="Mm", contrast=c(1,-1))
		print("---------------------------------------------------------------")

		# Wosczyna_2016_RNAseq
		print("---------------------------------------------------------------")
		spike <- NA
		rep.sub <- '_\\d+$'
		fle <- "/home/biter/Projects/biter_wosczyna_miRNAsInMuscleCellFate/results/2015-07-17/FC_Wosczyna_2015_RNAseq/featureCounts.txt.gz"
		tag <- sub("FC_", "", grep("FC_", strsplit(fle, "/")[[1]], value=T)) 
		print(paste("Doing ", tag, sep=""))
		caller.normalize(tag,fle=fle,rep.sub=rep.sub,spike=spike, species="Mm", contrast=c(1,-1))
		print("---------------------------------------------------------------")

		# Wosczyna_2015_PullDown
		print("---------------------------------------------------------------")
		spike <- list(spikeId="ERCC", spikeInp=10)
		rep.sub <- '_\\d+$'
		fle <- "/home/biter/Projects/biter_wosczyna_miRNAsInMuscleCellFate/results/2015-07-17/FC_Wosczyna_2015_PullDown/featureCounts.txt.gz"
		tag <- sub("FC_", "", grep("FC_", strsplit(fle, "/")[[1]], value=T)) 
		print(paste("Doing ", tag, sep=""))
		caller.normalize(tag,fle=fle,rep.sub=rep.sub,spike=spike, species="Mm", contrast=c(1,-1,-1,1))
		print("---------------------------------------------------------------")

		# Brett_2015_RNAseq
		print("---------------------------------------------------------------")
		spike <- list(spikeId="ERCC", spikeInp=10)
		rep.sub <- '_\\d+$'
		fle <- "/home/biter/Projects/biter_jbrett_DNAmethylation/results/2015-09-14/FC_Brett_2015_RNAseq/featureCounts.txt.gz"
		tag <- sub("FC_", "", grep("FC_", strsplit(fle, "/")[[1]], value=T)) 
		print(paste("Doing ", tag, sep=""))
		caller.normalize(tag,fle=fle,rep.sub=rep.sub,spike=spike, species="Mm", contrast=c(1,0,0,-1,0,0,0,0,0), removeSpikes=F)
		print("---------------------------------------------------------------")

		# Zebrafish
		print("---------------------------------------------------------------")
		spike <- list(spikeId="ERCC", spikeInp=10)
		rep.sub <- '\\d+$'
		tag <- "zebrafishRNASeq"
		caller.normalize(tag,rep.sub=rep.sub,spike=spike, contrast=c(-1,1))
		print("---------------------------------------------------------------")


	}


test.caller.profileNormalization.DE <- 
	function() {
		# 
		source("~/Projects/biter_biter_shared/bin/auxggplot2.R");
		print("---------------------------------------------------------------")
		spike <- list(id="ERCC", rm=F)
		gcol <- spike
		fle <- "/home/biter/Projects/biter_lingliu_DNAdamage/results/2015-10-12/all_Liu_RNAseq_featureCounts.txt.gz"
		tag <- sub("EG_", "", grep("EG_", strsplit(fle, "/")[[1]], value=T)) 
		tag <- "all_Liu_RNAseq_featureCounts"
		print(paste("Doing ", tag, sep=""))
		scale <- "none"
		scale <- "upper"
		species <- "Mm"
		nc <- 1
		param.DE <- list(method="edgeR", FDRcut=0.01, contrast=list(c(1,0,0,-1,0,0,0,0,0)))
		caller.profileNormalization.DE(tag,fle=fle, rep.sub=rep.sub, spike=spike, gcol=gcol, 
																	 species=species, param.DE=param.DE, scale=scale)
		print("---------------------------------------------------------------")

		# Zebrafish
		print("---------------------------------------------------------------")
		spike <- list(id="ERCC")
		gcol <- spike
		rep.sub <- '\\d+$'
		tag <- "zebrafishRNASeq"
		param.DE <- list(method="edgeR", FDRcut=0.01, contrast=list(c(-1,1)))
		scale <- "upper"
		species <- NA
		caller.profileNormalization.DE(tag, rep.sub=rep.sub, spike=spike, gcol=gcol, 
																	 species=species, param.DE=param.DE, scale=scale)
		print("---------------------------------------------------------------")

		# Brett_2015_RNAseq
		print("---------------------------------------------------------------")
		spike <- list(id="ERCC", rm=F)
		spike <- list(id="ERCC", rm=T)
		gcol <- spike
		fle <- "/home/biter/Projects/biter_jbrett_DNAmethylation/results/2015-09-14/EG_Brett_2015_RNAseq/featureCounts_GMM2expressed.txt.gz"
		tag <- sub("EG_", "", grep("EG_", strsplit(fle, "/")[[1]], value=T)) 
		print(paste("Doing ", tag, sep=""))
		scale <- "none"
		scale <- "upper"
		species <- "Mm"
		nc <- 1
		param.DE <- list(method="edgeR", FDRcut=0.01, contrast=list(c(1,0,0,-1,0,0,0,0,0)))
		caller.profileNormalization.DE(tag,fle=fle, rep.sub=rep.sub, spike=spike, gcol=gcol, 
																	 species=species, param.DE=param.DE, scale=scale)
		print("---------------------------------------------------------------")

		# vanVelthoven_2015_PullDownTrx2_2016_PullDownTrx
		print("---------------------------------------------------------------")
		spike <- list(id="ERCC", rm=T)
		gcol <- spike
		rep.sub <- '_\\d+$'
		fle <- "/home/biter/Projects/biter_cvanvelt_SCTrxTrl/data/vanVelthoven_2016/PullDownTrx/EG_FC_RM_TR_ada/featureCounts_GMM2expressed.txt.gz"
		tag <- sub("EG_", "", grep("EG_", strsplit(fle, "/")[[1]], value=T)) 
		print(paste("Doing ", tag, sep=""))
		scale <- "none"
		scale <- "upper"
		species <- "Mm"
		nc <- 1
		param.DE <- list(method="edgeR", FDRcut=0.01, contrast=list(c(0,0,0,1,0,-1)))
		param.DE <- list(method="edgeR", FDRcut=0.01, contrast=list(c(0,-1,0,1,0,0), c(0,0,0,1,0,-1)))
		caller.profileNormalization.DE(tag,fle=fle, rep.sub=rep.sub, spike=spike, gcol=gcol, 
																	 species=species, param.DE=param.DE, scale=scale)
		print("---------------------------------------------------------------")

		# LM20161025
		source("~/Projects/biter_biter_shared/bin/auxggplot2.R");
		print("---------------------------------------------------------------")
		spike <- list(id="ERCC", rm=T)
		rep.sub <- '_\\d+$'
		fle <- "./input4peakCountPerGene"
		tag <- paste0("LM")
		print(paste("Doing ", tag, sep=""))
		scale <- "none"
		scale <- "upper"
		param.DE <- list(method="edgeR", FDRcut=0.1, contrast=list(c(1,0,-1), c(0,1,-1)))
		caller.profileNormalization.DE(tag,fle=fle,rep.sub=rep.sub,spike=spike, species="Mm", param.DE=param.DE, scale=scale)

		print("---------------------------------------------------------------")
	}

# TODO check this for correlation matrix of trees
#http://www.sthda.com/english/wiki/hierarchical-clustering-essentials-unsupervised-machine-learning
caller.profileNormalization.DE <- 
	#	function(tag="zebrafishRNASeq", fle=NULL, rep.sub='\\d+$', id="ERCC", spikeInp=10, nc=1) {
	function(tag="zebrafishRNASeq", fle=NULL, dta=NULL, rep.sub='_(\\d+)$', spike=NA, nc=1, gcol=NA,
					 species=NA, param.DE=NA, scale=c("none", "upper"), k=1, ...) {
		# Risso et al Zebrafish test
		#ggplot <- set.ggplot(...)
		ggplot <- set.ggplot()

		# Read data
		if (tag=="zebrafishRNASeq") {
			library(zebrafishRNASeq)
			data(zfGenes)
			dta <- zfGenes
			filter <- apply(dta, 1, function(x) length(x[x>5])>=2)
			dta <- dta[filter,]
		} else {
			if (!is.null(fle)) 
				dta <- readBigTable(fle, header=T, row.names="Id") 
		}
		# make sure the data is sorted by colnames - useful for visualization
		dta <- as.data.frame(preprocess(dta=dta, param=list(sort.names=T), verbose=F))

		print(head(dta))

		# remove spikes if set so
		if ( is.list(spike) && any(grep("id", names(spike))) ) {
			if ( any(grep("rm", names(spike))) && spike$rm == T ) {
				spikei <- grepl(spike$id, rownames(dta))
				dta <- dta[!spikei,]
				print("All-gene-count-summary after spike removal:")
				print(summary(dta))
			}
		}

		# column group (eg. genes vs spikes etc)
		if ( is.list(gcol) && any(grep("id", names(gcol))) ) 
			gcol <- data.frame(gcol=ifelse(grepl(gcol$id, rownames(dta)), gcol$id, "Other"))

		# set param.DE for sample comparison
		if (is.list(param.DE)) {
			if (scale == "none") {
				param.DE$scale <- "none"
			} else if (scale == "upper") {
				param.DE$scale <- "upperquartile"
			} else {
				message("scale", scale, " is unknown")
			}
			param.DE$out <- "top"
			param.DE$species <- species
		}

		# transformation for EDA plots
		param.preprocess <- list(ps=1, logf='log', col.mean.center=T) 

		# Profiles of normalization parameters
		#		param.norms <- list(
		#												list(method="none", scale=scale),
		#												list(method="RUVg", scale=scale, k=1, contrast=param.DE$contrast, spikeId=spike$id),
		#												list(method="RUVg", scale=scale, k=1, contrast=param.DE$contrast, topN=5000, species=species),
		#												list(method="RUVr", scale=scale, k=1, contrast=param.DE$contrast, spikeId=spike$id, rmSpike=T),
		#												list(method="RUVs", scale=scale, k=1, contrast=param.DE$contrast, spikeId=spike$id, rmSpike=T)
		#												)
		param.norms <- list(
												#										list(method="none", scale=scale),
												list(method="RUVg", scale=scale, k=k, spikeId=spike$id),
												#										list(method="RUVg", scale=scale, k=k, topN=5000, species=species),
												#										list(method="RUVr", scale=scale, k=k, spikeId=spike$id, rmSpike=T),
												list(method="RUVs", scale=scale, k=k, spikeId=spike$id, rmSpike=T)
												)

		columns <- 1:ncol(dta)
		column.labels <- colnames(dta[, columns]) 
		group <- sub(rep.sub, '', column.labels, perl=T)
		#param.norm <- param.norms[[1]]
		for (param.norm in param.norms) {
			print("-------------------------------------------------")
			print("-------------------------------------------------")
			if (param.norm$method == "RUVg" && (!any(grep("topN", names(param.norm)))) && 
					is.list(spike) && any(grep("id", names(spike))) ) {
				ids <- grep(spike$id, rownames(dta), value=T)
				if (length(ids) == 0) {
					message("Skipping ", param.norm)	
					next
				}
			}
			# set ofle.tag
			ofle.tag <- paste(tag, "_sc", param.norm$scale, "_RS", param.norm$method, sep="")
			if (param.norm$method == "RUVg") { 
				ofle.tag <- paste(ofle.tag,  
													ifelse(any(grep("spikeId", names(param.norm))), 
																 "_spikeIn", 
																 paste("_topN", param.norm$topN, sep="")),
													sep="")
			}
			title <- ofle.tag

			# normalize
			nrm <- normalize.RUVseq(dta, group, column.labels, param.norm, gcol, title, ofle.tag)

			if (!is.list(nrm))
				next

			plot.new(); pdf(paste(ofle.tag,".pdf",sep=""),h=10,w=10)
			## remove after making sure they are reproduced in EDA plots
			#library(RColorBrewer)
			#colors <- brewer.pal(length(group), "Set2")
			#x <- as.numeric(as.factor(group))
			#plotRLE(nrm[["set"]], outline=F, ylim=c(-4, 4), col=colors[x], main=title)
			#plotPCA(nrm[["set"]], col=colors[x], cex=1.2, main=title)

			#plot

			get.exploratoryAnalysisPlots(nrm[["dta"]], rep.sub=rep.sub, title=title,
																	 param.preprocess=param.preprocess, gcol=nrm[["gcol"]])

			# DE
			if (param.DE$method == "edgeR") {
				prm <- param.DE
				for (i in 1:length(param.DE$contrast)) {
					prm$contrast <- param.DE$contrast[i]
					top <- caller.DE.edgeR(nrm[["set"]], group, param=prm,
																 tag=paste(c(ofle.tag,"_contrast", prm$contrast, "_tag", param.DE$tag[i]), collapse=""))
				}
			} else {
				warnings(param.DE$method, " is not implemented")
			}

			dev.off()
		}
	}

normalize.RUVseq <-
	function(dta, group, column.labels, param, gcol, title, ofle.tag) {
		ggplot <- set.ggplot()

		message("Doing normalize.RUVseq param: ", paste.namedList(param))

		set <- newSeqExpressionSet(as.matrix(dta),
															 phenoData = data.frame(group, 
																											row.names=column.labels))
		if (param$scale != "none") 
			set <- betweenLaneNormalization(set, which=param$scale)

		if (param$method == "RUVg") {
			if (any(grep("spikeId", names(param)))) {
				rows <- grepl(param$spikeId, rownames(dta))
			} else { # emprical
				#param2 <- list(out="topNinvert", topN=param$topN, contrast=param$contrast, species=param$species)
				param2 <- list(out="topNinvert", topN=param$topN, species=param$species)
				if (param$scale == "none") {
					param2$scale <- "none"
				} else if (param$scale == "upper") {
					param2$scale <- "upperquartile"
				} else {
					message("scale", param$scale, " is unknown")
				}
				rows <- caller.DE.edgeR(set, group, param=param2, tag=ofle.tag)
				gcol <- data.frame(group=ifelse (rows, "negCont", "other"))
			}
			ids <- rownames(dta)[rows]
			if (length(ids) == 0) {
				message("No negative control genes found")
				return(-1)
			}
			set <- RUVg(set, ids, k=param$k) 
		} else if (param$method == "RUVs") {
			# https://github.com/drisso/RUVSeq/issues/5
			makeGroups <- function(xs) {
				xs <- factor(xs)
				groups <- matrix(-1, nrow = length(levels(xs)), ncol = max(table(xs)))
				for (i in 1:length(levels(xs))) {
					idxs <- which(xs == levels(xs)[i])
					groups[i,1:length(idxs)] <- idxs
				}
				groups
			}
			differences <- makeGroups(group)
			print(differences)
			ids <- rownames(dta)
			if (any(grep("spikeId", names(param))) && 
					any(grep("rmSpike", names(param))) && 
					param$rmSpike == T) {
				ids <- grep(param$spikeId, rownames(dta), value=T, invert=T)
			}
			set <- RUVs(set, ids, k=param$k, differences)
		} else if (param$method == "RUVr") {
			param2 <- list(out="res")
			if (param$scale == "none") {
				param2$scale <- "none"
			} else if (param$scale == "upper") {
				param2$scale <- "upperquartile"
			} else {
				message("scale", param$scale, " is unknown")
			}
			res <- caller.DE.edgeR(set, group, param=param2)
			ids <- rownames(dta)
			if (any(grep("spikeId", names(param))) && 
					any(grep("rmSpike", names(param))) && 
					param$rmSpike == T) {
				ids <- grep(param$spikeId, rownames(dta), value=T, invert=T)
			}
			set <- RUVr(set, ids, k=param$k, res)
		}
		print(pData(set))

		if (!(param$method == "none" && param$scale == "none"))
			dta <- data.frame(set@assayData$normalizedCounts)

		print(summary(dta))
		writeData(dta, ofle=paste(ofle.tag, "_count.txt.gz", sep=""), 
							commentLine=title, row.names=T)
		list(dta=dta, gcol=gcol, set=set)
	}

annotate <- 
	function(top, species=NA, keytype="REFSEQ") {
		if (species == "Mm") {
			require(org.Mm.eg.db)
			db <- org.Mm.eg.db
		} else if (species == "Hs") {
			require(org.Hs.eg.db)
			db <- org.Hs.eg.db
		} else {
			stop("Species is undefined:", species)
		}
		top$SYMBOL <- select(db,keys=rownames(top),columns="SYMBOL", keytype=keytype)[,2]
		top$ENTREZID <- select(db,keys=rownames(top),columns="ENTREZID", keytype=keytype)[,2]

		top
	}



