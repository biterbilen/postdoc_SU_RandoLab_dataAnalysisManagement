library(TFBSTools)

source("~/Projects/biter_biter_shared/bin/auxGenomicFeatures.R")

is.null <- 
	function(x) {
		base::is.null(x) || x=="NULL"
	}


getSeqsFromGRanges <- 
	function(dta.gr, gnm="mm10", wnames=F, slop=slop) {
		# Extend the segment length
		if (! is.null(slop))
			dta.gr <- flank(dta.gr, slop, both=T)
		#TODO fix ends if necessary

		genSeq <- getBSgenome(gnm)
		DNAStringSet <- getSeq(genSeq, dta.gr)
		if (wnames)
			names(DNAStringSet) <- values(dta.gr)$name
		DNAStringSet	
	}

test.getSiteCountswConsensusStringAndDNAStringSet <-
	function() {
		source("./auxTFBSTools.R")
		load('getSeqFeatures.RData')
		cnts.mirs <- getSiteCountswConsensusStringAndDNAStringSet(consensusStringList, DNAStringSet, rev.comp=F, mc.cores=1)
		head(cnts.mirs[1:5,1:5])

	}

test.smithWatermanAlignment <-
	function() {
		source("./auxTFBSTools.R")
		#mir <- reverseComplement(RNAString("AGAGGUAGUAGGUGCUUUAGAAGAU"))
		mir <- RNAString("AGAGGUAGUAGGUGCUUUAGAAGAU")
		mir <- RNAString("ANNNNNNNUAGGUGCUUUAGAAGAU")
		mRNA <- RNAString("AAAAAUGAGGUACUUCGUAGUAUACAAUUUCGU")
		mRNA1 <- RNAString("UUAAAAAUGAGGUACUUCGUAGUAUACAAUUUCGUAA")
		mRNA2 <- RNAString("UUAAAAAUGAGGUACUUCCUAGUAUACAAUUUCGUAA")
		mRNA3 <- RNAString("UUAAAAAUGACGUACUUCCUAGUAUACAAUUUCGUAA")
		smithWatermanAlignment(mir, mRNA, gapOpening=1)
		smithWatermanAlignment(mir, mRNA, gapOpening=2)
		smithWatermanAlignment(mir, mRNA, gapOpening=5)
		smithWatermanAlignment(mir, mRNA, gapOpening=0, gapExtension=1)
		smithWatermanAlignment(mir, mRNA, gapOpening=1, gapExtension=1) # mimicks the Khorshid paper miRNA binding pattern pattern
		# mimicks the Khorshid paper miRNA binding pattern pattern
		for (m in list(mRNA, mRNA1, mRNA2, mRNA3)) {
			print(c("mir", as.character(mir)))
			print(c("mRNA", as.character(mRNA)))
			print(smithWatermanAlignment(mir, m, gapOpening=3, gapExtension=0))
			print("\n")
		}

	}

# default params mimick Khorshid paper miRNA-mRNA binding pattern
smithWatermanAlignment <- 
	function(s1, s2, mat=NULL, type="local", match=1, mismatch=-1, gapOpening=1, gapExtension=1 ) {
		# First use a fixed substitution matrix
		if (!is.null(mat))
			mat <- nucleotideSubstitutionMatrix(match = match, mismatch = mismatch, baseOnly = TRUE)
		pairwiseAlignment(s1, s2, type = type, substitutionMatrix = mat,
																		gapOpening = gapOpening, gapExtension = gapExtension)

	}

getSiteCountswConsensusStringAndDNAStringSet <- 
	function(consensusStringList, DNAStringSet, rev.comp=T,  mc.cores=1) 
	{
		cnts <- mclapply(consensusStringList, 
										 function(x, ds=DNAStringSet) {
											 x <- as.character(x)
											 cnts <- vcountPattern(x, ds)
											 if (rev.comp) { 
												 x.rc <- as.character(reverseComplement(DNAString(x)))
												 cnts <- cnts + vcountPattern(x.rc, ds)
											 }
											 cnts
										 }, mc.cores=mc.cores)
		head(cnts) %>% print
		cnts <- as.data.frame(cnts)
		colnames(cnts) <- names(consensusStringList)
		rownames(cnts) <- names(DNAStringSet)
		#cnts$name <- names(DNAStringSet)
		cnts
	}

getSiteCountswMatrixAndDNAString <- 
	function(MatrixList, DNAString, min.score="95%", rev.comp=T,  mc.cores=1) {
		res <- mclapply(1:length(MatrixList),
										function(j, ds=DNAString) { 
											pwm <- Matrix(MatrixList[[j]])
											cnts <- countPWM(pwm, ds, min.score=min.score)
											if (rev.comp) { 
												pwm.rc <- reverseComplement(pwm)
												cnts <- cnts + countPWM(pwm.rc, ds, min.score=min.score)
											}
											cnts
										}, mc.cores=mc.cores)
		cnts <- do.call(c,unname(res))

	}

convertPFM2PWMatrixList <- 
	function(MatrixList) {
		# convert to PWM since searchSeq accepts only PWM
		if (class(MatrixList) == "PFMatrixList") { 
			MatrixList <- lapply(MatrixList, toPWM)
			class(MatrixList) <- "PWMatrixList"
		}
		MatrixList
	}

# Takes 1 hour for 17K regions of 2.2K nc
getSiteCountswMatrixAndDNAStringSet <- 
	function(MatrixList, DNAStringSet, min.score="95%", rev.comp=T, add.seqname=T, mc.cores=1, ...) {

		MatrixList <- convertPFM2PWMatrixList(MatrixList)

		cnts <- lapply(1:length(DNAStringSet),
									 function(i) { 
										 if (i %% 500 == 0) print(i)
										 getSiteCountswMatrixAndDNAString(MatrixList, DNAStringSet[[i]], min.score=min.score, rev.comp=rev.comp, mc.cores=mc.cores);
									 })

		cnts <- data.frame(t(as.data.frame(cnts)))
		colnames(cnts) <- paste(gsub(" ", "_", name(MatrixList)), names(MatrixList), sep="|")
		if (add.seqname)
			#cnts$name <- names(DNAStringSet)
			row.names(cnts) <- names(DNAStringSet)

		cnts
	}
# TODO make more efficient
getSiteCountswMatrixAndBedCoords <- 
	function(MatrixList, DNAStringSet, mc.cores=1) {
		res <- mclapply(1:length(DNAStringSet),
										function(i) { x=as.character(DNAStringSet[i]); 
										searchSeq(MatrixList, x, seqname=names(x), min.score=min.score, strand=strand) 
										}, mc.cores=mc.cores)
		sitesetListAll <- do.call(c, unname(res))
	}

getSiteswMatrixAndBedCoords <- 
	function(MatrixList, DNAStringSet, mc.cores=1, 
					 outfileName=NULL, min.score="99%", strand="*") {

		# convert to PWM since searchSeq accepts only PWM
		MatrixList <- convertPFM2PWMatrixList(MatrixList)

		library(Biostrings)
		res <- mclapply(1:length(DNAStringSet),
										function(i) { x=as.character(DNAStringSet[i]); 
										searchSeq(MatrixList, x, seqname=names(x), min.score=min.score, strand=strand) 
										}, mc.cores=mc.cores)
		sitesetListAll <- do.call(c, unname(res))

		sites.gff3 <- writeGFF3(sitesetListAll)

		if (!is.null(outfileName))
			writeData(sites.gff3, outfileName, row.names=F, col.names=F)
		sites.gff3
	}

test.readMatrixSet <-
	function() {
		source("~/Projects/biter_biter_shared/bin/auxTFBSTools.R")
		# PFM matrix
		fle <- "~/PI_HOME/Data/casco/JASPAR/JASPAR_CORE_nonredundant_pfm_vertebrates.txt"
		mspfm <- readMatrixSet(fle)
		toPWM(mspfm[[1]])

		fle <- "~/PI_HOME/Data/casco/Homer/custom.motifs"
		mspwm <- readMatrixSet(fle, opts=list(matrixtype="PWM", ii=2, byrow=T))

		src <- "~/PI_HOME/Data/casco/CIS-DB/Homo_sapiens/all_pwms.txt"
		k3 <- getMatrixList(src=src, gnm="mm10", matrixtype="PWM", ii=1, ni=1, byrow=F)
		k3[[1]]

	}


readMatrixSet <- 
	function(fle, opts=NULL) {

		opts <- mergeSetwDefault(opts, list(matrixtype="PFM", ni=2, ii=1, byrow=F, split=" |\t"))
		f <- get(paste(opts$matrixtype,"atrix", sep=""))
		ff <- get(paste(opts$matrixtype,"atrixList", sep=""))
		MatrixSet <- ff()
		M <- NULL
		n <- NULL
		id <- NULL
		con <- file(fle, open="r")
		while( length(ln <- readLines(con, 1, warn=F))>0 ) {
			ch <- strsplit(ln, "")[[1]][1]
			if (ch == ">") {
				# current
				if (!is.null(M)) {
					# correct CIS-DB matrix entries
					M[M> 1] = 1.0 
					# in case there are negative values
					M[M< 0] = 0.0 
					m <- f(ID=id, name=n, strand="+",
								 bg=c(A=0.25, C=0.25, G=0.25, T=0.25),
								 #tags=list(...),
								 profileMatrix=matrix(M, nrow=4, byrow=opts$byrow, 
																		dimnames=list(c("A", "C", "G", "T"))))
					MatrixSet <- c(MatrixSet, list(m))
				}
				# next
				M <- NULL
				n <- gsub("^>", "", strsplit(ln, opts$split, perl=T)[[1]][opts$ni])
				id <- gsub("^>","", strsplit(ln, opts$split, perl=T)[[1]][opts$ii])
			} else {
				if (opts$matrixtype == "PWM") {
					M <- rbind(M, as.numeric(unlist(strsplit(ln, opts$split, perl=T)[[1]])))
				} else {
					# integer casting is necessary for TFBSTools PWM conversion
					M <- rbind(M, as.integer(unlist(strsplit(ln, opts$split, perl=T)[[1]])))
				}

			}
		}
		close(con)
		# last
		m <- f(ID=id, name=n, strand="+",
					 bg=c(A=0.25, C=0.25, G=0.25, T=0.25),
					 #tags=list(...),
					 profileMatrix=matrix(M, nrow=4, byrow=opts$byrow, 
																dimnames=list(c("A", "C", "G", "T"))))
		MatrixSet <- c(MatrixSet, list(m))
		names(MatrixSet) = ID(MatrixSet)
		print(MatrixSet)
		MatrixSet
	}

test.getMatrixList <-
	function() {
		source("~/Projects/biter_biter_shared/bin/auxTFBSTools.R")

		# test web file read
		src <- "http://homer.salk.edu/homer/custom.motifs"
		src <- "~/PI_HOME/Data/casco/Homer/custom.motifs"
		k2 <- getMatrixList(src=src, gnm="mm10", matrixtype="PWM", ii=2, byrow=T)


		# test all_versions
		source("~/Projects/biter_biter_shared/bin/auxTFBSTools.R")
		k1 <- getMatrixList(src="JASPAR2014", gnm="mm10", matrixtype="PWM", all_versions=T)
		k1 <- getMatrixList(src="JASPAR2014", gnm="mm10", matrixtype="PFM", all_versions=F)
		k1[[1]]

		source("~/Projects/biter_biter_shared/bin/auxTFBSTools.R")
		k1 <- getMatrixList(src="JASPAR2014", gnm="mm10", matrixtype="PWM", all_versions=T)
		k1[[1]]

		src <- "~/PI_HOME/Data/casco/Homer/custom.motifs"
		k2 <- getMatrixList(src=src, gnm="mm10", matrixtype="PWM", ii=2, ni=2, byrow=T)
		k2[[1]]

		src <- "~/PI_HOME/Data/casco/CIS-DB/Mus_musculus/all_pwms.txt"
		select <- "~/Projects/biter_biter_stemCellFateRegulation/results/2016-06-19/Motif_ID.TF_Names"
		ofle <- "CIS.selected.pwm"
		k3 <- getMatrixList(src=src, gnm="mm10", matrixtype="PWM", ii=1, ni=1, byrow=T, select=select, ofle=ofle)
		isTRUE(class(k3)[1]) 

		src <- "~/PI_HOME/Data/casco/CIS-DB/Homo_sapiens/all_pwms.txt"
		select <- "~/Projects/biter_biter_stemCellFateRegulation/results/2016-06-28/ivTFBS_test_hg19_Synapse_syn6131484/PWMPrep/TF.selected"
		#ofle <- "CIS.selected.pwm"
		#k4 <- getMatrixList(src=src, gnm="hg19", matrixtype="PWM", ii=1, ni=1, byrow=T)
		k4 <- getMatrixList(src=src, gnm="hg19", matrixtype="PWM", ii=1, ni=1, byrow=T, select=select)
		k4 <- getMatrixList(src=src, gnm="hg19", matrixtype="PWM", ii=1, ni=1, byrow=T, select=NA)
		isTRUE(class(k4)[1]) 
	}

# select is a file w/ columns Matrix_ID and TF_Name 
getMatrixList <- 
	function(src="JASPAR2014", gnm="mm10", select=NA, ofle=NA, ...) {

		if (gnm == "mm10") {
			spc <- 10090
		}
		else if (gnm == "hg19") {
			spc <- 9606
		} 
		else if (gnm == "hg38") {
			spc <- 9606
		} 
		else {
			stop("Set species for ", gnm)
		}

		opts <- mergeSetwDefault(list(...), list(species=spc, gnm=gnm))
		print(opts)
		if (src == "JASPAR2014") {
			library(JASPAR2014)
			MatrixList <- getMatrixSet(JASPAR2014, opts)
		} else {
			MatrixList <- readMatrixSet(src, opts)
		}

		if (!is.na(select)) {
			sdta <- readBigTable(select)
			MatrixList <- MatrixList[as.character(sdta$V1)]
			for (i in 1:length(MatrixList)) {
				MatrixList[[i]]@name <- as.character(sdta[sdta$V1 == ID(MatrixList[[i]]),2]); 
			}
		}

		if (!is.na(ofle))
			writeData(MatrixList, ofle)

		MatrixList
	}

countSegmentswMatrix <- 
	function(sites.gff3, nms) {
		unlist(lapply(nms, 
									function(x) {
										pat <- paste("TF=", x, ";", sep="");
										length(unique(sites.gff3[grep(pat, sites.gff3$attributes),"seqname"])) 
									})) 
	}

getEnrichment <- 
	function(cs, fun.type="fisher.test", ...) {
		fun <- get(fun.type)

		enrichment <- as.data.frame(t(apply(cs, 1, 
																				function(x) {
																					m <- matrix(x, nrow=2, byrow=T)
																					tst <- fun(m, ...)
																					c(tst$p.value, tst$estimate)
																				})))
		names(enrichment) <- paste(fun.type, c(".p.value", ".estimate"), sep="")
		enrichment	

	}

