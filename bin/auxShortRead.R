library(ShortRead)
library(parallel)

source(paste0(bindir, "/bigData.R"));

test.caller.qa <-
	function() {
		source(paste0(bindir, "auxShortRead.R"))

		# 1.
		dataDir <- "~/Projects/biter_wosczyna_miRNAsInMuscleCellFate/data/Wosczyna_2015/PullDown/rawData/"
		meta <- list(fle="~/Projects/biter_wosczyna_miRNAsInMuscleCellFate/data/Wosczyna_2015/PullDown/rawData/metadata",
								 sour="Alias", dest="Name", pair="Pair")
		fles <- caller.alias2name(dataDir=dataDir, meta=meta, pattern=".fastq.gz", nc=1)
		qa <- caller.qa(fles, outDir="qualityControl", n=1000, save.image=F, nc=1)
		#caller.qa(dataDir=dataDir, meta=meta, outDir="qualityControl", n=10000000, save.image=F, nc=2)

		# 2.
		dataDir <- "~/Projects/biter_cvanvelt_SCTrxTrl/data/Gonzalez_2014_JNeurosci/PullDownTrl/rawData/"
		meta <- list(fle="~/Projects/biter_cvanvelt_SCTrxTrl/data/Gonzalez_2014_JNeurosci/PullDownTrl/rawData/metadata",
								 sour="Alias", dest="Name", pair="Pair")
		fles <- caller.alias2name(dataDir=dataDir, meta=meta, pattern=".fastq.gz", nc=1)
		qa <- caller.qa(fles, outDir="qualityControl", n=1000, save.image=F, nc=1)

		# 3.
		dataDir <- "~/Projects/biter_demorree_PAS/data/Soleimani_2012_DevCell/ChIPseq/rawData/"
		meta <- list(fle="~/Projects/biter_demorree_PAS/data/Soleimani_2012_DevCell/ChIPseq/rawData/metadata",
								 sour="Alias", dest="Name", pair="Pair")
		fles <- caller.alias2name(dataDir=dataDir, meta=meta, pattern=".fastq.gz", nc=1)
		qa <- caller.qa(fles, outDir="qualityControl", n=1000, save.image=F, nc=1)

		

		# TODO test
		dataDir <- "~/Projects/biter_demorree_PAS/data/deMorree_2012/directRNAseq/rawData"
		meta <- list(fle="~/Projects/biter_demorree_PAS/data/deMorree_2012/directRNAseq/rawData/metadata",
								 sour="Alias", dest="Name", pair="Pair")
		#meta <- NULL
		fles <- caller.alias2name(dataDir=dataDir, meta=meta, pattern=".fasta.gz", nc=1, convert2fastq=T)
		qa <- caller.qa(fles, outDir="qualityControl", n=1000, save.image=F, nc=1)
		
		# ShortRead only works for single end ungapped alignments 
#		dataDir <- "~/Projects/biter_cvanvelt_SCTrxTrl/data/vanVelthoven_2015/PullDownTrl/RM_TR_ada/"
#		meta <- list(fle="~/Projects/biter_cvanvelt_SCTrxTrl/data/vanVelthoven_2015/PullDownTrl/RM_TR_ada/metadata",
#								 sour="Name", dest="Name", pair="Pair")
#		caller.qa(dataDir=dataDir, meta=meta, pattern=".sorted.bam", outDir="qualityControl", n=1000, save.image=F, nc=2, type="BAM")

		caller.qa(dir(pattern='.fastq.gz'), outDir='/home/biter/Projects/biter_lingliu_DNAdamage/results/2017-03-13/QA_Liu_2017_ATACseq', nc=1)

	}

# As a side effect converts fasta files to fastq
caller.alias2name <- 
	function(dataDir=".", meta=NULL, pattern=".fastq.gz", nc=1, convert2fastq=F) {

		if (is.null(meta)) {
			sfles <- dir(dataDir, pattern, full.names=T)
			fles  <- basename(sfles)
		} else if (file.exists(meta$fle)) {
			dtam  <- readBigTable(meta$fle, header=T)[,c(meta$sour, meta$dest, meta$pair)]
			sfles <- paste(dataDir, "/", dtam[,meta$sour], pattern, sep="");
			if (any(grep(".fast", pattern, fixed=T)) && !any(is.na(dtam[,meta$pair]))) {
				fles  <- paste(dtam[,meta$dest], "_", dtam[,meta$pair], pattern, sep="")
			} else {
				ndi <- !duplicated(sfles)
				sfles <- sfles[ndi]
				fles  <- paste(dtam[ndi,meta$dest], pattern, sep="")
			}
		} else {
			stop("ERROR: ", meta$fle, " does not exist")
		}

		if ( convert2fastq && (any(grep(".fasta", pattern, fixed=T))) ) {
			fles <- gsub(".fasta.", ".fastq.", fles)
			mclapply(1:length(fles), 
							 function(i, s=sfles, d=fles) { 
								r <- readDNAStringSet(s[i]); 
								writeXStringSet(r, d[i], format="fastq", compress=T) 
							 }, mc.cores=nc)
		} else {
			file.symlink(sfles, fles)
		}

		print(sfles)
		print(fles)
		fles
}

#Sequencing is done outside R; at the start of the work flow we have access to fastq files, in the fastqDir directory. A first step after sequencing might use ShortRead to produce a quality assessment report. Down-sample fastq files if they are big.
#TODO use nc parameter
#Takes ~3h for Ikeda_2014
caller.qa <- 
	function(fles, outDir=".", save.image=F, nc=1, pattern='.fastq.gz', ...) {
		qa <- qa(fles, BPPARAM=MulticoreParam(workers=nc), ...)
		#qa <- qa(fles, BPPARAM=MulticoreParam(workers=nc))
		url <- report(qa, dest=outDir)

		append <- F
		ofle <- file.path(outDir, "qa.txt")
		for (n in names(qa)) {
			print(n)
			df <- qa[[n]]
			if (class(df) == "data.frame") {
				dtal <- list(x=df)
				names(dtal) <- n
			} else {
				dtal <- df
				names(dtal) <- paste(n, names(dtal), sep=".")
			}
			lapply(1:length(dtal), 
						 function(i) { 
								df <- dtal[[i]]
								rownames(df) <- gsub(pattern, "", rownames(df))
								writeData(df, ofle, row.names=T, col.names=T, commentLine=names(dtal[i]), append=append)
							 })
			append <- T
		}

		library(tidyr)
		library(dplyr)
		df <- qa[["baseCalls"]]
		rownames(df) <- gsub(pattern,"", rownames(df))
		dfs <- data.frame(GC_o_ACTG = df %>% dplyr::select(one_of("C","G")) %>% rowSums / df %>% dplyr::select(one_of("A","C","G","T")) %>% rowSums)
		dfs
		writeData(dfs, file.path(outDir, "baseCalls.txt"), row.names=T, col.names=T)

		if (save.image)
			save(qa, file=file.path(outDir, "qa.rda"))

		print(sessionInfo())
		qa
	}
