library(GenomicRanges)
library(dplyr)
library(tidyr)
library(magrittr)

#TODO check package groups to make this package object oriented

#TODO implement 
test.joinFiles <- 
	function() {

		#flelist <- c("refGene.gz", "filtered.peaks.narrowPeak"),
		source("~/Projects/biter_biter_shared/bin/bigData.R")

		regions <- readBigTable("promoters.bed.gz", col.names=c("seqnames","start","end","name","score","strand"))
		pmeth <- readBigTable("./pmeth_promoters.tab.gz",header=T)
		TFBS <- readBigTable("./promotersJASPAR2014.txt.gz",header=T); names(TFBS)[1] <- "name"
		adta <- merge(merge(regions[,1:4], pmeth), TFBS) 
		
		iA <- grep("pmeth.A4", names(adta))
		iQ <- grep("pmeth.Q4", names(adta))
		nTFBS <- grep(".MA", names(adta), fixed=T,value=T)
		
		library(ggplot2)
		library(reshape2)
		ps <- 0.0001
		#plot.new(); pdf("LM.pdf",h=8,w=8)
		summary(adta[c(iA,iQ)])
		dta.pmeth <- data.frame(pmeth.A=rowMeans(adta[iA]), pmeth.Q=rowMeans(adta[iQ]))
		dta.pmeth$pmeth.AoQ <- (dta.pmeth[,1]+ps)/(dta.pmeth[,2]+ps)
		dta.pmeth$log2.pmeth.AoQ <- log2(dta.pmeth$pmeth.AoQ)
		summary(dta.pmeth)
		ggplot(melt(dta.pmeth[,1:2]), aes(x=value, color=variable)) + geom_density()
		ggplot(melt(dta.pmeth[,3]), aes(x=value)) + geom_density() + labs(x="pmeth.AoQ")
		ggplot(melt(dta.pmeth[,4]), aes(x=value)) + geom_density() + labs(x="log2.pmeth.AoQ")

		dta <- data.frame(log2.pmeth.AoQ=dta.pmeth$log2.pmeth.AoQ, adta[nTFBS])
		head(dta)
		frml <- paste(names(dta)[1], paste0(names(dta)[-1], collapse="+"), sep="~")
		mdl <- lm(as.formula(frml), data=dta)
		summary(mdl)

		dev.off()
	}

#TODO implement 
joinFiles <- 
	function(flelist, indlist) {
	nFiles <- length(flelist)
	if (nFiles < 2) 
		message("ERROR: Supply >1 files")

	for (i in 2:nFiles) {
		print(i)			
	}	
}

paste.namedList <- function(lst, sp=":", sp2=",") {
	  paste(lapply(names(lst),
					function(x) paste(x, sp, lst[x], sep="")), collapse=sp2)
}

test.getFileFormat <- function(filename) {

	source("~/Projects/biter_biter_shared/bin/bigData.R")
	getFileFormat("a.sorted.bam")
	getFileFormat("a.sorted.bam.txt")

}

getFileFormat <- function(filename) {
	if (grepl(".bam$", filename)) {
		return("bam")
	} else {
		return("text")
	}

}

#TODO add after format specifications
# is not robust
getDataFormat <- function(dta) {
	if (class(dta)[1] == "GRanges")
		return("GRanges")
	if (class(dta)[1] == "PWMatrixList")
		return("PWMatrixList")
	if (names(dta)[3]=="T" && names(dta)[4]=="M") 
		return("msr");
	if (names(dta)[1]=="chr" && names(dta)[2]=="start" && names(dta)[3]=="end") 
		return("bedplus");
	if (names(dta)[1]=="chr" && names(dta)[2]=="source" && 
			names(dta)[3]=="feature" && names(dta)[9]=="attributes")
		return("gtf");
	if (any(grep("n$", names(dta), perl=T)) && any(grep("y$", names(dta), perl=T))) 
		return("cov") 
	if (class(dta)[1] == "data.frame")
		return("data.frame")
	 
	return("unknown")
}

convertDataFormat <- function(dta, frmt="msr", nms=NULL, desc.tag=NA, posStrand=F) {

	if (!is.null(nms))
		names(dta) <- nms

	cdta <- NULL
	if (getDataFormat(dta) == "cov" && frmt == "msr") {
		cdta <- dta[,1]
		ctda <- dta[,2]+1 # cov files are 0-based
		cdta['T'] <- dta[4] + dta[5]
		cdta['M'] <- dta[4];
		names(cdta) <- c("chr","start","T","M")
	}
	if (getDataFormat(dta) == "bedplus" && frmt == "GRanges") {
		if (any(grep("name", names(dta)))) {
			name <- dta$name
		} else {
			name <- paste(dta$chr, ":", dta$start, "-", dta$end, sep="")
		}
		if (any(grep("score", names(dta)))) {
			score <- dta$score
		} else {
			score  <- rep(1000,length(name))
		}
		if (any(grep("strand", names(dta)))) {
			strand <- gsub("\\.", "*", dta$strand)
		} else {
			strand <- rep("*", length(name))
		}
		if (posStrand)
			strand <- gsub("\\*", "+", strand)
		cdta <- GRanges(seqnames=dta$chr, 
										ranges=IRanges(dta$start+1, dta$end),
										name=name, score=score, strand=strand)
		print(head(cdta))
	}
	if (getDataFormat(dta) == "GRanges" && frmt == "bedplus") {
		if (any(grep("name",names(mcols(dta))))) {
			i <- grep("name",names(mcols(dta)))
			name <- as.character(unlist(mcols(dta)[i]))
			mcols(dta)[[i]] <- NULL
		} else {
			name   <- 1:length(dta)
		}
		if (any(grep("score",names(mcols(dta))))) {
			i <- grep("score",names(mcols(dta)))
			score <- as.character(unlist(mcols(dta)[i]))
			mcols(dta)[[i]] <- NULL
		} else {
			score  <- rep(1000,length(dta))
		}
		strand <- gsub("\\*", ".", as.character(strand(dta)))
    cdta <- data.frame(chr=as.character(seqnames(dta)), start=start(dta)-1, #bed files are 0-based
											 end=end(dta), name=name, score=score, 
											 strand=strand, values(dta))
	}
	if (getDataFormat(dta) == "GRanges" && frmt == "gtf") {
		# TODO expand fields
		cdta <- GRanges(seqnames=dta$chr, ranges=IRanges(start=dta$start,end=dta$end),strand=dta$strand)
	}
	if (getDataFormat(dta) == "GRanges" && frmt == "msr") {
    cdta <- data.frame(chr=as.character(seqnames(dta)), start=start(dta), #bed files are 0-based
											 values(dta))
		print(names(cdta))
	}
	if (getDataFormat(dta) == "data.frame" && frmt == "gmt") {
		cdta <- NULL
		for (i in 1:ncol(dta)) {
			name <- paste0(desc.tag, ":", names(dta)[i])
			desc <- name
			ids <- dta[dta[,i]>0,] %>% rownames
			if (length(ids) > 0)
				cdta %<>%	bind_rows(data.frame(name, desc, paste(ids, collapse="\t")))
		}
	}
	
	return(cdta)
}

test.writeData <- 
	function () {
		source("~/Projects/biter_biter_shared/bin/bigData.R")

		# append
		dta <- data.frame(a=1:10, b=2:11); append <- F
		writeData(dta, "a.txt", row.names=T, append=append, col.names=T, commentLine="f1")
		writeData(dta[,1], "a.txt", row.names=T, append=append, col.names=F, commentLine="f1")
		dta <- data.frame(a=0:3, b=2:5, c=3:6); append <- T
		rownames(dta) <- c("a","b","c","d")
		writeData(dta, "a.txt", append=T, row.names=T, col.names=T, commentLine="f2")
		writeData(dta, "a.gmt", frmt="gmt", row.names=F, col.names=F, desc.tag="dummy")

	}

writeData <- function(dta, ofle, frmt=NULL, desc.tag=NA, commentLine=NULL, verbose=T, ... ) {

	if ((!is.null(frmt)) && getDataFormat(dta) != frmt) {
		dta <- convertDataFormat(dta=dta, frmt=frmt, desc.tag=desc.tag)
	}

	params <- list(...)
	if (isTRUE(frmt=="GRanges")) { # TODO add here others too
		saveRDS(dta, ofle)
		len <- length(dta)
	} else {
		fmode <- 'w'
		if (any(grepl("append", names(params))) && params$append == T) 
			fmode <- 'a'
		if (any(grep(".gz$", ofle, perl=T))) {
			con <- gzfile(ofle, fmode, compression=9);
		} else {
			con <- file(ofle, fmode)
		}
		on.exit(close(con))

		if (isTRUE(class(dta)[1] == "PWMatrixList")) {
			len <- length(dta)
			writeLines("MEME version 4\nALPHABET= ACGT\nstrands: + -\nBackground letter frequencies\nA 0.25 C 0.25 G 0.25 T 0.25", con=con)
			for (i in 1:length(dta)) {
				writeLines(paste("MOTIF ", ID(dta[[i]]), " ", name(dta[[i]]), sep=""), con=con)
				writeLines(paste("letter-probability matrix: alength=4 w=", ncol(dta[[i]]@profileMatrix), sep=""), con=con)
				#writeLines(paste(">", ID(dta[[i]]), " ", name(dta[[i]]), sep=""), con=con)  
				# this is ACTG as rows
				#write.table(dta[[i]]@profileMatrix, file=con, append=T, row.names=F, col.names=F, sep="\t") 
				# ACTG as column
				write.table(t(dta[[i]]@profileMatrix), file=con, append=T, row.names=F, col.names=F, sep="\t") 
			}
		} else {
			if (! is.null(commentLine))
				writeLines(paste("#", commentLine), con=con)
			# rownames adjustment
			if (any(grepl("row.names", names(params))) && params$row.names == T)  {
				dta <- cbind(rownames(dta), dta)
				names(dta)[1] <- "Id"
				params$row.names <- F
			}
			do.call("write.table", c(list(x=dta, file=con, sep="\t", quote=F), params)) 

			len <- dim(dta)[1]
		}
	}
	if (verbose == T) {
		print(paste(ofle, " entries= ", len, sep='')) 
	}
}


#rowSumGroupData <- function(dta, colNames, verbose=F) {
#	return(rowSums(dta[,colNames]))
#}

dummy.sumGroupData <- function(dta, mdta, ofleTag="group", pattern=c("y","n"),
															 format="cov", verbose=F) {
	grps <- as.character(unique(mdta$GroupName))
	for (g in grps) {
		b <- NULL
		smps  <- as.character(mdta[mdta$GroupName==g, "LibraryName"])

		#cov specific format
		if (format == "cov") {
			b <- dta[,1:3]
		}

		#sum group data
		sums <- 1:dim(dta)[2]*0;
		for (p in pattern) {
			smpNames <- as.character(lapply(smps, function(x, m=p) { paste(x, m, sep="") }))
			b[paste(g, p, sep="")] <- rowSums(dta[,smpNames])       
			sums <- sums + b[paste(g, p, sep="")]
		}
		#b[paste(g, "n", sep="")] <- rowSums(a[,smpns])       

		#filter 0 coverage regions
		b <- b[sums>0, ]

		#dimension
		if (verbose) {
			print(paste("group=", g, " libraries=", paste(smps,collapse=","), " size=", dim(b)[1], sep='')) 
		}

		#write to file
		ofle = paste(ofleTag, "_", g, ".gz", sep="");
		gz <- gzfile(ofle,'w',compression=9);
		on.exit(close(gz))
		write.table(b, file=gz, sep="\t", quote=F, row.names=F) 
	}
}



getFreq <- function(dta, breaks=seq(0.0,1,0.1), percent=T, include.lowest=F, mc.cores=1) {
	if (typeof(dta[1,1]) == "double") {
		if (include.lowest) {
			dta[ dta == breaks[1] ] <- dta[ dta == breaks[1] ] + ( (breaks[1] + breaks[2]) / 2)
		}
		dta <- data.frame(apply(dta, 2, function(x) cut(x, breaks=breaks)))
	}
	Reduce(rbind,
				 mclapply(1:dim(dta)[2],
									function(i) {
										dta.f <- melt(table(dta[,i]))
										if (percent)
											dta.f[,2] <- dta.f[,2] / sum(dta.f[,2])
										names(dta.f) <- c("annotation", "count")
										dta.f$sample <- names(dta)[i]
										dta.f
									}, mc.cores=mc.cores))
}
test.readBigTable <- 
	function() {
		source("~/Projects/biter_biter_shared/bin/bigData.R")
		fle <- "~/PI_HOME/Data/casco/ENCODE/hg19/Synapse_syn6131484/labels/MYC.train.labels.tsv.gz"
		readBigTable(fle, nrows=10, verbose=F,header=T) 
		#readBigTable(fle, nrows=10, verbose=T,header=T, col.names=NULL) ; does not run
		readBigTable(fle, nrows=10, verbose=F, header=T,col.names=c("chr","start","end"))
		readBigTable(fle, nrows=10, verbose=F, skip=1, col.names=c("chr","start","end"), frmt="GRanges")
		readBigTable(fle, nrows=10, verbose=F, skip=1, col.names=c("chr","start","end","d1","d2","d3","d4","d5"))
		fle <- "~/PI_HOME/Data/casco/ENCODE/hg19/Synapse_syn6131484/ChIPseq/ChIPseq.HeLa-S3.CEBPB.relaxed.narrowPeak.gz"
		readBigTable(fle, nrows=10, verbose=T, col.names=c("chr","start","end")) 

		readBigTable(fle, nrows=10, verbose=T, col.names=c("chr","start","end"), frmt="GRanges") 

		fle <- "~/Projects/biter_biter_shared/results/2016-04-22/benchmark/a.bam"
		dta.p <- readBigTable(fle, paired=T)
		dta.s <- readBigTable(fle, paired=F)
		source("~/Projects/biter_biter_shared/bin/bigData.R")
		# TODO does not work
		readBigTable(fle, strandMode=1)
		readBigTable(fle, paired=T, strandMode=1)

}


#7GB
readBigTable <- function(filename, filterStr=NA, frmt=NULL, posStrand=F, nrows=NULL, verbose=F, ...) {
	params <- list(...)
	col.names <- params$col.names
	params$col.names <- NULL

	if (getFileFormat(filename) == "bam") {
		library(Rsamtools)
		library(GenomicAlignments)
		flags = ScanBamParam(flag=scanBamFlag(isUnmappedQuery=F,isDuplicate=F,isNotPassingQualityControls=F),what=c('mapq'))
		#print(flags)
		if (params$paired) {
			readBam <- readGAlignmentPairs 
		} else {
			readBam <- readGAlignments 
		}
		dta <- readBam(filename,param=flags)
		# TODO add strandMode argument now works with strandMode=1 default setting
		#dta = do.call(readBam, c(list=c(file=filename,param=flags), params))

	} else if (grepl(".fa", gsub(".gz","", filename), fixed=T) || 
						 grepl(".fasta", gsub(".gz","", filename), fixed=T) ) {
		library(Biostrings)
		#filename <- "/home/biter/PI_HOME/Data/casco/refGeneSeq/mm10/mm10_refGene_canonical.v2.fiveUTRsByTranscript.fa.gz"
		dta <- Biostrings::readDNAStringSet(filename, "fasta") 
		
	} else {

		#tab5rows <- read.table(gzfile(filename), nrows = 100000, ...)
		tab5rows <- do.call("read.table", 
											  c(list(file=gzfile(filename), nrows = 100000), params))
		classes <- sapply(tab5rows, class)

		# get the number of rows
		if (is.null(nrows)) 
			nrows <- as.numeric(system(paste("zless", filename, "| wc -l"), intern = T))

		#dta <- read.table(gzfile(filename), colClasses = classes, nrows = nrows, ...)
		dta <- do.call("read.table", 
									 c(list(file=gzfile(filename), colClasses = classes, nrows = nrows), params))
	}

	if (verbose) {
		print(classes)
		print(dim(dta))
	}		

	if (!is.null(col.names)) {
		n <- min(ncol(dta) , length(col.names))
		names(dta)[1:n] <- col.names[1:n]
	}
	if (!is.na(filterStr)) 
		dta %<>% filter_(filterStr)

	if (verbose)
		print(head(dta))

	if ((!is.null(frmt)) && getDataFormat(dta) != frmt) # if format change required 
		dta <- convertDataFormat(dta, frmt=frmt, posStrand=posStrand);

	return(dta)
}

#14GB
getCoverage <- function(dta) {
	yi <- grep("y$", names(dta), perl = T) 
	ni <- grep("n$", names(dta), perl = T) 
	coverage <- dta[yi] + dta[ni] 
	names(coverage) <- sub("y$", "", names(coverage), perl = T) 

	return(coverage)
}

getFilteredRatio <- function(dta, ps = 1, aveDenom = 0, type = "mean", 
														 verbose = F) {
	#select columns to get the ratio
	yi <- grep("y$", names(dta), perl = T)
	ni <- grep("n$", names(dta), perl = T)

	#select rows
	coverage <- dta[yi] + dta[ni]
	if (type == "mean") {
		indices <- rowMeans(coverage) > aveDenom
	} 
	if (type == "allMin") {
		indices <- apply(coverage, 1, min) > aveDenom
	}

	ratio <- (dta[indices, yi] + ps) / (coverage[indices, ] + ps + ps)
	names(ratio) <- sub("y$", "", names(ratio), perl = T)

	if (verbose == T) {
		dtaSize <- dim(coverage)[1]
		fdtaSize <- dim(ratio)[1]
		s <- sprintf('dta N=%.f filtered dta N=%.f(%.2f%%) w type=%s aveDenom=%.f & ps=%.f', 
								 dtaSize, fdtaSize, (100 * fdtaSize / dtaSize), type, aveDenom, ps)		
		print(paste("TAG FILTERING:", s))
	}
	return(ratio)
}

test.mergeSetwDefault <- 
	function() {
		source("./bigData.R")
		# 0.
		default <- NULL
		set     <- list(a=10, b=20)
		res     <- mergeSetwDefault(set, default)
		res
		# 0.
		default <- list(a=10, b=20)
		set     <- NULL
		res     <- mergeSetwDefault(set, default)
		res
		# 1
		default <- list(a=10, b=20)
		set     <- list(a=9)
		res     <- mergeSetwDefault(set, default)
		res
		# 2 list && non-list
		default <- list(aa=11)
		set     <- "set"
		res     <- mergeSetwDefault(set, default)
		res
		# 2 non-list && list
		default <- "kk"
		set     <- list(aa=11)
		res     <- mergeSetwDefault(set, default)
		res
		# 3 nested
		default <- list(a=list(aa=1, bb=2, cc=list(ccc=2),  dd="a"),       b=20)
		set     <- list(a=list(aa=11,      cc=list(ccc2=4), dd=list(ddd=5))    )
		res     <- mergeSetwDefault(set, default)
		res  
		# 3 nested 2
		set     <- list(a=list(aa=1, bb=2, cc=list(ccc=2),  dd="a"),       b=20)
		default <- list(a=list(aa=11,      cc=list(ccc2=4), dd=list(ddd=5))    )
		res     <- mergeSetwDefault(set, default)
		res  
	}
# precedence: atomic > list > NULL
mergeSetwDefault <-
	function(list.s, list.d) {
		# base case -non-list
		if (is.atomic(list.s) && !is.null(list.s)) {
			res <- list.s
		} else {
			names.d <- names(list.d)
			names.s <- names(list.s)
			names.i <- intersect(names.d, names.s)
			names.d.unique <- setdiff(names.d, names.s)
			names.s.unique <- setdiff(names.s, names.d)
			res <- c(list.s[names.s.unique], list.d[names.d.unique])
			n <- names.i[1]
			for (n in names.i) {
				if (is.list(list.s[n][[1]]) || is.list(list.d[n][[1]])) {
					k <- mergeSetwDefault(list.s[n][[1]], list.d[n][[1]])
					res[n] <- list(k)
					#res[n] <- mergeSetwDefault(list.s[n][[1]], list.d[n][[1]])
				} else {
					res <- c(res, list.s[n])
				}
			}
		}
		res
	}

test.preprocess <-
	function() {
		source("bigData.R")
		dta <- iris[,1:4]
		head(dta)
		dim(dta)
		param <- NA 
		param <- NULL
		param <- list(melt=F, subset.cond="Species == 'setosa'", col.sel.pat="Length")
		param <- list(col.sel.pat="Length", ps=1, logf='log2')
		dta2 <- preprocess(dta=dta,param=param)
		param <- NULL
		dta2 <- preprocess(dta=dta, param=param)
		# RLE
		param.rle <- list(ps=1, logf='log2', col.median.center=T)
		dta2 <- preprocess(dta=dta[,1:4], param=param.rle)
		# round
		param <- list(ps=1, logf='log2', round=3, ceiling=T )
		dta2 <- preprocess(dta=dta[,1:4], param=param, verbose=T)
		head(dta2)
		dim(dta2)

		# header 
		param <- list(col.sel.pat='.srt.bam', colname.split.rank=2, colname.split.pat='.', mk.uniq.names=T, rm.col.sel.pat=".", sort.names=T)
		fle <- "~/Projects/biter_lingliu_DNAdamage/results/2015-10-12/EG_Liu_2015_2016_RNAseq/featureCounts.txt.gz"
		dta <- preprocess(fle, param=param)
		names(dta)

	}

# TODO add normalization
#	normFactor <- as.data.frame(matrix(rep(as.numeric(dta[1,-1]),dim(dta)[1]),nrow=dim(dta)[1], byrow=T))
#	dta[,-1] <- dta[,-1] / normFactor
# TODO replacement
#		dta$group <- gsub("(.+)_.+", "\\1", dta.m$variable)
#			dta$group <- gsub(group.col.pattern, "\\1", dta[,varname]); 
# TODO renaming
preprocess <- 
	function(fle, dta=NULL, param=NULL, verbose=F) {

		if (is.null(dta))  
			dta <- readBigTable(fle, header=T) 
		if (verbose) print(head(dta,2))

		if (any(grepl("partial",names(param))) && param$partial == T)
			dta.all <- dta

		for (p in names(param))  {
			pv <- param[p][[1]]
			if (verbose) print(c(p, pv))
			dta <- switch(p, 
										scale = { if (pv) {dta / data.frame(matrix(rep(colSums(dta), nrow(dta)),byrow=T,nrow=nrow(dta))) * 1000000 }},
										sort.names = { if (pv) { dta <- dta[sort(names(dta))] } },
										mk.uniq.names = { if (pv) { names(dta) <- make.names(names(dta), unique=T) }; dta},
										col.sel.pat = subset(dta, select=grep(pv, names(dta))),
										col.rm.pat = subset(dta, select=grep(pv, names(dta), invert=T)), # | for OR
#										rm.col.sel.pat = {names(dta) <- gsub(pv, "", names(dta), fixed=T); dta}, # does it make sense???
										colname.split.pat = { names(dta) <- sapply(strsplit(names(dta), pv, fixed=T), "[", param$colname.split.rank); dta; },
										melt = if (pv) { melt(dta) } else { dta },
										ps = (dta + param$ps),
										logf = do.call(pv, list(x=(dta))),
										min.rowMean = subset(dta, subset=rowMeans(dta)>pv),
										#subset.cond <- "value > 0.01"
										#subset.cond = subset(dta, subset=eval(parse(text=pv))),
										subset.cond = filter_(dta, pv),
										col.median.center = if (pv) { m <- apply(dta, 1, median); dta - m } else { dta },
										row.median.center = if (pv) { m <- apply(dta, 2, median); r <- t(t(dta) - m); rownames(r) <- rownames(dta); r } else { dta },
										col.mean.center = if (pv) { m <- apply(dta, 1, mean); dta - m } else { dta },
										row.mean.center = if (pv) { m <- apply(dta, 2, mean); r <- t(t(dta) - m); rownames(r) <- rownames(dta); r } else { dta },
										ceiling = if (pv) { ceiling(dta) } else { dta }, 
										round = round(dta, pv),
										t = if (pv) t(dta),
										dta)
			if (verbose) print(head(dta,1))
		}

		if (any(grepl("partial",names(param))) && param$partial == T && any(grepl("partial",names(param)))) {
			if (verbose) print(head(dta,1))
			if (verbose) print(head(dta.all,1))
			dta.all[,param$col.sel.pat] <- dta[, param$col.sel.pat] 
			dta <- dta.all
		}
		if (verbose) print(head(dta,1))

		dta

	}

# ratio (summarize) drops unused categories
my_transform <- 
	function(data, columns=1:ncol(data), dcolumns=NA, rep.sub='_(\\d+)$', filter=NA, 
					 ratio=F, parse.rowname=T, combine=F, round=F, long=T, 
					 verbose=T, ...) {

		data %<>% data.frame

		cdata_cols <- colnames(data)[columns]
		if (is.atomic(dcolumns) && is.na(dcolumns))
			ddata_cols <- unname(select_vars_(colnames(data), args=colnames(data), exclude=cdata_cols))
		else
			ddata_cols <- colnames(data)[dcolumns] 
		print(c("cdata_cols:::", cdata_cols))
		print(c("ddata_cols:::", ddata_cols))

		# conditional probability: P (ddata_cols-last | ddata_cols-notlast)
		if (ratio) {
			data %<>%
				group_by_(.dots=ddata_cols) %>% 
				summarize_each_("sum", vars=cdata_cols) %>% 
				mutate_each_(funs(./sum(.)), vars=cdata_cols) %>% 
				data.frame
		}
		if (parse.rowname) {
			data %<>% mutate(KeyR=rownames(.))
			data %<>%	mutate(NameR=gsub(rep.sub, "", data$KeyR), ReplicateR=getRepNumber(data$KeyR, rep.sub=rep.sub))
#			tidyr::extract(KeyR, c("NameR", "ReplicateR"), paste0('(.+)',rep.sub), remove=F)
		}
		# TODO replace DUMMY mutate w NSE: following does not work yet
		# https://cran.r-project.org/web/packages/dplyr/vignettes/nse.html
		# mutate(.dots=setNames(rowMeans(.[, grep(name, names(data))]), c(name))) 
		if (combine) { 
			group_names <- rle(gsub(rep.sub, "", cdata_cols))$values
			for (name in group_names) {
				coli <- grep(name, names(data))
				if (length(coli) > 1)
					data %<>%	mutate(DUMMY=rowMeans(.[, coli]))
				else
					data %<>%	mutate(DUMMY=(.[, coli]))
				if (round)
					data$DUMMY <- round(data$DUMMY)
				names(data)[ncol(data)] <- name
			}
			cdata_cols <- c(cdata_cols, group_names) #update
		}
		if (long) {
			data %<>% gather_("Key", "Value", cdata_cols)
				#tidyr::extract(Key, c("Name", "Replicate"), paste0('(.+)',rep.sub), remove=F)
			data %<>%	mutate(Name=gsub(rep.sub, "", data$Key), Replicate=getRepNumber(data$Key, rep.sub=rep.sub))
		}

		if (! is.na(filter))
			data %<>% filter_(filter)

		data %<>% data.frame

		if (verbose)
			print(head(data,25))

		data
	}

test.my_transform <- 
	function() {
		source("auxggplot2.R")
		data <- iris	
		names(data) <- c("sepal_1", "sepal_2", "petal_1", "petal_2", "Species")
		data.l <- my_transform(data, columns=1:4, filter="Species == 'setosa'", verbose=T)
		dim(data.l)
		data.l <- my_transform(data, columns=1:4, verbose=T)
		dim(data.l)

		# w/ rownames 
		data <- iris	
		names(data) <- c("sepal11", "sepal22", "sepal33", "petal222", "Species")
		my_transform(data, columns=1:4, combine=T, verbose=F, rep.sub="([[:digit:]]+)$")
		# combine
		data <- iris	
		names(data) <- c("sepal_1", "sepal_2", "sepal_3", "petal_22", "Species")
		data <- data.frame(cor(data[,1:4]))
		#data$dummy = c("a","a","b","b") 
		#data$dummy2 = c("aa","aa","bb","bb") 
		source("auxggplot2.R")
		data
		my_transform(data, columns=1:4, combine=T, verbose=F)
		my_transform(data, columns=1:4, combine=T, verbose=F, filter='is.na(Replicate)')
		my_transform(data, columns=1:4, combine=T, verbose=F, filter='!is.na(Replicate)')
		my_transform(data, columns=1:4, combine=T, round=T, verbose=F, filter='is.na(Replicate)')
		# ratio
		data <- iris	
		names(data) <- c("sepal_1", "sepal_2", "sepal_3", "petal_22", "Species")
		data <- data.frame(abs(cor(data[,1:4])))
		data <- my_transform(data, columns=1:4, ratio=F, parse.rowname=T, long=F, combine=F, verbose=T)
		source("auxggplot2.R")
		my_transform(data, columns=1:4, dcolumn=6, ratio=T, parse.rowname=F, long=F, combine=F, verbose=F)
		my_transform(data, columns=1:4, dcolumns=6, ratio=T, parse.rowname=F, long=T, combine=F, verbose=F)
		# ratio && others
		my_transform(data, columns=1:4, dcolumns=6, ratio=T, parse.rowname=F, long=T, combine=F, verbose=F)
		my_transform(data, columns=1:4, dcolumns=6, ratio=T, parse.rowname=T, long=T, combine=F, verbose=F)
		my_transform(data, columns=1:4, dcolumns=6, ratio=T, parse.rowname=T, long=T, combine=T, verbose=F)
	}

getRepNumber <- 
	function(s, rep.sub='_(\\d+)$') {
		sg <- as.character(sub(rep.sub,"", s))
		x <- substr(as.character(s), nchar(sg)+1, nchar(as.character(s)))
		x <- sub('[^[:digit:]]*', '', x)
		x[nchar(x) == 0 ] = NA
		x

#		if (all(sapply(sg, nchar) < sapply(as.character(s), nchar))) {
#			sub('[^[:digit:]]*', '', x)
#		} else {
#			s
#		}
	}

test.getRepNumber <- 
	function() {
		source("./auxggplot2.R")
		rep.sub <- "_(\\d*)$"
		s <- c("AAAa_1","hdhdh_WIDTH")
		getRepNumber(s, rep.sub=rep.sub)
		rep.sub <- '(\\d+)$'
		getRepNumber("Ct111rl_11", rep.sub=rep.sub)
		getRepNumber(c("Ct111rl_11","aa_WIDTH"), rep.sub=rep.sub) # not removed since pattern not found in all
		getRepNumber(c("Ct111rl_11","aa_22"), rep.sub=rep.sub)
		rep.sub <- "(\\d*)$"
		s <- c("AAAa1","hdhdh1")
		getRepNumber(s, rep.sub=rep.sub)
	}

