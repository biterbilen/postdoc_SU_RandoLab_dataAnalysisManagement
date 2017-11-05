library(edgeR)
library(EDASeq)
#library(rtracklayer)

source("~/Projects/biter_biter_shared/bin/bigData.R");
source("~/Projects/biter_biter_shared/bin/auxggplot2.R");
source("~/Projects/biter_biter_shared/bin/auxML.R");
source("~/Projects/biter_biter_shared/bin/auxGenomicFeatures.R");

test.caller.getSeqFeatures <- 
	function() {
		source("./auxslurm.R")
		d <- "~/Projects/biter_wosczyna_miRNAsInMuscleCellFate/data/Wosczyna_2015/"
		expressedTx.file <- paste0(d, "RNAseq/ER/expressedFeatureCounts.id") 
		expressedTx.file <- "head"
		mirSeq.file <- "~/Projects/biter_jbrett_epigenetics/data/miRBase/mature.fa.gz"
		miRExp.file <- "mmu-miR-206-3p"
		miRExp.file <- "~/Projects/biter_wosczyna_miRNAsInMuscleCellFate/data/Wosczyna_2013/selected3"
		exp.filter <- F
		regions.str <- "threeUTRsByTranscript"
		debug <- F
		gnm <- "mm10"
		mc.cores <- 1
		caller.getSeqFeatures(expressedTx.file=expressedTx.file, mirSeq.file=mirSeq.file, miRExp.file=miRExp.file, exp.filter=exp.filter,regions.str=regions.str,gnm=gnm,mc.cores=mc.cores, debug=debug)

		expressedTx.file="~/PI_HOME/Data/casco/UCSC_tracks/mm10/a.bed.gz"
		regions.str="self"
		expressedTx.colInd=1
		motifs.str="TFBS"
		caller.getSeqFeatures(expressedTx.file=expressedTx.file, regions.str=regions.str, expressedTx.colInd=expressedTx.colInd, motifs.str=motifs.str)

		# prtfa, OK 
		source("./auxslurm.R")
		d <- "~/Projects/biter_wosczyna_miRNAsInMuscleCellFate/data/Wosczyna_2015/"
		#expressedTx.file <- paste0(d, "RNAseq/DE/expressedFeatureCounts.id") 
		expressedTx.file <- paste0(d, "RNAseq/DE/edgeR_DE.txt.gz") 
		regions.str <- "promoters"; rfle <- NA
		regions.str <- "self"; rfle <- "~/Projects/biter_wosczyna_miRNAsInMuscleCellFate/data/Wosczyna_2015/RNAseq/DSWu2017BBAstateStrongActiveEnh/self.fa.gz"
		caller.getSeqFeatures(expressedTx.file=expressedTx.file, regions.str=regions.str, rfle=rfle, gmt=T, debug=T, tag="FF_Wosczyna_2015_RNAseq")

		# cdsmira 3utrmira, OK 
		d <- "~/Projects/biter_wosczyna_miRNAsInMuscleCellFate/data/Wosczyna_2015/"
		expressedTx.file <- paste0(d, "PullDown/DE/edgeR_DE.txt.gz") 
		mirSeq.file <- "~/Projects/biter_jbrett_epigenetics/data/miRBase/mature.fa.gz"
		miRExp.file <- NA
		miRExp.file <- "~/Projects/biter_wosczyna_miRNAsInMuscleCellFate/data/Wosczyna_2013/selected.txt.gz"
		#miRExp.file <- "mmu-miR-206-3p"
		regions.str <- "threeUTRsByTranscript"
		regions.str <- "cdsBy"
		debug <- F
		source("./auxslurm.R")
		caller.getSeqFeatures(expressedTx.file=expressedTx.file,mirSeq.file=mirSeq.file,miRExp.file=miRExp.file,exp.filter=F,regions.str=regions.str,gnm='mm10',mc.cores=1, gmt=T, debug=debug)

		# TODO old arguments to function caller.getSeqFeatures
		#expressedTx.file="~/Projects/biter_jbrett_epigenetics/results/2014-11-28/back.bs20_edgeRnorm_Cheung_A3_Q3/STAR_Cheung_A3_Q3_DE.txt",
		#mirSeq.file="~/PI_HOME/Data/casco/miRBase/mature.fa.gz",
		#miRExp.file="~/Projects/biter_jbrett_epigenetics/data/Cheung_2014_Nature/TaqManRodentmiRNA_QSC_ASC.txt.gz",

	}

caller.getSeqFeatures <- 
	function(expressedTx.file, exp.filter=T, mirSeq.file=NA, miRExp.file=NA, rfle=NA,
					 regions.str="promoters,threeUTRsByTranscript,cdsBy",
					 matrixSource="JASPAR2014", min.score="90%", seed.start=2, seed.end=8,
					 gnm="mm10", tag="", title="", 
					 gmt=F, mc.cores=1, debug=F, ...)
	{
		print(date())

		regions <- strsplit(regions.str, ",")[[1]]
		ofle.tag <- gsub(", ","_", regions.str)

		expressedTx <- readBigTable(expressedTx.file, header=T)

		aggt <- function(dta, tag) {
			data.frame(count=apply(dta, 1, function(x) length(x[x>0])), tag=tag)
		}

		title1 <- paste0(title,";FF:",ifelse(is.null(ofle.tag),"NULL",ofle.tag))
		pdta <- NULL
		for (r in regions) {
			if (is.na(rfle))
				rfle <- paste0("/share/PI/casco/Data/casco/refGeneSeq/",gnm,"/",gnm,"_refGene_canonical.v2.",r,".fa.gz")
			DNAStringSet <- readBigTable(rfle)
			if (r != "self")
				DNAStringSet <- DNAStringSet[names(DNAStringSet) %in% expressedTx$Id,] 
			if (r == "threeUTRsByTranscript" || r == "cdsBy") {
				mir.df <- get.revCompSeedsOfExpressedMirFamilies(mirSeq.file=mirSeq.file, mirExp.file=miRExp.file, exp.filter=exp.filter, seed.start=seed.start, seed.end=seed.end, mc.cores=mc.cores, ofle=paste0(r, "_revCompSeeds.txt.gz"))
				print(dim(mir.df))
				consensusStringList <- mir.df$seedRevComp 
				names(consensusStringList) <- mir.df$Geneid
				dta <- getSiteCountswConsensusStringAndDNAStringSet(consensusStringList,DNAStringSet,rev.comp=F,mc.cores=mc.cores)
				title1 <- paste0(title1,"\n;seedStart:",seed.start,";seedEnd=",seed.end)
				desc.tag <- paste0(r,"_seedS",seed.start,"E",seed.end,":")
			} else if (r == "promoters" || r == "self") {
				MatrixList <- getMatrixList(src=matrixSource, gnm, matrixtype="PWM", all_versions=T)
				dta <- getSiteCountswMatrixAndDNAStringSet(MatrixList,DNAStringSet,min.score=min.score,rev.comp=T,mc.cores=mc.cores)
				title2 <- paste0(title1,"\n;minScore=",min.score,";matrixSource=",matrixSource)
				desc.tag <- paste0(r,"_",matrixSource,":" )
			}
			dta$Id <- rownames(dta)
			cdta <- aggt(dta, r)
			print(head(cdta))
			title1 <- paste0(title1,";summary:", paste.namedList(summary(cdta$count)))
			pdta <- bind_rows(pdta, cdta)

			if (r != "self")
				dta %<>% inner_join(expressedTx %>% dplyr::select(Id, ENTREZID))
			writeData(dta %>% mutate(ENTREZID=NULL) %>% tibble::column_to_rownames("Id"), 
								paste0(ofle.tag, ".txt.gz"), row.names=T)
			if (gmt)
				if (r != "self") {
					writeData(dta %>% mutate(Id=NULL) %>% tibble::column_to_rownames("ENTREZID"),
										paste0(r,".gmt.gz"),frmt="gmt",desc.tag=desc.tag,row.names=F,col.names=F)
				} else {
					writeData(dta %>% mutate(Id=NULL),
										paste0(r,".gmt.gz"),frmt="gmt",desc.tag=desc.tag,row.names=F,col.names=F)
				}
		}

		print(title1)

		pdf(paste0(ofle.tag, ".pdf"), h=12, w=12)
		p <- my_bar_hist(pdta, aes(x=count,color=tag), 
											 facet="tag", long=F, stat="bin", title=title1)
		print(p)
		dev.off()

		print(sessionInfo())
		dta
	}

test.caller.GSEA <- function() {
	source("./auxslurm.R")

	d <- "~/Projects/biter_wosczyna_miRNAsInMuscleCellFate/data/Wosczyna_2015/"
	fle <- paste0(d, "RNAseq/DE/edgeR_DE.txt.gz") 
	#flec <- dir(d, pattern="PullDown/DE*/edgeR_DE.gmt.gz")
	flec <- Sys.glob(paste0(d,"/PullDown/DE*/*gmt.gz")); flec
	flec <- Sys.glob(paste0(d,"../Wosczyna_2015/PullDown/DE*/*gmt.gz")); flec
	flec <- Sys.glob(paste0(d,"/*/{DE,DEnulltotRNA}/*gmt.gz")); flec
	flec <- Sys.glob(paste0(d,"../Wosczyna_2015/PullDown/{DE,DEnulltotRNA}/*gmt.gz")); flec
	caller.GSEA(fle, flec=flec, plot2file=T)
	caller.GSEA(fle, flec=flec, plot2file=T, types=c("goana"), tag="GSEA_GO_wfisher")
	caller.GSEA(fle, flec=flec, plot2file=T, types=c("kegga"), tag="GSEA_KEGG_wfisher")
	caller.GSEA(fle, flec=flec, plot2file=T, types=c("customa"), tag="GSEA_CUS_wfisher")
	caller.GSEA(fle, flec=flet, plot2file=T, types=c("customa"), tag="GSEA_TF_wfisher")
	caller.GSEA(fle, flec=flem, plot2file=T, types=c("customa"), tag="GSEA_MIR_wfisher")

	top <- readBigTable(fle, header=T, row.names="Id")
	GSEA(top, species="Mm")
}

# taylored for edgeR top data frame 
caller.GSEA <-
	function(fle, title="", tag="GSEA", param=NULL, types=NULL, flec=NULL, 
					 plot2file=F) {

		top <- readBigTable(fle, header=T, row.names="Id")
		if (is.null(param))
			#param <- list(FDRcut=0.01, GSEApcut=0.01, species="Mm", stat="auto")
			param <- list(FDRcut=0.01, GSEApcut=0.05, species="Mm", stat="fisher")
		if (is.null(types))
			#types = c("customa")
			types = c("goana", "kegga", "customa")

		if (plot2file) 
			pdf(paste0(tag, ".pdf"), h=12, w=12)
		
		print("param:")
		print(param)
		print("flec:")
		print(flec)
		# dplyr discards rownames
		for (type in types) {
			top4GSEA <- 
				top %>% 
				tibble::rownames_to_column("Id") %>% 
				tibble::column_to_rownames("Id")
			print("top4GSEA:")
			print(dim(top4GSEA))
			GSEA(top4GSEA, FDRcut=param$FDRcut, type=type, 
					 flec=flec,
					 species=param$species, pcutoff=param$GSEApcut, stat=param$stat, 
					 ofle.tag=tag, title=title)
		} 

		if (plot2file)
			dev.off()
	}

customa <- 
	function(de, universe, species, gene.pathway, pathway.names) {
		kegga(de, universe=universe, species=species, gene.pathway=gene.pathway, pathway.names=pathway.names)
	}
# universe is restricted when FDRcut != NULL and universe is not set
GSEA <- 
	function(top, FDRcut=0.01, universe=NULL, keytype="ENTREZID",  
					 type="goana", stat="auto", flec=NULL,
					 species=NA, pcutoff=0.01, ofle.tag=NULL, 
					 valuetype="SYMBOL", title="") {
		#	library(KEGGREST)
		#	listDatabases()
		#	org <- keggList("organism")
		#	hsa <- keggList("hsa")
		#	query <- keggGet(c("hsa:10458"))
		if (is.na(species)) 
			stop("Set species!")

		ggplot <- set.ggplot()

		sorts <- list(updown=list(sort=NULL, eq="logFC>0 || logFC<=0"),
									up=list(sort="up", eq="logFC>0"), 
									down=list(sort="down", eq="logFC<=0"))

		# set functions
		if (type == "goana") {
			if (species == "Mm") {
				library(org.Mm.eg.db)
			} else if (species == "Hs") {
				library(org.Hs.eg.db)
			} else {
				stop(species, " library set")
			}
			GO2ALLEGS <- paste0("org.", species, ".egGO2ALLEGS")
			EG <- AnnotationDbi::toTable(get(GO2ALLEGS))
			EGn <- NULL
		} else if (type == "kegga") {
			species.KEGG <- switch(species, "Hs"="hsa", "Dm"="dme", "Mm"="mmu", "Rn"="rno", "Pt"="ptr")
			EG  <- getGeneKEGGLinks(species.KEGG) 
			EGn <- getKEGGPathwayNames(species.KEGG, remove=TRUE)
		} else {
			library(GSA)
			print("flec:")
			print(flec)
			flesa <- flec
			EG  <- data.frame(GeneID=NULL, PathwayID=NULL)
			EGn <- data.frame(PathwayID=NULL, Description=NULL)
			for (f in flesa) {
				EGt <- GSA.read.gmt(f)
				for (i in 1:length(EGt$genesets)) {
					EG <- bind_rows(EG, data.frame(GeneID=EGt$genesets[[i]], PathwayID=EGt$geneset.names[i]))
				}
				EGn <- bind_rows(EGn, data.frame(PathwayID=EGt$geneset.names, Description=EGt$geneset.descriptions))  
			}
		}

		unambig <- !duplicated(top[,keytype])

		if (is.null(universe))
			universe <- top %>% filter(unambig) %>% select_(keytype) 

		for (srt in names(sorts)) {
			print(srt)

			# reduce top keys to a list
			topl <- 
				top %>% 
				filter(unambig & FDR < FDRcut) %>% 
				filter_(sorts[[srt]][["eq"]]) %>% 
				dplyr::select_(keytype)

			topGSEA <- 
				do.call(type, list(de=topl[,1], universe=universe[,1], species=species, gene.pathway=EG, pathway.names=EGn)) %>%
				tibble::rownames_to_column(var="Id")

			if (stat == "fisher") {
				topGSEA$P.DE <-	
					lapply(1:nrow(topGSEA), 
								 function (i) fisher.test(matrix(c(topGSEA[i,"DE"], topGSEA[i,"N"], nrow(topl), nrow(universe)), nrow=2),
																					alternative="greater")$p.value)
				topGSEA %<>% mutate(P.DE=as.numeric(P.DE))
			}
			topGSEA %<>% arrange(P.DE) 

			head(topGSEA) %>% print

			if (type != "goana") 
				topGSEA %<>% rename(Term=Pathway)

			topGSEAwpcut <- topGSEA %>% filter(P.DE < pcutoff)

			names(EG)[1:2] <- c("gid", "GSid")
			EG.topl <- filter(EG, gid %in% topl[,1] & GSid %in% topGSEAwpcut$Id)

			getCommon <- 	
				function(GSid) {
					k <- top[,keytype] %in% unique(EG.topl[EG.topl$GSid %in% GSid, "gid"])
					k2 <- 
						filter(top, k) %>%
						dplyr::select_(valuetype) %>% 
						arrange_(valuetype)
					paste0(k2[,valuetype], collapse="|")
				}

			topGSEAwpcut$DEgenes <-	as.character(lapply(topGSEAwpcut$Id, getCommon))

			nP.DE <- nrow(topGSEAwpcut)
			topGSEA$pcutoff <- pcutoff
			topGSEA$rank <- 1:nrow(topGSEA)
			topGSEA$label <- paste0(topGSEA$rank,":",topGSEA$Term)
			title1 <- paste0(title,";tag:",ifelse(is.null(ofle.tag),"NULL",ofle.tag),";GSEA:",type,";species:",species,";sort=",srt,
											"\n;FDRcut=",FDRcut, ";universeN=",nrow(universe),";sigN=",nrow(topl),
											"\n;pcutoff=",pcutoff,";gsN=", nrow(topGSEA),";siggsN=",nP.DE)

			if (nrow(topGSEA) > 20) {
				# for ordering
				p <- my_bin(topGSEA, aes(x=log(N), y=DE/N, color=P.DE < pcutoff, label=label),
										textRepelFilter=paste0("rank < ", min(nP.DE,20)), smooth=F,
										long=F, slope=0, intercept=-1, size=rel(2), bins=nrow(topGSEA)) + labs(title=title1)
				print(p)
			}
			topGSEA$Term <- factor(topGSEA$Term, levels=topGSEA$Term)
			p <- my_bar_hist(topGSEA %>% head(min(20,nrow(topGSEA))),
											 aes(x=-log10(P.DE), y=Term, color=P.DE < pcutoff, 
													 label=paste0(";N=",N,";DE=",DE,";DE/N=",round(DE/N,2))), 
											 long=F, title=title1, size=rel(2))
			print(p)

			top$pcutoff <- NULL
			top$rank <- NULL
			top$label <- NULL

			if (!is.null(ofle.tag)) {
				commentLine <- gsub("\n",";", title1)
				writeData(topGSEAwpcut, 
									ofle=paste0(ofle.tag, "_", type, "_", srt, "_wcutoff.txt.gz"), 
									commentLine=commentLine, row.names=F)
			}
		}
	}



test.caller.DE.edgeR <- 
	function() {

		# Mike's data
		rep.sub <- '_(\\d+)$'
		d <- "~/Projects/biter_wosczyna_miRNAsInMuscleCellFate/data/Wosczyna_2015/RNAseq/ER/"
		fle <- paste0(d, "featureCounts.txt.gz") 
		efle <- paste0(d, "expressedFeatureCounts.id") 
		source("~/Projects/biter_biter_shared/bin/auxslurm.R")
		top <- caller.DE.edgeR(fle=fle, efle=efle, rep.sub=rep.sub, tag="RNAseq", plot2file = T)

		d <- "~/Projects/biter_wosczyna_miRNAsInMuscleCellFate/data/Wosczyna_2015/PullDown/ER/"
		fle <- paste0(d, "featureCounts.txt.gz") 
		efle <- paste0(d, "expressedFeatureCounts.id") 
		source("~/Projects/biter_biter_shared/bin/auxslurm.R")
		param <- list(out="top", scale="upperquartile", contrast=c(1,-1,-1,1), FDRcut=0.01, species="Mm")
		top <- caller.DE.edgeR(fle=fle, efle=efle, rep.sub=rep.sub, param=param, tag="PullDown", plot2file = T)

	}

#CURR EDA plot call for normalized counts
caller.DE.edgeR <-
	function(set, group, fle=NA, efle=NA, rep.sub='_(\\d+)$', col.sel.pat='.nsrt.bam', suf.pat=col.sel.pat,
					 param=NULL, tag="edgeR_DE", plot2file=F) {

		if (is.null(param))
			param <- list(out="top", scale="upperquartile", contrast=c(1,-1), FDRcut=0.01, species="Mm")

		print(paste("Doing caller.DE.edgeR param:", paste(param, collapse=" "), sep=""))

		if (plot2file) {
			pdf(paste(tag, ".pdf", sep=""), h=12, w=12)
		}

		ggplot <- set.ggplot()

		if (!is.na(fle)) {
			dta  <- readBigTable(fle, header=T)
			if (!is.na(efle)) {
				dta2 <- readBigTable(efle, header=T)
				dta %<>% inner_join(dta2, by="Geneid")
			}
			rownames(dta) <- dta$Geneid

			# select tag count fields and sort them
			dots <- sort(grep(col.sel.pat, colnames(dta), value=T))
			dta %<>% dplyr::select_(.dots=dots)
 			colnames(dta) <- gsub(paste0('.+\\.(.+)', suf.pat), '\\1', colnames(dta)) %>% 
				make.names(unique=T)

			column.labels <- colnames(dta)
			group <- sub(rep.sub, '', column.labels, perl=T)
			set <- newSeqExpressionSet(as.matrix(dta),
																 phenoData = data.frame(group, row.names=column.labels))
			print(set)
		}

		ret <- NULL
		# If multiple groups are present and topNinvert is set, get negative controls by correlation to 75% quantile
		if (param$out == "topNinvert" && length(unique(group)) > 2) {
			print("Correlation based negative control gene selection")
			# TODO add scaling based colSums
			dta <- counts(set)
			if (!any(is.na(normCounts(set))))
				dta <- normCounts(set)
			#cs <- colSums(dta)
			cs <- apply(dta, 2, quantile, probs=0.75)

			# negative controls are at the end; correlations are sorted in acsending order
			top <- data.frame(cor=sort(apply(dta, 1, function(x) { cor(cs, x) })))
			topN <- floor(nrow(top) *0.95) 
			p <- my_bar_hist(top, aes(x=cor), long=F, bins=100, fill="black", size=rel(2)) +
				labs(title="Pearson's R of gene expressions for negative control genes")
			print(p)

			top$negCon <- c(rep(F, topN), rep(T, nrow(dta) - topN ))

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
			print("contrast")
			print(contrast)

			lrt <- glmLRT(fit, contrast=contrast) # Likelihood ratio test
			top <- topTags(lrt, n=nrow(set))$table
			if (param$out == "topNinvert") {
				ret <- !(rownames(set) %in% rownames(top)[1:param$topN])
			} else if (param$out == "top") {

				dt <- data.frame(type=decideTestsDGE(lrt, adjust.method="BH", p.value=param$FDRcut)) %>% group_by(type) %>%
					count %>% data.frame 
				dt %<>% set_rownames(dt$type)

				print(paste("DE at FDR<", param$FDRcut, sep=""))
				print(dt)

				#plotQLDisp(fit)
				print(plotBCV(y))		
				# TODO what is it plotting?
				plotMDS(y) # noprint
				#plotMDS(y,col=as.numeric(targets$Genotype)) # noprint
				title <- paste(tag, paste.namedList(param),
											 paste0(colnames(lrt$design), collapse=","), 
											 paste.namedList(c(down=unname(dt["-1","n"]), up=unname(dt["1","n"]), all=sum(dt[,"n"])), sp="="),
											 sep="\n")

				top <- annotate(top, param$species)

				# TODO plot as layers
				#				http://stackoverflow.com/questions/15706281/controlling-order-of-points-in-ggplot2-in-r
				top$FDRcut <- param$FDRcut # patch for ggplot
				top$rank <- 1:nrow(top)
				top$label <- paste0(top$rank,":", top$SYMBOL)
				N <- nrow(top[top$FDR < param$FDRcut,])
				p <- my_bin(top, 
										aes(x=logCPM, y=logFC, size=-log10(FDR), color=FDR<FDRcut, label=label),
										textRepelFilter=paste0("rank < ", min(N, 60)), 
										long=F, slope=0, intercept=0, size=rel(2), bins=75) + labs(title=title)
				print(p)
				p <- my_bin(top, 
										aes(x=logFC, y=-log10(FDR), color=FDR<FDRcut, label=label),
										textRepelFilter=paste0("rank < ", min(N, 60)), smooth=F, 
										long=F, slope=0, intercept=0, size=rel(2), bins=75) + labs(title=title)
				print(p)
				top$rank <- NULL
				top$FDRcut <- NULL
				top$label <- NULL
					
				commentLine <- gsub("\n",";", title)

				writeData(top, ofle=paste(tag, ".txt.gz", sep=""), 
									commentLine=commentLine, row.names=T)

				ret <- top
			}
		} else {
			stop(param$out, " is not defined in function caller.DE.edgeR")
		}

		if (plot2file) 
			dev.off()

		ret
	}


test.get.EDAfeatureCounts <- 
	function() {
		source("~/Projects/biter_biter_shared/bin/auxslurm.R")
		frac <- 0.05
		fle <- "~/Projects/biter_lingliu_DNAdamage/results/2017-03-13/FC_Liu_2017_ATACseq/featureCounts.txt.gz"
		get.EDA4featureCounts(fle, plotfrac=frac, ofle.tag='his2_featureCounts')

		frac <- 1
		fle <- "~/Projects/biter_wosczyna_miRNAsInMuscleCellFate/results/2017-04-05/FC_Wosczyna_2015_PullDown/featureCounts.txt.gz"
		ddata_cols <- "ERCC"
		get.EDA4featureCounts(fle, plotfrac=frac, ddata_cols=ddata_cols)
		get.EDA4featureCounts(fle, plotfrac=frac, ddata_cols=ddata_cols, ofle.tag='mir2_featureCounts')
		get.EDA4featureCounts(fle, plotfrac=frac, ddata_cols=ddata_cols, ofle.tag='mir2_featureCountsGMM2', GMM2expressed=T, plotRatio=F)

		frac <- 1
		fle <- "~/Projects/biter_demorree_Fosl1/results/2017-06-26/FC_Wold_2011_ChIPseq/featureCounts.txt.gz"
		ddata_cols <- "ERCC"
		get.EDA4featureCounts(fle, plotfrac=frac, ddata_cols=ddata_cols)
		#get.EDA4featureCounts(fle, frac=frac, ofle.tag='fosl1_featureCounts')

	}

get.EDA4featureCounts <- 
	function(fle, plotfrac=1, ofle.tag=NA, ddata_cols=NA, col.sel.pat='.nsrt.bam', 
					 ps=1, plotRatio=F, GMM2expressed=F, trimmedMeanScaled=F, trim=0.2,
#					 UQscaled=F, probs=0.75, 
					 ...) {

    data <- readBigTable(fle, header=T, row.names='Geneid')
		ids <- rownames(data)

		# rename unify column names before dplyr operations
		cdata_cols <- grep(col.sel.pat, colnames(data)) #index
#colnames(data)[cdata_cols[1]] <- colnames(data)[cdata_cols[8]] #change name
		colnames(data)[cdata_cols] <- 
 			gsub(paste0('.+\\.(.+)', col.sel.pat), '\\1', colnames(data)[cdata_cols]) %>%
			make.names(unique=T)

		# discrete data MAKE everything factor
		if (!is.na(ddata_cols[1]))
			for (f in ddata_cols) 
				if (f == "ERCC") 
					data %<>% mutate(ERCC=as.factor(paste0("ERCC=", grepl("ERCC", data$Chr))))

		# sort cdata_cols and ddata_cols
		dots <- sort(colnames(data)[cdata_cols])
		if (!is.na(ddata_cols))
			dots <- c(ddata_cols, dots)

		# select relevant fields
		data %<>% dplyr::select_(.dots=dots)
		columns <- (1+length(ddata_cols)):ncol(data)

		if (!is.na(ofle.tag)) {
			pdf(paste(ofle.tag,".pdf",sep=""), h=12, w=12)
		}

		title <- paste0("N=", nrow(data))
		dataa <- NA
		if (!is.na(ddata_cols)) {
			dataa <- my_transform(data, columns=columns, ratio=T, parse.rowname=F, long=F, combine=F, verbose=T)
			print("data.ratio")
			print(dataa)
		}
		if ((!is.atomic(dataa)) && nrow(dataa) > 1 && plotRatio) { 
			(my_boxplot(dataa, mapping=aes(x=Name, color=Name, shape=Name, label=Replicate),
									columns=columns,  
									panel="dot", facet=ddata_cols[1], ncol=1, size=rel(6)) + labs(y="Ratio", title=title)) %>% print
			dataa <- dataa[ddata_cols]
		}

		# log transformation
		data %<>% mutate_at(columns, funs(log(.+ps)))
		title <- paste0(title, ";ps=", ps, ";logf:", "log");

		if (GMM2expressed) {
			# plot
			param.Mclust=list(G=2) # for log transformed data
			datac <- mGMM(data, columns=columns, param.Mclust=param.Mclust)
			# write expressed genes to file
			idse <- data.frame(Geneid=ids[datac$NotClass1==2])
			writeData(idse, ofle=paste0(ofle.tag, ".id"), row.names=F)
			print(paste("Expressed genes:\t", length(idse)))
			data <- bind_rows((data %>% mutate(NotClass1=as.factor(0))), datac)

			my_dens(data,
							mapping=aes(color=NotClass1), 
							columns=columns,
							facet="Key", size=rel(2), title=title, nrow=2) %>% print 
			data <- filter(datac, NotClass1==2) %>% dplyr::select(-NotClass1)
			title <- paste(title,";N=",nrow(data), sep="")
		}

		if (trimmedMeanScaled) {			
			my_mean <- function(x) mean(x, trim=trim)	
			data %<>% mutate_at(columns, funs(. - my_mean(.)))
			title <- paste0(title, ";trimmedMeanCentered_trim=", trim)
		}
#		if (UQscaled) {			
#			my_UQ <- function(x) quantile(x, probs=probs)	
#			data %<>% mutate_at(columns, funs(. - my_UQ(.)))
#			title <- paste0(title, ";UQCentered")
#		}

		title <- paste0(title, ";plotfrac=", plotfrac)

		if (plotfrac < 1) {
			data %<>% sample_frac(plotfrac)
		}

		print(paste("Doing", title))
		get.EDA(data, dataa=dataa, columns=columns, size=rel(2), title=title, ...)

		if (!is.na(ofle.tag))
			dev.off()
	}	


