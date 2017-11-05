#calculates the mRNA fraction of a makerefmetdfb-created dataframe (calcmrnafrageneral is more accomomdating of other input files)
test.calcmrnafrageneral <- 
	function() {
	source("erccdashboard.R")
	#source("~/Projects/biter_biter_shared/bin/auxlattice.R")
	#ercc <- readBigTable("ercc.txt", header=T, sep="\t")
	source("~/Projects/biter_biter_shared/bin/bigData.R")
	dta <- readBigTable("./FC_Wosczyna_2015_PullDown/featureCounts.txt.gz", header=T)
	names(dta) <- gsub(".sorted.bam", "", names(dta))
	names(dta) <- gsub("RM.", "", names(dta))
	ei <- grepl("ERCC", dta$Chr)
	rownames(dta) <- c(as.character(dta[!ei, "Geneid"]), as.character(dta[ei, "Chr"]))
	head(dta)
	cols <- grep("FAP", names(dta))
	dat <- as.matrix(dta[,cols])
	head(dat)
	spikeID <- "ERCC-"
	spikemassfraction <- c(1e-6, 1e-6, 1e-5, 1e-5, 1e-6, 1e-6, 1e-5, 1e-5)
	countcolumns <- 1:dim(dat)[2]
	annotcolumn <- NA
	head(dat)
	dim(dta)
	mfrac <- calcmrnafracgeneral(dat, spikeID, spikemassfraction, countcolumns) 
	mfrac <- calcmrnafracgeneral(dat, spikeID=spikeID, spikemassfraction=spikemassfraction) 
	mfrac
}

calcmrnafracgeneral <- 
	function(dat, spikeID="ERCC-", spikemassfraction=.1,  countcolumns=NULL)
	{
		if (is.null(countcolumns)) 
			countcolumns <- 1:dim(dat)[2]

		#one way to identify spikes, if row names = spikeID
		ercc <- rownames(dat)[which(substr(rownames(dat),1,5)==spikeID)] 

		nonercc <- !(rownames(dat) %in% ercc)
		#determines the counts for spikes and non-spikes.
		count <- rbind(colSums(dat[ercc,countcolumns]),
									 colSums(dat[nonercc,countcolumns])) 
		#defines the "targeted" mass fraction for spikes : Either a vector with length = #columns,or a scalar
		ercc.targ <- spikemassfraction 
		mRNA.frac <- ercc.targ*count[2,]/count[1,] #calculates an mRNA fraction based on those available data
		#this part doesn't normalize to one, but that's not exactly complicated.
		return(mRNA.frac)
	}

# TODO merge with ERCC file
test.GeneralLMest <- 
	function() {
	source("erccdashboard.R")
	source("~/Projects/biter_biter_shared/bin/bigData.R")
	ercc <- readBigTable("ercc.txt", header=T, sep="\t")
	head(ercc)
	dta <- readBigTable("./FC_Wosczyna_2015_PullDown/featureCounts.txt.gz", header=T)
	names(dta) <- gsub(".sorted.bam", "", names(dta))
	names(dta) <- gsub("RM.", "", names(dta))
	ei <- grepl("ERCC", dta$Chr)
	rownames(dta) <- c(as.character(dta[!ei, "Geneid"]), as.character(dta[ei, "Chr"]))
	head(dta)
	cols <- grep("FAP", names(dta))
	dat <- dta[,cols]
	head(dat)
	infile <- dat
	spikeID <- "ERCC-"
	components <- c("FAP_mir206Tx_MWmodPD","FAP_mir206Tx_QiagenRNEasy","FAP_NCTx_MWmodPD","FAP_NCTx_QiagenRNEasy")
	mixes <- "Mix.1"
	spikemassfraction <- c(1e-6, 1e-6, 1e-5, 1e-5, 1e-6, 1e-6, 1e-5, 1e-5)
	}

GeneralLMest <- 
	function(infile=dat, spikeID="ERCC-",components=c("FAP_mir206Tx_MWmodPD","FAP_mir206Tx_QiagenRNEasy","FAP_NCTx_MWmodPD","FAP_NCTx_QiagenRNEasy"), mixes=c("Mix.1"), spikemassfraction=c(1e-6, 1e-6, 1e-5, 1e-5, 1e-6, 1e-6, 1e-5, 1e-5)) {
		##formatting of infile:
		# Rownames: gene <- ids, with a phrase ("spikeID") that identifies Spike-In controls
		# Columns: 1 column per component, 1 column per mix
		#must input list of components and list of mixes: For example: the BLM mix contains components= c("bep","lep","mep") mixes=c("a1","a2")
		#spikemassfraction is the mass proportion you targeted with your spike-in controls
		#input should already be averaged across replicates unless i want to try to find a way to guess...
		#I don't know what else to do with replicates right now other than average them, but i imagine i'll come up with something more interesting later.
		#find all replicates and turn them into means
		RepMeans <- data.frame(matrix(0, dim(infile)[1], length(components)))
		for(i in 1:length(components)) {
			replicates <- grep(components[i],colnames(infile))
			RepMeans[,i] <- rowMeans(infile[,replicates])
		}
		rownames(RepMeans) <- rownames(infile)
		colnames(RepMeans) <- components

		#testphase code trying to handle SEQCDFs and arbitrary ones with replicates. Because I really need to spend a day making this code extensible instead of just remaking the figure in a more ad-hoc way...
		mfrac <- calcmrnafracgeneral(infile, spikeID, spikemassfraction)

		mixval <- NULL;
		rdf    <- NULL
		#convincing LM to handle arbitrary columns is not simple: The following block tries to do that.
		for(mixtext in mixes) {
			modeltext <- NULL
			for(i in 1:length(components)){
				modeltext <- paste(modeltext,(paste("+I(",mfrac[I],"*", components[i],")+0")),sep=" ")
				print(modeltext)
			}
			typeof(infile)
			lm(data=infile, as.formula(paste(mixtext,"~",modeltext)))
			mixval <- coefficients(lm(data=infile, as.formula(paste(mixtext,"~",modeltext))))/sum(coefficients(lm(data=infile,as.formula(paste(mixtext,"~",modeltext)))))
			rdf<-rbind(rdf,mixval)
		}
		rownames(rdf) <- mixes
		colnames(rdf) <- components
		return(rdf)
	}

do.EDA <- function(fle="./FC_Wosczyna_2015_PullDown/featureCounts.txt.gz", ps = 1) {
	source("erccdashboard.R")
	source("~/Projects/biter_biter_shared/bin/auxML.R")

	dta <- caller.getExpressed(fle=fle)

	return(dta)

	ercc <- readBigTable(erccfle, header=T, sep="\t")
	dta.ercc <- merge(dta, ercc, by.x="Geneid", by.y="ERCC.ID")[cols]
	dta.ercc.back <- dta.ercc
	cols.counts <- c(grepl("FAP", names(dta.ercc)))
	dta.ercc[,cols.counts] <- dta.ercc[,cols.counts] + ps
	cols.4log2 <- c(grep("FAP|concentration", names(dta.ercc),value=T))
	dta.ercc[,cols.4log2] <- log2(dta.ercc[,cols.4log2])
	cols <- c("Length", grep("FAP|concentration", names(dta.ercc),value=T))
	get.exploratoryAnalysisPlots(dta.ercc[,cols], ofle=paste("EDA_expressedGene_ERCCreadCounts",sep=""))
	dim(dta.ercc)
	head(dta.ercc)

	# Negative binomial fitting of log2 counts data with concentration and Length: Length does not have hgh predictive value
	library(MASS)
	nb1 <- glm.nb(ceiling(FAP_mir206Tx_QiagenRNEasy_1) ~ concentration.in.Mix.1..attomoles.ul. + Length, data = dta.ercc, subset=FAP_mir206Tx_QiagenRNEasy_1>0)
	warnings(nb1)
	summary(nb1)
	nb2 <- glm.nb(ceiling(log2(FAP_mir206Tx_QiagenRNEasy_1)) ~ concentration.in.Mix.1..attomoles.ul. + Length, data = dta.ercc.back, subset=FAP_mir206Tx_QiagenRNEasy_1>0)
	summary(nb2)
	nb3 <- glm.nb(ceiling(log2(FAP_mir206Tx_QiagenRNEasy_1)) ~ concentration.in.Mix.1..attomoles.ul., data = dta.ercc.back, subset=FAP_mir206Tx_QiagenRNEasy_1>0)
	summary(nb3)
	anova(nb1, nb2, nb3)
	anova(nb2, nb3)

	plot(predict(nb3, type="response"), residuals(nb, type="deviance"))

}


polishData <- function() {
	source("erccdashboard.R")
	source("./bigData.R")
	source("./auxlattice.R")
	ercc <- readBigTable("ercc.txt", header=T, sep="\t")
	dta <- readBigTable("./FC_Wosczyna_2015_PullDown/featureCounts.txt.gz", header=T)
	names(dta) <- gsub(".sorted.bam", "", names(dta))
	names(dta) <- gsub("RM.", "", names(dta))
	dim(dta)
	head(dta)

	dta.ercc <- subset(dta, subset=grepl("ERCC", Chr), select=grep("FAP", names(dta), value=T), drop=T)
	dim(dta.ercc)
	hist(log2(dta.ercc$FAP_mir206Tx_QiagenRNEasy_1), n=20)
	dta.nercc <- subset(dta, subset=!grepl("ERCC", Chr), select=grep("FAP", names(dta), value=T), drop=T)
	dim(dta.nercc)
	hist(log2(dta.nercc$FAP_mir206Tx_QiagenRNEasy_1), n=20)
	tot <- rbind(colSums(dta.ercc), colSums(dta.nercc))
	apply(tot, 2, function(x) { x[1]/x[2]; } )

	filter <- apply(zfGenes, 1, function(x) length(x[x>5])>=2)
	filtered <- zfGenes[filter,]

	dta.s <- subset(dta, subset=grepl("ERCC", Chr), select=c("Chr","Length",grep("FAP",names(dta), value=T)), drop=T)
	hist(log2(dta.s$FAP_mir206Tx_MWmodPD_1), n=20)
	hist(log2(dta.s$FAP_mir206Tx_QiagenRNEasy_1), n=20)
	hist(log2(dta.s$concentration.in.Mix.1..attomoles.ul.), n=20)
	dim(dta.s)
	head(dta.s)
	i <- (rowSums(dta.s[,-c(1,2)] > 0)+0) > (dim(dta.s)[2]/2)
	i <- (rowSums(dta.s[,-c(1,2)] > 0)+0) > 0
	sort(rowSums(dta.s[,-c(1,2)] > 0)+0)
	i
	dta.s <- merge(dta.s[i,], ercc, by.x="Chr", by.y="ERCC.ID")
	dim(dta.s)
	head(dta.s)

	# negative binomial
	head(dta.s[,c(13,5)])
	plot(log2(dta.s[,c(13,5)] + 1))

	library(edgeR)
	summary(m1 <- glm.nb(daysabs ~ math + prog, data = dat))
	hist(log2(dta.s$FAP_mir206Tx_QiagenRNEasy_2), n=20)

	nb <- glm.nb(FAP_mir206Tx_QiagenRNEasy_1 ~ concentration.in.Mix.1..attomoles.ul. + Length, data = dta.s, subset=FAP_mir206Tx_QiagenRNEasy_1>0)
	summary(nb)
	warnings()
	nb1 <- glm.nb(FAP_mir206Tx_QiagenRNEasy_1 ~ concentration.in.Mix.1..attomoles.ul. + Length, data = dta.s, subset=FAP_mir206Tx_QiagenRNEasy_1>0, link="log")
	warnings()
	summary(nb1)
	nb2 <- glm.nb(ceiling(log2(FAP_mir206Tx_QiagenRNEasy_1)) ~ concentration.in.Mix.1..attomoles.ul. + Length, data = dta.s, subset=FAP_mir206Tx_QiagenRNEasy_1>0)
	summary(nb2)
	nb3 <- glm.nb(ceiling(log2(FAP_mir206Tx_QiagenRNEasy_1)) ~ concentration.in.Mix.1..attomoles.ul., data = dta.s, subset=FAP_mir206Tx_QiagenRNEasy_1>0)
	summary(nb3)

	head(dta.s)
	par(mfrow=c(2,2))
	plot(dta.s$Length, log2(dta.s$FAP_mir206Tx_QiagenRNEasy_1))
	plot(log2(dta.s$concentration.in.Mix.1..attomoles.ul.), log2(dta.s$FAP_mir206Tx_QiagenRNEasy_1))
	summary(nb2)
	warnings()
	nb1
	summary(nb)
	coef(nb)
	plot(predict(nb, type="response"), residuals(nb, type="deviance"))
	hatvalues(nb)
	library(car)
	influencePlot(nb)
	warnings()

	ps <- 1
	cols <- grep("FAP|attomoles", names(dta.s), value=T, perl=T) 
	dta.p <- dta.s[,cols]
	cols <- grep("FAP", names(dta.p), value=T, perl=T) 
	dta.p[,cols] <- dta.p[,cols] + ps
	dta.p <- log2(dta.p)
	dta.p
	head(dta.p)
	get.exploratoryAnalysisPlots(dta.p, "ERCC_log2_count")
	layout(matrix(c(2,2),2,2,byrow=T))
	par(mar = c(2,2,1,1))
	dev.new()
	par(mfrow=c(2,4))
	plot(dta.p[,c(9,1)])
	plot(dta.p[,c(9,2)])
	plot(dta.p[,c(9,3)])
	plot(dta.p[,c(9,4)])
	plot(dta.p[,c(9,5)])
	plot(dta.p[,c(9,6)])
	plot(dta.p[,c(9,7)])
	plot(dta.p[,c(9,8)])
	dev.off()

	head(dta.p)
	#dta.p[,cols] <- dta.p[,cols] / dta.s[,"Length"] 
	#dta.p
	#get.exploratoryAnalysisPlots(dta.p, "ERCC_log2_count_per_length")

	library(MASS)
	head(dta.p)
	head(quine, 30)
	hist(quine$Days, n=50)
	quine.nb1 <- glm.nb(Days ~ Sex/(Age + Eth*Lrn), data = quine)
	quine.nb2 <- update(quine.nb1, . ~ . + Sex:Age:Lrn)
	quine.nb3 <- update(quine.nb2, Days ~ .^4)
	anova(quine.nb1, quine.nb2, quine.nb3)
}

test.ruvseq <- function() {
	fle <- "/home/biter/Projects/biter_wosczyna_miRNAsInMuscleCellFate/results/2015-07-17/FC_Wosczyna_2015_PullDown/featureCounts.txt.gz"
	spikeIn <- list(spikeId="ERCC")
	ps <- 1

}

ruvseq <- function(fle="./FC_Wosczyna_2015_PullDown/featureCounts.txt.gz", spikeIn=list(spikeId="ERCC"), ps=1) {
	## ----style-knitr, eval=TRUE, echo=FALSE, results="asis"------------------
	BiocStyle::latex()
	source("~/Projects/biter_biter_shared/bin/auxML.R")
	dta <- caller.getExpressed(fle=fle, spikeIn=spikeIn)

	head(dta)

	## ----setup, echo=FALSE---------------------------------------------------
	library(knitr)
	opts_chunk$set(dev="pdf", fig.align="center", cache=TRUE, message=FALSE, out.width=".49\\textwidth", echo=TRUE, results="markup", fig.show="hold")
	options(width=80)


	## ----data----------------------------------------------------------------
	library(RUVSeq)
	library(zebrafishRNASeq)
	data(zfGenes)
#	head(zfGenes)
#	tail(zfGenes)


	## ----filter--------------------------------------------------------------
	filter <- apply(zfGenes, 1, function(x) length(x[x>5])>=2)
	filtered <- zfGenes[filter,]

	filtered <- data.frame(dta[,-c(1,2)], row.names=dta$Geneid)
	i <- grepl(spikeIn$spikeId, dta$Geneid)
	genes <- rownames(filtered)[!i]
	spikes <- rownames(filtered)[i]

	## ----store <- data----------------------------------------------------------
	x <- as.factor(gsub("[0-9]*$","", names(filtered), perl=T))
	set <- newSeqExpressionSet(as.matrix(filtered),
														 phenoData = data.frame(x, row.names=colnames(filtered)))

	## ----rle, fig.cap="No normalization.",fig.subcap=c("RLE plot","PCA plot")----
	library(RColorBrewer)
	colors <- brewer.pal(length(unique(x)), "Set2")
	colors <- brewer.pal(length(unique(x)), "Paired")

	plotRLE(set, outline=FALSE, ylim=c(-4, 4), col=colors[x], las=2)
	plotPCA(set, k=4, col=colors[x], cex=1.2, las=2)
	get.pcaplot(filtered, col=colors, pch=19)

	library(yeastRNASeq)
	data(geneLevelData) 
	mat <- as.matrix(geneLevelData) 
	data <- newSeqExpressionSet(mat,phenoData=AnnotatedDataFrame(data.frame(conditions=factor(c("mut", "mut", "wt", "wt")),row.names=colnames(geneLevelData))))
	prc <- princomp(filtered)
	EDASeq::plotPCA(data, col=rep(1:2, each=2), k=3)
	get.exploratoryAnalysisPlots(dta=log2(as.data.frame(counts(set))+ps), ofle=paste("set_counts",sep=""))
	
	## ----uq, fig.cap="Upper-quartile normalization.", fig.subcap=c("RLE plot","PCA plot")----
	set <- betweenLaneNormalization(set, which="upper")
	plotRLE(set, outline=FALSE, ylim=c(-4, 4), col=colors[x], las=2)
	plotPCA(set, k=3, col=colors[x], cex=1.2)
	get.pcaplot(filtered, comp="rotation", col=colors, pch=19)
	get.exploratoryAnalysisPlots(dta=log2(as.data.frame(normCounts(set))+ps), ofle=paste("set_counts_upper",sep=""))


	## ----ruv <- spikes, fig.cap="RUVg normalization based on spike-in controls.", fig.subcap=c("RLE plot","PCA plot")----
	set1 <- RUVg(set, spikes, k=2)
	pData(set1)
	plotRLE(set1, outline=FALSE, ylim=c(-4, 4), col=colors[x],las=2)
	plotPCA(set1, col=colors[x], cex=1.2)


	## ----edger, eval=FALSE---------------------------------------------------
	design <- model.matrix(~x + W_1, data=pData(set1))
	y <- DGEList(counts=counts(set), group=x)
	y <- calcNormFactors(y, method="upperquartile")
	y <- estimateGLMCommonDisp(y, design)
	y <- estimateGLMTagwiseDisp(y, design)
	
	fit <- glmFit(y, design)
	lrt <- glmLRT(fit, coef=2)
	 
	topTags(lrt)


	## ----empirical-----------------------------------------------------------
	design <- model.matrix(~x, data=pData(set))
	y <- DGEList(counts=counts(set), group=x)
	y <- calcNormFactors(y, method="upperquartile")
	y <- estimateGLMCommonDisp(y, design)
	y <- estimateGLMTagwiseDisp(y, design)

	fit <- glmFit(y, design)
	lrt <- glmLRT(fit, coef=2)

	top <- topTags(lrt, n=nrow(set))$table
	empirical <- rownames(set)[which(!(rownames(set) %in% rownames(top)[1:5000]))]


	## ----emp <- ruvg, fig.cap="RUVg normalization based on empirical controls.", fig.subcap=c("RLE plot","PCA plot")----
	set2 <- RUVg(set, empirical, k=2)
	pData(set2)
	plotRLE(set2, outline=FALSE, ylim=c(-4, 4), col=colors[x])
	plotPCA(set2, col=colors[x], cex=1.2)


	## ----diff----------------------------------------------------------------
	differences <- matrix(data=c(1:4, 5:8), byrow=F, nrow=2)
	differences


	## ----ruvs, eval=FALSE----------------------------------------------------
	set3 <- RUVs(set, genes, k=1, differences)
	pData(set3)


	## ----res, eval=FALSE-----------------------------------------------------
	design <- model.matrix(~x, data=pData(set))
	y <- DGEList(counts=counts(set), group=x)
	y <- calcNormFactors(y, method="upperquartile")
	y <- estimateGLMCommonDisp(y, design)
	y <- estimateGLMTagwiseDisp(y, design)
	 
	fit <- glmFit(y, design)
	res <- residuals(fit, type="deviance")


	## ----ruvr, eval=FALSE----------------------------------------------------
	set4 <- RUVr(set, genes, k=1, res)
	pData(set4)


	## ----sessionInfo, results="asis"-----------------------------------------
	toLatex(sessionInfo())

}


erccdashboard <- function() {

### R code from vignette source 'erccdashboard.Rnw'

###################################################
### code chunk number 1: erccdashboard.Rnw:15-19
###################################################
options(width=60, continue = "  ")
#options(SweaveHooks=list(fig=function()
#                par(mar=c(5.1,4.1,1.1,2.1))))
library( "erccdashboard" )


###################################################
### code chunk number 2: loadExampleData
###################################################
data(SEQC.Example)


###################################################
### code chunk number 3: defineInputData
###################################################
datType = "count" # "count" for RNA-Seq data, "array" for microarray data
isNorm = FALSE # flag to indicate if input expression measures are already
# normalized, default is FALSE 
exTable = MET.CTL.countDat # the expression measure table
filenameRoot = "RatTox" # user defined filename prefix for results files
sample1Name = "MET" # name for sample 1 in the experiment
sample2Name = "CTL" # name for sample 2 in the experiment
erccmix = "RatioPair" # name of ERCC mixture design, "RatioPair" is default
erccdilution = 1/100 # dilution factor used for Ambion spike-in mixtures
spikeVol = 1 # volume (in microliters) of diluted spike-in mixture added to 
#   total RNA mass
totalRNAmass = 0.500 # mass (in micrograms) of total RNA 
choseFDR = 0.05 # user defined false discovery rate (FDR), default is 0.05



###################################################
### code chunk number 4: inspectRatCount
###################################################
head(MET.CTL.countDat)
dim(MET.CTL.countDat)


###################################################
### code chunk number 5: runDashboardRatcount
###################################################
library(gtable)
exDat <- runDashboard(datType="count", isNorm = FALSE,
											exTable=MET.CTL.countDat,
											filenameRoot="RatTox", sample1Name="MET",
											sample2Name="CTL", erccmix="RatioPair",
											erccdilution=1/100, spikeVol=1,
											totalRNAmass=0.500, choseFDR=0.1)


###################################################
### code chunk number 6: initializeData
###################################################
summary(exDat)


###################################################
### code chunk number 7: ratPlotA
###################################################
grid.arrange(exDat$Figures$dynRangePlot)


###################################################
### code chunk number 8: ratPlotB
###################################################
grid.arrange(exDat$Figures$rocPlot)


###################################################
### code chunk number 9: ratPlotC
###################################################
grid.arrange(exDat$Figures$lodrERCCPlot)


###################################################
### code chunk number 10: ratPlotD
###################################################
grid.arrange(exDat$Figures$maPlot)


###################################################
### code chunk number 11: SEQCMainArray
###################################################
exDat <- runDashboard(datType="array", isNorm = FALSE,
											exTable=UHRR.HBRR.arrayDat,
											filenameRoot = "Lab13.array",
											sample1Name = "UHRR", sample2Name="HBRR",
											erccmix = "RatioPair", erccdilution = 1, 
											spikeVol = 50, totalRNAmass = 2.5*10^(3), choseFDR=0.01)


###################################################
### code chunk number 12: viewDashboardOrder
###################################################
runDashboard


###################################################
### code chunk number 13: sessionData
###################################################
sessionInfo()

}
