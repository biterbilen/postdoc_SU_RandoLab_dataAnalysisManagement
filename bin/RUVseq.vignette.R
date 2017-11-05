> library(zebrafishRNASeq)
> data(zfGenes)
> head(zfGenes)
                   Ctl1 Ctl3 Ctl5 Trt9 Trt11 Trt13
ENSDARG00000000001  304  129  339  102    16   617
ENSDARG00000000002  605  637  406   82   230  1245
ENSDARG00000000018  391  235  217  554   451   565
ENSDARG00000000019 2979 4729 7002 7309  9395  3349
ENSDARG00000000068   89  356   41  149    45    44
ENSDARG00000000069  312  184  844  269   513   243
> tail(zfGenes)
            Ctl1 Ctl3 Ctl5  Trt9 Trt11 Trt13
ERCC-00163   204   59  183   152   104    59
ERCC-00164     6    1   74    11   206    21
ERCC-00165   140  119   93   331    52    38
ERCC-00168     0    0    0     0     2     0
ERCC-00170   216  145  111   456   196   552
ERCC-00171 12869 6682 7675 47488 24322 26112
> 
> ## ----filter--------------------------------------------------------------
> filter <- apply(zfGenes, 1, function(x) length(x[x>5])>=2)
> filtered <- zfGenes[filter,]
> genes <- rownames(filtered)[grep("^ENS", rownames(filtered))]
> spikes <- rownames(filtered)[grep("^ERCC", rownames(filtered))]
> 
> ## ----store <- data----------------------------------------------------------
> x <- as.factor(rep(c("Ctl", "Trt"), each=3))
> set <- newSeqExpressionSet(as.matrix(filtered),
														 +                            phenoData = data.frame(x, row.names=colnames(filtered)))
> set
SeqExpressionSet (storageMode: lockedEnvironment)
assayData: 20865 features, 6 samples 
  element names: counts, normalizedCounts, offset 
protocolData: none
phenoData
  sampleNames: Ctl1 Ctl3 ... Trt13 (6 total)
  varLabels: x
	  varMetadata: labelDescription
	featureData: none
	experimentData: use 'experimentData(object)'
	Annotation:  
	> 
	> ## ----rle, fig.cap="No normalization.",fig.subcap=c("RLE plot","PCA plot")----
	> library(RColorBrewer)
	> colors <- brewer.pal(3, "Set2")
	> plotRLE(set, outline=FALSE, ylim=c(-4, 4), col=colors[x])
	> plotPCA(set, col=colors[x], cex=1.2)
	> 
	> ## ----uq, fig.cap="Upper-quartile normalization.", fig.subcap=c("RLE plot","PCA plot")----
	> set <- betweenLaneNormalization(set, which="upper")
	> plotRLE(set, outline=FALSE, ylim=c(-4, 4), col=colors[x])
	> plotPCA(set, col=colors[x], cex=1.2)
	> 
	> ## ----ruv <- spikes, fig.cap="RUVg normalization based on spike-in controls.", fig.subcap=c("RLE plot","PCA plot")----
	> set1 <- RUVg(set, spikes, k=1)
	> pData(set1)
	        x         W <- 1
	Ctl1  Ctl -0.04539413
	Ctl3  Ctl  0.50347642
	Ctl5  Ctl  0.40575319
	Trt9  Trt -0.30773479
	Trt11 Trt -0.68455406
	Trt13 Trt  0.12845337
	> design <- model.matrix(~x, data=pData(set))
	> y <- DGEList(counts=counts(set), group=x)
	> y <- calcNormFactors(y, method="upperquartile")
	> y <- estimateGLMCommonDisp(y, design)
	> y <- estimateGLMTagwiseDisp(y, design)
	> 
	> fit <- glmFit(y, design)
	> lrt <- glmLRT(fit, coef=2)
	> 
	> top <- topTags(lrt, n=nrow(set))$table
	> empirical <- rownames(set)[which(!(rownames(set) %in% rownames(top)[1:5000]))]
	> 
	> ## ----emp <- ruvg, fig.cap="RUVg normalization based on empirical controls.", fig.subcap=c("RLE plot","PCA plot")----
	> set2 <- RUVg(set, empirical, k=1)
	> pData(set2)
	        x         W <- 1
	Ctl1  Ctl -0.10879677
	Ctl3  Ctl  0.23066424
	Ctl5  Ctl  0.19926266
	Trt9  Trt  0.07672121
	Trt11 Trt -0.83540924
	Trt13 Trt  0.43755790
	> differences <- makeGroups(x)
	Error: could not find function "makeGroups"
	> differences
	Error: object 'differences' not found
	> makeGroups <- function(xs) {
		+ xs <- factor(xs)
		+ groups <- matrix(-1, nrow = length(levels(xs)), ncol = max(table(xs)))
		+ for (i in 1:length(levels(xs))) {
			+ idxs <- which(xs == levels(xs)[i])
			+ groups[i,1:length(idxs)] <- idxs
			+ }
		+ groups
		+ }
	> differences <- makeGroups(x)
	> differences
	     [,1] [,2] [,3]
	[1,]    1    2    3
	[2,]    4    5    6
	> set3 <- RUVs(set, genes, k=1, differences)
	> pData(set3)
	        x        W <- 1
	Ctl1  Ctl  0.1860825
	Ctl3  Ctl  0.4917394
	Ctl5  Ctl  0.4926805
	Trt9  Trt  0.4279274
	Trt11 Trt -0.5784460
	Trt13 Trt  0.7289145
	> design <- model.matrix(~x, data=pData(set))
	> y <- DGEList(counts=counts(set), group=x)
	> y <- calcNormFactors(y, method="upperquartile")
	> y <- estimateGLMCommonDisp(y, design)
	> y <- estimateGLMTagwiseDisp(y, design)
	> fit <- glmFit(y, design)
	> res <- residuals(fit, type="deviance")
	set4 <- RUVr(set, genes, k=1, res)
	> set4 <- RUVr(set, genes, k=1, res)
	> pData(set4)
	        x        W <- 1
	Ctl1  Ctl -0.2392330
	Ctl3  Ctl  0.1430670
	Ctl5  Ctl  0.1449671
	Trt9  Trt  0.1018808
	Trt11 Trt -0.7384969
	Trt13 Trt  0.5878151
	> q()

