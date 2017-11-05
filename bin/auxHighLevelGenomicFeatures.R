source("~/Projects/biter_biter_shared/bin/bigData.R")

# FIXME doesn't plot; use the data from the example and test plotting
# Does not scale well
test.caller.getSeq <- 
	function(alg="bag", nseq=500, nc = 1, ii=4, train.rate=0.8) {
		library(doMC); registerDoMC(cores = nc);

		#plot.new(); pdf(paste(alg,".pdf",sep=""),h=8,w=8)
		source("~/Projects/biter_biter_shared/bin/auxHighLevelGenomicFeatures.R")
		fle <- "~/PI_HOME/Data/casco/ENCODE/hg19/Synapse_syn6131484/labels/TAF1.train.labels.tsv.gz"

		# read and convert to GRanges 
		#seq.tab <- readBigTable(fle, header=T, nrow=5000000, col.names=c("chr","start","end"))
		# 30K gave segmentation fault for DNAshapeR
		nrow <- 20000
		#nrow <- 2000
		seq.tab <- readBigTable(fle, header=T, nrow=nrow, col.names=c("chr","start","end"))
		#seq.test <- seq.tab[1:(nrow(seq.tab)-1),]
		seq.tab <- seq.tab[(nrow(seq.tab)-nseq): nrow(seq.tab),]
		levels(seq.tab[,ii])
		seq.gr <- convertDataFormat(seq.tab, frmt="GRanges")
		names(seq.gr) <- mcols(seq.gr)$name
		seq.gr

		# TODO parallelize rest
		# get sequences
		library(BSgenome.Hsapiens.UCSC.hg19)
		#seqnms <- levels(seqnames(seq.gr))
		#chr <- unmasked(Hsapiens[[seqnms[1]]])
		#dss <- DNAStringSet(Views(chr, start=start(seq.gr), end=end(seq.gr), names=mcols(seq.gr)$name))
		dss <- getSeq(Hsapiens, names=seq.gr);
		dss

		# Get shapes
		library(Biostrings)
		fle.tag <- paste(alg,"_nseq", nseq, "_", basename(fle),sep="")
		fle.fa <- paste(fle.tag, ".fa", sep="")
		Biostrings::writeXStringSet(dss, file=fle.fa, format="fasta", width=201)
		library(DNAshapeR)
		#featureType <- c("1-mer", "1-shape")
		featureType <- c("1-shape")
		pred <- getShape(fle.fa, "All", parse=T)
		#pred <- getShape(fle.fa) # 
		featureVector <- encodeSeqShape(fle.fa, pred, featureType, normalize=T)
		# TODO variant normalize=F
		#featureVector <- encodeSeqShape(fle.fa, pred, featureType, normalize=F)
		head(featureVector[1:3,1:5])
		rownames(featureVector) <- names(seq.gr)

		writeData(featureVector, paste(fle.tag, ".DNAshape.txt.gz", sep=""), row.names=T)

		# ML
		#prepare something for Tom by plotting DNA shape features :)
		# get metrics from base method implementation
		library(caret)
		library(caretEnsemble)
		folds <- 11
		repeats <- 3
		dataset <- data.frame(Class=seq.tab[,ii], featureVector)
		dataset <- droplevels(subset(dataset, Class != "A"))
		dim(dataset)


		#Split train/test
		train <- runif(nrow(dataset)) <= train.rate

		control <- trainControl(method="repeatedcv", number=folds, repeats=repeats, 
														savePredictions=T, classProbs=T, 
														#index=createMultiFolds(dataset[, "Class"], k=folds, times=repeats))
														index=createMultiFolds(dataset[train, "Class"], k=folds, times=repeats))
		#seed <- 7
		metric <- "Accuracy"

		#set.seed(seed)
		if (alg == "bag")
			mList <- c('treebag', 'rf'); # Bagging algorithms
		if (alg == "boost")
			mList <- c('C5.0', 'gbm') ; # Boosting algorithms
		if (alg == "stack") 
			mList <- c('lda', 'rpart', 'glm', 'knn', 'svmRadial'); # Stacking algorithms
		print(alg)
		models <- caretList(Class~ ., data = dataset[train,], metric = metric,
											 trControl=control, methodList=mList, preProcess=NULL, continue_on_fail=T)

		save(models, file=paste(fle.tag, "_models.rda", sep=""))

		library(caTools); #for colAUC
		preds <- data.frame(sapply(models, function(x){predict(x, dataset[!train,-1], type='prob')[,2]}))
		print(sort(data.frame(colAUC(preds, dataset[!train,1]))))
		#print(models)

		# summarize results
		# debug stack resamples
		#resamples <- lapply(models, function(x) x[["Resample"]])
		#names(resamples) <- names(names)
		#print(unique(resamples))
		#boosting_results <- resamples(list(c5.0=fit.c50, gbm=fit.gbm))
		results <- resamples(models)
		print(summary(results))
		dotplot(results)
		ggsave(paste(fle.tag, "_dotplot.pdf", sep=""))

		# correlation between results
		print(modelCor(results))
		splom(results)
		ggsave(paste(fle.tag, "_splom.pdf", sep=""))

		#dev.off()

		# stack using glm
		stackControl <- control
		#set.seed(seed)
		stack.glm <- caretStack(models, method="glm", metric=metric, trControl=stackControl)
		print(stack.glm)

		# stack using random forest
		#set.seed(seed)
		stack.rf <- caretStack(models, method="rf", metric=metric, trControl=stackControl)
		print(stack.rf)


		# FIXME; did not run 
#		preds$stack.glm <- predict(stack.glm, newdata=dataset[!train,-1])[,2]
#		preds$stack.rf <- predict(stack.rf, newdata=dataset[!train,-1], type='prob')[,2]
#		print(sort(data.frame(colAUC(preds, dataset[!train,1]))))

		print(warnings())
		print(lsos())

		return()
		# ----------------------------------------------------------------------------------

		trainControl <- trainControl(method="cv", number=3, savePredictions=T, classProbs=T) 
		#model <- train(class~ ., data = df, trControl=trainControl, method="parRF", preProcess=NULL) 
		mList <- c('lda', 'rpart', 'glm', 'knn', 'svmRadial') ; #stacking algorithms
		models <- caretList(class~ ., data = droplevels(subset(df, class != "A")), 
											 trControl=trainControl, methodList=mList, preProcess=NULL) 

		results <- resamples(models)
		summary(results)
		dotplot(results)
		# correlation between results
		modelCor(results)
		splom(results)

		#Make a greedy ensemble - currently can only use RMSE
		greedy <- caretEnsemble(models, iter=1000L)
		sort(greedy$weights, decreasing=TRUE)
		greedy$error

		#Make a linear regression ensemble
		linear <- caretStack(models, method='glm', trControl=trainControl(method='cv'))
		linear$error

		# ---------------------
		#Data
		library(mlbench)
		dat <- mlbench.xor(500, 2)
		X <- data.frame(dat$x)
		Y <- factor(ifelse(dat$classes=='1', 'Yes', 'No'))

		#Split train/test
		train <- runif(nrow(X)) <= .66

		library(caretEnsemble)
		#Setup CV Folds
		#returnData=FALSE saves some space
		folds=5
		repeats=1
		myControl <- trainControl(method='cv', number=folds, repeats=repeats, 
															returnResamp='none', classProbs=TRUE,
															returnData=FALSE, savePredictions=TRUE, 
															verboseIter=TRUE, allowParallel=TRUE,
															summaryFunction=twoClassSummary,
															index=createMultiFolds(Y[train], k=folds, times=repeats))
		PP <- c('center', 'scale')

		#Train some models
		model1 <- train(X[train,], Y[train], method='gbm', trControl=myControl,
										                tuneGrid=expand.grid(n.trees=500, interaction.depth=15, shrinkage = 0.01, n.minobsinnode=10 ))
		model2 <- train(X[train,], Y[train], method='blackboost', trControl=myControl)
		model3 <- train(X[train,], Y[train], method='parRF', trControl=myControl)
		model4 <- train(X[train,], Y[train], method='mlpWeightDecay', trControl=myControl, trace=FALSE, preProcess=PP)
		model5 <- train(X[train,], Y[train], method='knn', trControl=myControl, preProcess=PP)
		model6 <- train(X[train,], Y[train], method='earth', trControl=myControl, preProcess=PP)
		model7 <- train(X[train,], Y[train], method='glm', trControl=myControl, preProcess=PP)
		model8 <- train(X[train,], Y[train], method='svmRadial', trControl=myControl, preProcess=PP)
		model9 <- train(X[train,], Y[train], method='gam', trControl=myControl, preProcess=PP)
		model10 <- train(X[train,], Y[train], method='glmnet', trControl=myControl, preProcess=PP)
		algorithmList <- c('lda', 'rpart', 'glm', 'knn', 'svmRadial')
		model10.0 <- caretList(class~., data=df, trControl=myControl, methodList=algorithmList)
		results <- resamples(model10.0)
		summary(results)
		dotplot(results)

		#Make a list of all the models
		all.models <- list(model1, model2, model3, model4, model5, model6, model7, model8, model9, model10)
		names(all.models) <- sapply(all.models, function(x) x$method)
		sort(sapply(all.models, function(x) min(x$results$ROC)))

		#Make a greedy ensemble - currently can only use RMSE
		greedy <- caretEnsemble(caretList(methodList(all.models)), iter=1000L)
		greedy <- caretEnsemble(all.models, iter=1000L)
		sort(greedy$weights, decreasing=TRUE)
		greedy$error

		#Make a linear regression ensemble
		linear <- caretStack(all.models, method='glm', trControl=trainControl(method='cv'))
		linear$error

		#Predict for test set:
		library(caTools)
		preds <- data.frame(sapply(all.models, function(x){predict(x, X[!train,], type='prob')[,2]}))
		preds$ENS_greedy <- predict(greedy, newdata=X[!train,])[,2]
		preds$ENS_linear <- predict(linear, newdata=X[!train,], type='prob')[,2]
		sort(data.frame(colAUC(preds, Y[!train])))


	}

caller.getSeq <-
	function(dta, fle=NULL, gnm="mm10", ...) {

		if (!is.null(fle))
			dta <- readBigTable(fle, frmt="GRanges", ...)

		dta

#		library(BSgenome.Hsapiens.UCSC.hg19)
#		genSeq <- Hsapiens
#		chr1 <- unmasked(genSeq$chr1)
#		r1 <- IRanges(start=c(800,850), width=c(200,200))
#		myseqs <- Views(chr1, start=start(r1), end=end(r1))

	}

# http://stackoverflow.com/questions/1358003/tricks-to-manage-the-available-memory-in-an-r-session
.ls.objects <- function (pos = 1, pattern, order.by,
												 decreasing=FALSE, head=FALSE, n=5) {
	napply <- function(names, fn) sapply(names, function(x)
																			 fn(get(x, pos = pos)))
	names <- ls(pos = pos, pattern = pattern)
	obj.class <- napply(names, function(x) as.character(class(x))[1])
	obj.mode <- napply(names, mode)
	obj.type <- ifelse(is.na(obj.class), obj.mode, obj.class)
	obj.size <- napply(names, object.size)
	obj.dim <- t(napply(names, function(x)
											as.numeric(dim(x))[1:2]))
	vec <- is.na(obj.dim)[, 1] & (obj.type != "function")
	obj.dim[vec, 1] <- napply(names, length)[vec]
	out <- data.frame(obj.type, obj.size, obj.dim)
	names(out) <- c("Type", "Size", "Rows", "Columns")
	if (!missing(order.by))
		out <- out[order(out[[order.by]], decreasing=decreasing), ]
	if (head)
		out <- head(out, n)
	out
}
# shorthand
lsos <- function(..., n=10) {
	.ls.objects(..., order.by="Size", decreasing=TRUE, head=TRUE, n=n)
}
