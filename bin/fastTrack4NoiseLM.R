source("~/Projects/biter_biter_shared/bin/auxggplot2.R");

# TODO set 
test.combineFiles <- function(pat="Brett_2015_RNAseq") {
	source("fastTrack4NoiseLM.R")
	combineFiles(pat)

}


combineFiles <- function(pat="Brett_2015_RNAseq.summary", rep.pair.sub="_(\\d+)$", 
													ofle.tag=paste(pat, "_combined", sep="")) {

	ggplot <- set.ggplot()
	plot.new(); 
	pdf(paste(ofle.tag,".pdf",sep=""),h=10,w=10)  

	files <- dir(pattern=paste(pat, sep=""))
	files <- grep(".summary$", files, value=T)

	# Prep for join
	fnames <- sub("^([^.]*).*", "\\1", files)
	dtaa <- (lapply(files, function(f) readBigTable(f,header=T)))
	dtaa[[1]] <- dtaa[[1]] %>% 
		mutate(Name=gsub(paste(fnames[1], "_X*",sep=""), "", Name))
	dtaa[[2]] <- dtaa[[2]] %>% 
		mutate(Name=gsub(paste(fnames[2], "_X*",sep=""), "", Name))
	dtaa[[3]] <- dtaa[[3]] %>% 
		mutate(Name=gsub(rep.pair.sub, "", Id), Pair=paste("Pair", getRepNumber(Id), sep="")) %>%
		dplyr::select(-Id) %>% spread(Pair, GC_o_ACTG) %>% 
		mutate(GC=(Pair1+Pair2)/2)

	# join list of data frames
	dta <- plyr::join_all(dtaa) %>% mutate(Replicate=getRepNumber(Name), Name=gsub(rep.pair.sub,"",Name))
	dta %>% head

	head(dta)

	fcoef <- lm(Pair2 ~ Pair1, data = dta) %>% coef
	p <- my_labdpoint(dta,aes(x=Pair1, y=Pair2, colour=Name, 
														label=sprintf("%s:%s;%.1f", Name, Replicate, 100*GC)),
										size=(rel(4))) +
		geom_abline(color="gray", linetype=2, size=rel(1), intercept=fcoef[1], slope=fcoef[2]) +
		labs(title="Nucleotide GC / ACGT of read pairs")
	p %>% print 

	fcoef <- lm(FeatureReadsPercent ~ GC, data = dta) %>% coef
	p <- my_labdpoint(dta,aes(x=GC, y=FeatureReadsPercent, colour=Name, 
														label=paste(Name,":",Replicate,";",floor(FeatureReads/1000000),"M", sep=""), shape=Name), 
										size=(rel(4))) +
		geom_abline(color="gray", linetype=2, size=rel(1), intercept=fcoef[1], slope=fcoef[2]) +
		labs(title="Reads count of feature / effective")
	p %>% print 

	fcoef <- lm(ExpressedFeaturesPercent ~ FeatureReadsPercent, data = dta) %>% coef
	p <- my_labdpoint(dta,aes(x=FeatureReadsPercent, y=ExpressedFeaturesPercent, colour=Name, 
														label=paste(Replicate, ExpressedFeatures, sep=":"), shape=Name), 
										size=(rel(4))) +
		geom_abline(color="gray", linetype=2, size=rel(1), intercept=fcoef[1], slope=fcoef[2]) +
		labs(title="ExpressedFeatures after GMM2 fitting")
	p %>% print 

	fcoef <- lm(ExpressedFeaturesPercent ~ GC, data = dta) %>% coef
	p <- my_labdpoint(dta,aes(x=GC, y=ExpressedFeaturesPercent, colour=Name, 
														label=paste(Name,":",Replicate,";",ExpressedFeatures,"\n",floor(FeatureReads/1000000),"M",sep=""), shape=Name), 
							 size=(rel(4))) +
		geom_abline(color="gray", linetype=2, size=rel(1), intercept=fcoef[1], slope=fcoef[2]) +
		labs(title="ExpressedFeatures after GMM2 fitting wrt GC content")
	p %>% print 

	dev.off()
}
