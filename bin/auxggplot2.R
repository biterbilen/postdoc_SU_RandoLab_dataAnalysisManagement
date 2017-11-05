source("~/Projects/biter_biter_shared/bin/bigData.R");
#library(multidplyr) # TODO integrate later
library(magrittr)
library(dplyr)
library(tidyr)
library(parallel)
library(gridExtra) #grid.arrange
library(ggthemes)
library(ggplot2)
library(GGally)
library(ggrepel)

# TODO LIST:
# replace redundant code? 
# test long data
test.all <- 
	function() {
		source("auxggplot2.R")
		test.my_wrapp()
		test.my_bar_hist()
		test.my_ecdf()
		test.my_boxplot() 
		test.my_dens()
		test.my_labdpoint()
		test.my_bin()
		test.my_tile()
		test.my_seg()
		test.get.corrheat_seg_boxplot()
		test.get.splom() 
		test.get.svdplot() 
		test.setformula() 
		test.get.plot() 
	}

set.ggplot <- 
	function(...) {
		N <- 20 # default gdocs color pal size is 20  
		cols <- ggthemes_data$gdocs[1:N]
		#cols <- unname(ggthemes_data$colorblind); cols[1] <- "#999999" #make black gray
		ggplot <- 
			function(...) 
				ggplot2::ggplot(...) + 
					theme_bw() +
					theme(text = element_text(size=rel(4), angle = 0),
						line = element_line(size=rel(1)),
						rect = element_rect(size=rel(1)),
						plot.title = element_text(size=rel(4)),
						#legend.title = element_blank(),
						legend.direction = "horizontal", legend.box = "vertical",
						legend.text = element_text(size=rel(2)),
						legend.position = "top", 
						legend.key = element_rect(colour=NA),
						legend.key.size = unit(0.3,"cm"), # cannot be rel
						strip.text = element_text(size=rel(1)),
						strip.background = element_blank()) +
					scale_color_manual(values=c(cols,cols)) + 
					scale_fill_manual(values=c(cols,cols)) + 
					scale_shape_manual(values=as.factor(c(1:100))) #+
		unlockBinding("ggplot",parent.env(asNamespace("GGally")))
		assign("ggplot",ggplot,parent.env(asNamespace("GGally")))
		ggplot
	}

# adapted from Ggally package
is_discrete <- function(data, mapping, val = "y") {
		yData <- eval(mapping[[val]], data)
		is.factor(yData) || is.character(yData)
	}

mapping_swap_x_y <- 
	function(mapping, x="x", y="y") {
		tmp <- mapping[[x]]
		mapping[[x]] <- mapping[[y]]
		mapping[[y]] <- tmp
		mapping
	}

test.my_wrapp <-
	function() {
		lst <- "blank"
		lst <- list(discrete="na", continuous=list(f="my_dens", param=list(size=rel(2))))
		lst <- list(continuous=list(f="my_bin",   param=list(size=rel(2))),
								discrete=list(f="my_boxplot", param=list(size=rel(2))))
		source("auxggplot2.R")
		my_wrapp(lst)
	}

my_wrapp <-
	function(lst, verbose=T, ...) {
		if (length(lst) == 0)
			res = lst
		else
			res = list()
		if (verbose)
			print(lst)
		for (i in names(lst)) {
			if (is.atomic(lst[[i]])) {
				res[[i]] = lst[[i]]
			} else {
				res[[i]] = wrapp(get(lst[[i]]$f), lst[[i]]$param)
			}
		}
		res
	}

test.my_bar_hist <-
	function() {
		data <- iris
		names(data) <- c("sepal_1", "sepal_2", "petal_1", "petal_2", "Species")
		data$Group <- factor(rnorm(nrow(data)) > 0)
		head(data)
		source("auxggplot2.R")
		my_bar_hist(data, columns=1:4, size=rel(2)) # C 
		my_bar_hist(data, columns=4:5, size=rel(2)) # D 
		my_bar_hist(data, columns=1:4, facet="Species", size=rel(2)) # C w facet
		# discrete -actually geom_bar 
		my_bar_hist(data, long=F, mapping=aes(x=Group, y=Species, color=Species), size=rel(2)) # D C w color
		my_bar_hist(data, long=F, mapping=aes(x=Group, color=Species), size=rel(2))

		# stat = identity, OK
		data <- data.frame(sepal=iris[1:3,1], label=c("a","b","c"))
		source("auxggplot2.R")
		my_bar_hist(data, mapping=aes(x=label, y=sepal, color=sepal>5), long=F, stat="identity", size=rel(2))
		my_bar_hist(data, mapping=aes(x=sepal, y=label, color=sepal>5, label=label), long=F, stat="identity", size=rel(2)) # for GSEA
		# stat = count wo/y, OK
		data <- data.frame(sepal=iris[1:4,1], label=c("a","b","c","c"))
		source("auxggplot2.R")
		my_bar_hist(data, mapping=aes(x=label, color=sepal>5), long=F, stat="count", size=rel(2))
	
	}

# discrete mappings in diag and lower ggpairs
my_bar_hist <- 
	function(data, mapping=NULL, columns=1:ncol(data), rep.sub='_(\\d+)$', filter=NA, 
					 long=T,
					 facet=NA, title="", nrow=NULL, ncol=NULL, 
					 stat="identity", fill=NA, verbose=T, ...) {

		if (long) { 
			data <- my_transform(data, columns=columns, rep.sub=rep.sub, filter=filter, verbose=verbose)
			mapping <- structure(mergeSetwDefault(mapping, aes(x=Value, color=Name, group=Key)), class="uneval")
		}

		if (verbose) {
			print("dim data:")
			print(dim(data))
			print(head(data))
			print(mapping)
		}

		discretex <- is_discrete(data, mapping, val="x")
		discretey <- is_discrete(data, mapping, val="y")

		if (discretey)
			mapping <- mapping_swap_x_y(mapping, "x", "y")

		if (verbose) {
			print("updated:")
			print(mapping)
		}
		ggplot <- set.ggplot()
		p <- ggplot(data = data, mapping = mapping)

		p <- p + geom_histogram(stat=stat, alpha=0.5, fill=fill, ...)

		# facets same
		if (! is.na(facet)) 
			p <- p + facet_wrap(as.formula(paste0("~", facet)), nrow=nrow, ncol=ncol)

		if (discretey) {
			p <- p + coord_flip()
			if ("label" %in% names(mapping))
				p <- p + geom_text(aes(y=0.01), hjust="inward", vjust="inward", show.legend=F)
		}

		p <- p + labs(title=title) #same
		p 
	}

test.my_ecdf <-
	function() {
		data <- iris
		names(data) <- c("sepal_1", "sepal_2", "petal_1", "petal_2", "Species")
		head(data)
		source("auxggplot2.R")
		my_ecdf(data, columns=1:4, size=rel(2), title="combine=F") 
		my_ecdf(data, columns=1:4, combine=T, size=rel(2), title="combine=T") 
		my_ecdf(data, columns=1:4, facet="Species", size=rel(2))
		my_ecdf(data, columns=1:4, combine=T, facet="Species", size=rel(2))
	}

my_ecdf <- 
	function(data, mapping=NULL, columns=1:ncol(data), rep.sub='_(\\d+)$', filter=NA, 
					 long=T, combine=F, round=F,
					 facet=NA, title="", nrow=NULL, ncol=NULL, verbose=T, ...) {

		if (long) { 
			data <- my_transform(data, columns=columns, rep.sub=rep.sub, filter=filter, combine=combine, round=round, verbose=verbose)
			mapping <- structure(mergeSetwDefault(mapping, aes(x=Value, color=Name, group=Key)), class="uneval")
			if (combine) 
				mapping2 <- structure(mergeSetwDefault(aes(color=Key), mapping), class="uneval")
		}

		if (verbose) {
			print(mapping)
		}

		ggplot <- set.ggplot()
		p <- ggplot(data = data, mapping = mapping)

		if (combine) {
			p <- p +
				stat_ecdf(data=filter_(data, '!is.na(Replicate)'), alpha=0.3, size=rel(1), show.legend=F) + 
				stat_ecdf(data=filter_(data, 'is.na(Replicate)'), mapping=mapping2, ...)
		} else {
			p <- p + stat_ecdf(...)
		}

		# facets #same
		if (! is.na(facet)) 
			p <- p + facet_wrap(as.formula(paste0("~", facet)), nrow=nrow, ncol=ncol)

		p <- p + labs(title=title, y="Empirical CDF") #same

		p 
	}

test.my_boxplot <- 
	function() {
		source("auxggplot2.R")
		data <- iris
		names(data) <- c("sepal_1", "sepal_2", "petal_1", "petal_2", "Species")
		head(data)
		my_boxplot(data, columns=1:4, panel="dot", facet="Species") 
		my_boxplot(data, columns=1:4, wjitter=T, panel="box") 
		my_boxplot(data, columns=1:4, wjitter=T, panel="violin") 
		my_boxplot(data, columns=1:4, wjitter=T, panel="violin", facet="Species") 

#		data <- dta.sum %>% data.frame; names(data)[1] <- "Spike"
#		head(data)
#		my_boxplot(data, columns=2:8, filter='Spike==T', mapping=aes(x=Name,label=Replicate), panel="dot") 
#		my_boxplot(data, columns=2:8, filter='Spike==T', mapping=aes(x=Name,label=Replicate), panel="box", wjitter=T) 

		# dotplots of the cor coefficients
		data <- iris
		names(data) <- c("sepal_1", "sepal_2", "petal_1", "petal_2", "Species")
		data <- cor(data[,1:4])
		my_boxplot(data, wjitter=T, panel="box", mapping=aes(x=Name))
		my_boxplot(data, wjitter=T, panel="box", mapping=aes(x=Name), filter="NameR==Name & ReplicateR>Replicate")

		# horizontal
		data <- iris
		names(data) <- c("sepal_1", "sepal_2", "petal_1", "petal_2", "Species")
		source("auxggplot2.R")
		my_boxplot(data, columns=1:5, mapping=aes(y=Species, x=petal_1), long=F, wjitter=T, panel="box")# + coord_flip()

	}

my_boxplot <- 
	function(data, mapping=NULL, columns=1:ncol(data), rep.sub='_(\\d+)$', filter=NA,
					 long=T,
					 facet=NA, panel="violin", wjitter=F, title="", nrow=NULL, ncol=NULL, verbose=T, ...) {
		if (long) { 
			data <- my_transform(data, columns=columns, rep.sub=rep.sub, filter=filter, verbose=verbose)
			mapping <- structure(mergeSetwDefault(mapping, aes(x=Key, y=Value, color=Name)), class="uneval")
		}

		flip <- is_discrete(data, mapping, val="y")
		if (flip)
			mapping <- mapping_swap_x_y(mapping, "x", "y")

		if (verbose) 
			print(mapping)

		ggplot <- set.ggplot()
		p <- ggplot(data = data, mapping = mapping)

		if (panel == "violin") {
			p <- p + geom_violin(draw_quantiles=c(0.05,0.25,0.5,0.75,0.95), scale="count", ...)
		} else if (panel == "dot") {
			p <- p + geom_jitter(shape="+", stroke=rel(2), size=rel(6), width=0.1, height=0, ...)
			if ("label" %in% names(mapping))
				p <- p + geom_text_repel(show.legend=F)
		} else if (panel == "box") {
			p <- p + geom_boxplot(outlier.shape=NA, ...)
		}

		if (panel != "dot" && wjitter) 
			p <- p + geom_jitter(shape="-", color="black", stroke=rel(2), size=rel(6), width=0.1, height=0, alpha=0.5)

		# facets
		if (! is.na(facet)) 
			p <- p + facet_wrap(as.formula(paste0("~", facet)), nrow=nrow, ncol=ncol)

		# labels
		if (mapping$x == "Key" && "Replicate" %in% names(data))
			p <- p +	scale_x_discrete(labels=getRepNumber) + labs(x="Replicate")

		p <- p + labs(title=title)

		if (flip) 
			p <- p + coord_flip()

		p 
	}

test.my_dens <-
	function() {
		source("auxggplot2.R")
		data <- iris
		names(data) <- c("sepal_1", "sepal_2", "petal_1", "petal_2", "Species")
		head(data)
		my_dens(data, columns=1:4, size=rel(2)) 
		my_dens(data, columns=1:4, facet="Species", size=rel(2))
	}

my_dens <- 
	function(data, mapping=NULL, columns=1:ncol(data), rep.sub='_(\\d+)$', filter=NA,
					 long=T,
					 facet=NA, title="", nrow=NULL, ncol=NULL, verbose=T, ...) {
		if (long) { 
			data <- my_transform(data, columns=columns, rep.sub=rep.sub, filter=filter, verbose=verbose)
			mapping <- structure(mergeSetwDefault(mapping, aes(x=Value, color=Name)), class="uneval")
		}

		if (verbose) 
			print(mapping)

		ggplot <- set.ggplot()
		p <- ggplot(data = data, mapping = mapping) + 
			geom_density(...) + labs(y="Density")
								 
		if (! is.na(facet)) 
			p <- p + facet_wrap(as.formula(paste0("~", facet)), nrow=nrow, ncol=ncol)

		p <- p + labs(title=title)

		p 
	}	

test.my_text <-
	function() {
		data <- iris
		names(data) <- c("sepal_1", "sepal_2", "petal_1", "petal_2", "Species")
		head(data)
		source("auxggplot2.R")
		my_text(data, columns=1:2, mapping=aes(x=sepal_1))
		my_text(data, label="MY_LABEL")
	}


my_text <- 
	function(data, mapping=NULL, columns=1:ncol(data), rep.sub='_(\\d+)$', filter=NA, textRepelFilter=NA,
					 column.labels=if(length(columns)<5) names(data)[columns] else unname(abbreviate(names(data)[columns])),
					 long=F,
					 x=0, y=0, label=NULL,
					 facet=NA, title="", nrow=NULL, ncol=NULL, verbose=T, ...) {

		if (long) { 
			data <- my_transform(data, columns=columns, rep.sub=rep.sub, filter=filter, verbose=verbose)
			mapping <- structure(mergeSetwDefault(mapping, aes(x=Value, y=Value, color=Name, label=Value)), 
													 class="uneval")
		}

		if (is.null(label))
			label <- column.labels[names(data) == as.character(mapping$x)]
		data <- data.frame(x=x, y=x, label=label)
		mapping <- aes(x,y,label=label)


		if (verbose) {
			print(mapping)
			print(data)
		}

		ggplot <- set.ggplot()
		p <- ggplot(data = data, mapping = mapping) + 
			geom_text(size=rel(4), ...) + theme_void()

		#p <- p + geom_text(aes(x = -Inf, y = Inf), hjust="inward", vjust="inward", show.legend=F)

		if (! is.na(facet)) 
			p <- p + facet_wrap(as.formula(paste0("~", facet)), nrow=nrow, ncol=ncol)

		p <- p + labs(title=title)

		p
	}

test.my_labdpoint <-
	function() {
		data <- iris
		names(data) <- c("sepal_1", "sepal_2", "petal_1", "petal_2", "Species")
		head(data)
		source("auxggplot2.R")
		my_labdpoint(data, mapping=aes(x=sepal_1, y=sepal_2, label=petal_2, color=Species == 'setosa'),
								 textRepelFilter="Species == 'setosa'",
								 size=rel(2)) 
	}

my_labdpoint <- 
	function(data, mapping=NULL, columns=1:ncol(data), rep.sub='_(\\d+)$', filter=NA, textRepelFilter=NA,
					 long=F,
					 facet=NA, title="", nrow=NULL, ncol=NULL, verbose=T, ...) {
		if (long) { 
			data <- my_transform(data, columns=columns, rep.sub=rep.sub, filter=filter, verbose=verbose)
			mapping <- structure(mergeSetwDefault(mapping, aes(x=Value, y=Value, color=Name, label=Value)), 
													 class="uneval")
		}

		if (verbose) 
			print(mapping)

		ggplot <- set.ggplot()
		p <- ggplot(data = data, mapping = mapping) + 
			geom_point(stroke=rel(2), ...)

		if ("label" %in% names(mapping)) {
			if (is.na(textRepelFilter)) {
				p <- p + geom_text_repel(show.legend=F, color="black", ...)
			} else {
				p <- p + geom_text_repel(data=filter_(data, textRepelFilter), show.legend=F, color="black", ...)
			}
		}
								 
		if (! is.na(facet)) 
			p <- p + facet_wrap(as.formula(paste0("~", facet)), nrow=nrow, ncol=ncol)

		p <- p + labs(title=title)

		p
	}

test.my_bin <-
	function() {
		data <- iris
		names(data) <- c("sepal_1", "sepal_2", "petal_1", "petal_2", "Species")
		head(data)
		source("auxggplot2.R")
		my_bin(data, mapping=aes(x=sepal_1, y=sepal_2, label=petal_1, color= Species == 'setosa'),
					 textRepelFilter="Species == 'setosa'", smooth=F, bins=150,
					 size=rel(2)) 
	}


my_bin <- 
	function(data, mapping=NULL, columns=1:ncol(data), rep.sub='_(\\d+)$', filter=NA, textRepelFilter=NA,
					 long=F,
					 smooth=T, se=F, method="loess", bins=15, slope=1, intercept=0,
					 facet=NA, title="", nrow=NULL, ncol=NULL, verbose=T, ...) {
		if (long) { 
			data <- my_transform(data, columns=columns, rep.sub=rep.sub, filter=filter, verbose=verbose)
			mapping <- structure(mergeSetwDefault(mapping, aes(x=Value, y=Value, color=Name, label=Value)), 
													 class="uneval")
		}

		if (verbose) 
			print(mapping)

		ggplot <- set.ggplot()
		p <- ggplot(data = data, mapping = mapping) +
			#geom_density_2d() +
			#geom_density2d() +
			geom_bin2d(bins=c(bins,bins), size=rel(1), ...) +
			scale_fill_gradientn(colours=c("darkgray", "black")) +
			geom_abline(slope=slope, intercept=intercept, color="darkgray", linetype=2, size=rel(1))

		if (smooth)
			p <- p + geom_smooth(method=method, se=se, ...)
		
		# TODO how to replace redundant code? here/ my_labdpoint?
		# p <- p + my_labdpoint
		if (is.na(textRepelFilter)) {
			if ("label" %in% names(mapping)) {
				p <- p + 
					geom_point(stroke=rel(2), alpha=0.5, ...) +
					geom_text_repel(show.legend=F, color="black", ...)
			}
		} else {
			if ("label" %in% names(mapping)) {
				p <- p + 
					geom_point(data=filter_(data, textRepelFilter), stroke=rel(2), alpha=0.5, ...) +
					geom_text_repel(data=filter_(data, textRepelFilter), show.legend=F, color="black", ...)
			}
		}
								 
		if (! is.na(facet)) 
			p <- p + facet_wrap(as.formula(paste0("~", facet)), nrow=nrow, ncol=ncol)

		p <- p + labs(title=title)

		p 
	}

test.my_tile <-
	function() {
		data <- iris
		names(data) <- c("sepal_1", "sepal_2", "petal_1", "petal_2", "Species")
		data <- cor(data[,1:4])
		head(data)
		source("auxggplot2.R")
		my_tile(data, size=rel(2)) 
		
	}

my_tile <- 
	function(data, mapping=NULL, columns=1:ncol(data), rep.sub='_(\\d+)$', filter=NA, 
					 long=T,
					 facet=NA, title="", nrow=NULL, ncol=NULL, verbose=T, ...) {

		if (long) { 
			data <- my_transform(data, columns=columns, rep.sub=rep.sub, filter=filter, verbose=verbose)
			mapping <- structure(mergeSetwDefault(mapping, aes(x=KeyR, y=Key, fill=Value)), class="uneval")
		}

		if (verbose) {
			print(mapping)
		}

		ggplot <- set.ggplot()

		p <- ggplot(data = data, mapping = mapping) +
			geom_tile(...) + 
			scale_fill_gradientn(colours=c("gray30", "white", "gold")) #+
			#scale_fill_gradientn(colours=c("darkblue", "white", "darkred")) #+

		# facets #same
		if (! is.na(facet)) 
			p <- p + facet_wrap(as.formula(paste0("~", facet)), nrow=nrow, ncol=ncol)

		p <- p + labs(title=title) #same

		p
	}

test.my_seg <-
	function() {
		data <- iris
		names(data) <- c("sepal_1", "sepal_2", "petal_1", "petal_2", "Species")
		data <- cor(data[,1:4])
		head(data)
		source("auxggplot2.R")
		my_seg(data, size=rel(2)) 
		
	}

my_seg <- 
	function(data, mapping=NULL, columns=1:ncol(data), rep.sub='_(\\d+)$', filter=NA, 
					 long=F,
					 dend=T, linkage.method="complete",
					 facet=NA, title="", ylab="", nrow=NULL, ncol=NULL, verbose=T, ...) {

		if (long) { 
			data <- my_transform(data, columns=columns, rep.sub=rep.sub, filter=filter, verbose=verbose)
			mapping <- structure(mergeSetwDefault(mapping, aes(x=KeyR, y=Key, fill=Value)), class="uneval")
		}

		# could be generalized for other trees -begin
		if (dend) {
			corr <- data %>% cor(method="pearson", use="pairwise.complete.obs")
			data <- (1 - corr) %>% as.dist %>% hclust(method=linkage.method) %>% as.dendrogram		
			title <- paste(title, ";", linkage.method, " linkage",sep="")
			ylab <- paste(ylab, ";", "Dissimilarity\n(1 - Pearson's R)", sep="")
		}

		library(ggdendro) # for dendro_data
		data %<>% dendro_data
		data$labels %<>%
			tidyr::extract(label, c("Name", "Replicate"), paste0('(.+)', rep.sub), remove=F)
		mapping  <- aes(x=x, y=y, xend=xend, yend=yend)
		mapping2 <- aes(x=x, y=y, color=Name, label=Replicate, angle=0)

		if (verbose) {
			print(data)
			print(mapping)
		}

		ggplot <- set.ggplot()

		# trick for aligning the dendrogram is to add xlim
		# TODO does not work for large trees; try putting units in xlim
		p <- ggplot() +
			geom_segment(data=data$segments, mapping=mapping, ...) + 
			geom_text(data=data$labels, mapping=mapping2, vjust=0.49, size=rel(6), ...) + 
			xlim(0.99, nrow(data$labels) + 0.01) #+ 
#			ylim(range(data$segment$y)) +
#			theme(axis.text.x=element_blank())

		# facets #same
		if (! is.na(facet)) 
			p <- p + facet_wrap(as.formula(paste0("~", facet)), nrow=nrow, ncol=ncol)

		p <- p + labs(title=title, y=ylab, x="") #same

		p

	}


test.get.corrheat_seg_boxplot <-
	function() {
		data <- iris[,1:4]
		data <- bind_cols(data,data,data,data)
		data2 <- rnorm(nrow(data)*ncol(data)) %>% matrix(nrow=nrow(data), ncol=ncol(data)) %>% data.frame
		data <- data + data2
		names(data) <- c("sepal_1", "sepal_2", "petal_1", "petal_2", "sepal_11", "sepal_22", "petal_11", "petal_22", "sepal_111", "sepal_222", "petal_111", "petal_222", "sepal_1111", "sepal_2222", "petal_1111", "petal_2222")
		head(data)
		source("auxggplot2.R")
		get.corrheat_seg_boxplot(data) 
	}	

get.corrheat_seg_boxplot <- 
	function(data, mapping=NULL, columns=1:ncol(data), rep.sub='_(\\d+)$', filter=NA,
					 column.labels=if(length(columns)<5) names(data)[columns] else unname(abbreviate(names(data)[columns])),
					 long=F,
					 linkage.method="complete", nodendrogram=F, noheatmap=F, noboxplot=F,
					 title="", verbose=T, ofle.tag=NA, ...) {

		if (long) { 
			data <- my_transform(data, columns=columns, rep.sub=rep.sub, filter=filter, verbose=verbose)
			mapping <- structure(mergeSetwDefault(mapping, aes(x=Value, y=Value, color=Name)), class="uneval")
		}

		# data
		corr <- dplyr::select(data, columns) %>% cor(method="pearson", use="pairwise.complete.obs")
		dend <- (1-corr) %>% as.dist %>% hclust(method=linkage.method) %>% as.dendrogram		
		ord <- dend %>% order.dendrogram

		if (verbose) {
			print(mapping)
			print(corr)
		}

		# plot data
		if (!is.na(ofle.tag)) {
#			plot.new(); 
			pdf(paste(ofle.tag, ".pdf", sep=""), h = 8, w = 8);
		}

		title <- paste(title, ";", linkage.method, " linkage",sep="")
		ylab <- paste("Dissimilarity\n(1 - Pearson's R)", sep="")

		ggplot <- set.ggplot()

		b <- my_boxplot(corr, mapping=aes(x=Name), filter="NameR==Name & ReplicateR>Replicate",
										wjitter=T, panel="box", title=title, ...)
		m <- my_tile(corr[ord,ord], ...)
		s <- my_seg(dend, dend=F, size=rel(2), title=title, ylab=ylab, ...) 

		if (!noboxplot)
			print(b + labs(y="Within-group similarity\n(Pearson's R)"))

		if (!nodendrogram)
			print(s + theme(legend.position="top", panel.border=element_blank(), axis.text.x=element_blank(), panel.grid.minor=element_blank()) + labs(x="Replicate"))

		if (!noheatmap) 
			grid.arrange(s + theme_void() + theme(legend.position="none", strip.placement="outside"), 
									 grab_legend(m + theme(legend.direction="vertical")),
									 m + theme(axis.title=element_blank(), axis.text=element_blank(), legend.position="none"), 
									 grab_legend(s + theme(legend.direction="vertical")), 
									 ncol=2, nrow=2, widths = c(0.75, 0.25), heights = c(0.25, 0.75)) %>% print

		if (!is.na(ofle.tag))
			dev.off()
	}

test.get.splom <- 
	function() {
		data <- iris
		names(data) <- c("sepal_1", "sepal_2", "petal_1", "petal_2", "Species")
		data$Group <- factor(rnorm(nrow(data)) > 0)
		data$Group2 <- factor(data$petal_2  > 0.2)
		#data$Group <- rnorm(nrow(data)) > 0
		head(data)
		source("auxggplot2.R")
		# default no group
		columns=c(4:5)
		columns=c(4:6)
		columns=c(5,6,4)
		get.splom(data, columns=columns, detailed=F, title="title") 
		# default w/ group
		columns=c(3:6)
		columns=c(5,6,3,4)
		get.splom(data, columns=columns, mapping=aes_string(color="Species"))
		# default w/group 
		source("auxggplot2.R")
		columns=c(3:6)
		columns=c(5,6,3,4)
		columns=c(6,5)
		columns=c(1:6)
plot.new(); 
pdf(paste("seg",".pdf",sep=""), h=12, w=12)
		get.splom(data, columns=columns, mapping=aes_string(color="Group2"))
dev.off()

		# labeled points -in svd
		get.splom(data[96:104,], columns=2:4, mapping=aes_string(color="Species", shape="Group", label="Group"), 
						  lower=list(continuous=list(f="my_labdpoint", param=list(size=rel(4)))),
						  diag=list(continuous=list(f="my_text", param=list(label="", size=rel(4)))),
							detailed=F, title="title")
	}

# lower triangualar scatterplot-like & upper triangular legend & diagonal density or text
get.splom <- 
	function(data, mapping=NULL, columns=1:ncol(data), rep.sub='_(\\d+)$', filter=NA,
					 column.labels=if(length(columns)<5) names(data)[columns] else unname(abbreviate(names(data)[columns])),
					 long=F,
					 title="", diag=NULL, lower=NULL, detailed=F, verbose=T, ofle.tag=NA, ...) {
		# set defaults
		upper = "blank"
		lower = my_wrapp(mergeSetwDefault(lower, 
																			list(continuous=list(f="my_bin",   param=list(size=rel(2))),
																	         discrete=list(f="my_bar_hist", param=list(long=F, size=rel(3))),
																	         combo=list(f="my_boxplot", param=list(long=F, panel="violin", size=rel(2)))
																					 )))
		diag  = my_wrapp(mergeSetwDefault(diag,
																			list(continuous=list(f="my_dens", param=list(long=F, size=rel(2))),
																	         discrete=list(f="my_bar_hist", param=list(long=F, size=rel(2)))
																					 )))
		ggplot <- set.ggplot()
		if ( !is.na(ofle.tag)) {
#			plot.new();
			pdf(paste(ofle.tag, ".pdf", sep=""), h = 8, w = 8);
		}

		if (verbose) {
			print("HEAD data") 
			print(head(data))
			print("MAPPING:", str(mapping))
			print("DIAG:", str(diag))
			print("LOWER:", str(lower))
		}

		# lower and diagonal
		p <- ggpairs(data, mapping=mapping, columns=columns, columnLabels=column.labels,
								 upper=upper, diag=diag, lower=lower, title=title) #, axisLabels = c("show", "show", "show"))

		# upper legend
		for (j in p$ncol:2)  {
			for (i in 1:(j-1)) {
				p[i,j] <- grab_legend(p[j,i] + 
															guides(
																		 colour = guide_legend(order = 1),
																		 shape = guide_legend(order = 1),
																		 fill = guide_legend(order = 2),
																		 size = guide_legend(order = 3)
																		 ) +
															theme(legend.direction="vertical", legend.box="horizontal"))
			}
		}
		if (p$ncol < 3) 
			print(p[2,1] + 
						labs(title=title, x=gsub("\n"," ", column.labels[1]), y=gsub("\n", " ", column.labels[2])))
		else
			print(p)

		if (detailed)
			for (j in 1:p$ncol)
				for (i in 1:j)
					print(p[i,j]) 


		if (!is.na(ofle.tag))
			dev.off()

	}

test.get.svdplot <- 
	function() {
		source("./auxggplot2.R")
		data <- iris
		names(data) <- c("sepal_1", "sepal_2", "petal_1", "petal_2", "Species")
		get.svdplot(data, columns=1:4, ks=2, verbose=T)
		get.svdplot(data, columns=1:4, ks=3, verbose=T)

		library(zebrafishRNASeq)
		data(zfGenes)
		filter <- apply(zfGenes, 1, function(x) length(x[x>5])>=2)
		dta <- zfGenes[filter,]
		"
		library(RUVSeq)
		pd <- as.factor(sub('\\d+','',colnames(dta),perl=T))
		set <- newSeqExpressionSet(as.matrix(dta), 
															 phenoData=data.frame(pd, row.names=colnames(dta)))
		library(RColorBrewer)
		colors <- brewer.pal(3, 'Set2')
		plotRLE(set, outline=FALSE, ylim=c(-4, 4), col=colors[pd])
		plotPCA(set, col=colors[pd], cex=1.2)
		"

		param.preprocess <- list(ps=1, logf='log', col.mean.center=T) # RUVseq: Risso, 2014, Nat Biot
		dta.risso <- preprocess(dta=dta, param=param.preprocess)
		head(dta.risso)
		get.svdplot(dta.risso, ks=2, rep.sub='(\\d+)$', title='RUVseq: Risso, 2014, Nat Biot', verbose=T)
		get.svdplot(dta.risso, ks=4, rep.sub='(\\d+)$', title='RUVseq: Risso, 2014, Nat Biot', verbose=T)

		get.svdplot(dta.risso, ks=NULL, rep.sub='(\\d+)$', title='RUVseq: Risso, 2014, Nat Biot')
		get.svdplot(dta.risso, ks=c(2:5), rep.sub='(\\d+)$', title='RUVseq: Risso, 2014, Nat Biot', verbose=T)
		get.svdplot(dta.risso, ofle.tag="Risso_svd", ks=c(5,6,10), rep.sub='\\d+$', title='RUVseq: Risso, 2014, Nat Biot', verbose=T)

		# PAS_20161121
#		a <- read.table("./polyA5p_5p50nc.nuc", header=T);
#		b <- read.table("./polyA5p_3p50nc.nuc", header=T);
#		dta <- merge(a, b, by="X4_usercol", suffix=c(".polyA5p_5p", ".polyA5p_3p"));
#		dta <- dta[,grep("num_[ACGT]", names(dta))] / 50;
#		N <- nrow(dta)
#		get.svdplot(dta, ks=2, title=paste0("polyA5p", "N=", N), verbose=T);
#		head(dta)

	}

# http://genomicsclass.github.io/book/pages/pca_svd.html
# prcomp : svd
# dta    : t(dta)
# u : rotation
# v : x
# d^2/(ncol(dta)-1) : sdev
# d^2/sum(d^2)      : variance (?)
get.svdplot <- 
	function(data, columns=1:ncol(data), rep.sub='_(\\d+)$',
					 column.labels=names(data)[columns], param.preprocess=NULL,
					 ks=NULL,  max.pc=5, 
					 title="", verbose=F, ofle.tag=NA, ...) { 
		print(column.labels)

		if (!is.na(ofle.tag)) {
#			plot.new(); 
			pdf(paste(ofle.tag, ".pdf", sep=""), h = 8, w = 8);
		}

		# select columns
		data <- data[,columns]

		if (!is.null(param.preprocess)) 
			data <- preprocess(dta=data, param=param.preprocess)

		# perform svd
		s <- svd(t(data))

		ev <- preprocess(dta=s$d^2 / sum(s$d^2) * 100, param=list(round=2))
		max.pc <- min(max.pc, length(ev))

		if (is.null(ks))
			ks <- max.pc

		# adjust 
		if (any(which(ks >= max.pc)))  
			ks <- c(ks[ks < max.pc], max.pc)
		if (any(which(ks <= 2)))  
			ks <- c(2, ks[ks > 2])

		sp <- as.data.frame(s$u, row.names=column.labels)
		colnames(sp) <- paste0('PC', seq_along(s$d)) # col.names did not work above
		spw <- my_transform(sp, rep.sub=rep.sub, long=F)

		#redundant
		column.labels <- unlist(lapply(1:length(ev), 
																	 function(i) paste("PC", i, "\nExpl.var.(%):\n", as.character(ev[i]), sep='')))

		if (verbose) {
			print(c("ks:", ks))
			print(c("Variance explained(%):"));
			print(ev)
			print(sp)
			print(spw)
		}

		for (k in ks) {
			get.splom(spw, detailed=F, columns=1:k, column.labels=column.labels[1:k],
								mapping=aes_string(color="NameR", shape="NameR", label="ReplicateR"),
								lower=list(continuous=list(f="my_labdpoint", param=list(size=rel(4)))),
								diag=list(continuous=list(f="my_text", param=list(label="", size=rel(4)))),
								verbose=verbose,title=title)
		}

		if (!is.na(ofle.tag))
			dev.off()

	}

# TODO is it used???
test.setformula <- 
	function() {
		source("./auxggplot2.R")
		param <- list(valname="value")
	}

setformula <-
	function(param=NULL) {
		if (is.null(param))
			param <- list(valname="value", varname="variable")
		if (!any(grep("varname", names(param)))) {
			frml <- '~get(param$valname)'
		} else {
			frml <- 'get(param$valname) ~ get(param$varname)'
		}
		if (any(grep("panname", names(param)))) {
			frml <- paste(frml, ' | as.factor(get(param$panname))')
		}

		as.formula(frml)
	}

test.get.plot <- 
	function() {
		source("./auxggplot2.R")
		dta <- rbind(data.frame(data=c(rnorm(9, 6, 1),10), class="class1"), 
								 data.frame(data=rnorm(50, 15, 1), class="class2"))
		dta$split <- rep(c("panel1","panel2"),30)  
		head(dta)

		get.plot(dta=dta, param.formula=list(varname="class", valname="data", panname="split"), 
						 param.plot=list(chart="bwplot", cex=1, ylab="data", scales=list(x=list(rot=60)))) 
		get.plot(dta=dta, param.formula=list(varname="class", valname="data"),
						 param.plot=list(chart="bwplot", col="blue", ylab="data", scales=list(x=list(rot=60)))) 
		get.plot(dta=dta, param.formula=list(varname="class", valname="data"),
						 param.plot=list(chart="bwplot", col="blue", ylab="data", jitter.x=T,factor=0.1, scales=list(x=list(rot=60))), 
						 w.xyplot=T)

	}

# TODO make a generic ggplot function
get.plot <- 
	function(p, ofle=NULL, h=8, w=8) {

		if (! is.null(ofle)) {
			if (any(grep("pdf$", ofle))) {
				f <- get('pdf')	
			} else {
				f <- get('png')
			}
			plot.new(); f(ofle, h=h, w=w)
		}
		print(p)
		if (! is.null(ofle))
			dev.off()
	}

# TODO minimize data processing in function calls
# TODO generalize multiple catrgorial data columns; use color and shape and maybe fill
# TODO scale axis than log transforming data
# TODO check the pattern in Ling's data
get.EDA <-
	function(data, dataa=NA, columns=1:ncol(data), filter=NA, 
					 title="", ...) {

		ggplot <- set.ggplot()
		# remove
		mapping <- NULL
		facet <- NA
		title <- paste0(title, ";N=", nrow(data))
		if ((!is.na(dataa)) && nrow(dataa) > 1) { 
			mapping <- aes_string(color=names(dataa)[1])
			facet <- names(dataa)[1]
#			if (nrow(dataa) > 3) {
#				mapping <- structure(mergeSetwDefault(mapping, aes(shape=names(dataa)[2])))
#			}
		}

		print(paste("Doing", title))

		print("my_boxplot")
		my_boxplot(data, columns=columns, ncol=1, title=title, ...) %>% print
		if (!is.na(facet))
			my_boxplot(data, columns=columns, facet=facet, ncol=1, title=title, ...) %>% print
		print("my_ecdf")
		my_ecdf(data, columns=columns, ncol=1, ..., combine=T, round=F, title=title) %>% print
		if (!is.na(facet))
			my_ecdf(data, columns=columns, facet=facet, ncol=1, ..., combine=T, round=F, title=title) %>% print
		print("get.splom")
		get.splom(data, columns=columns, mapping=mapping, detailed=F, title=title, ...)

		# ungrouped data
		print("get.svdplot")
		get.svdplot(data, columns=columns, ks=c(2,4), size=rel(6), title=title, ...)
		print("get.corrheat_seg_boxplot")
		get.corrheat_seg_boxplot(data, columns=columns, title=title, ...)
		if ((!is.na(dataa)) && nrow(dataa) > 1) { 
			for (i in 1:nrow(dataa)) {
				filterv <- dataa[i, names(dataa)[1]]
				filter <- paste0(names(dataa)[1], "==", "'", filterv, "'")
				print(i)
				print(filter)
				data %>% filter_(filter) %>% head %>% print
				datac <- data %>% filter_(filter)
				title1 <- paste0(title,";",filterv,";N=", nrow(datac))
				print("get.svdplot")
				get.svdplot(datac, columns=columns, ks=c(2,4), size=rel(6), title=title1, ...)
				print("get.corrheat_seg_boxplot")
				get.corrheat_seg_boxplot(datac, columns=columns, title=title1, ...)
			}
		}
	}


