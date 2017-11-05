source("~/Projects/biter_biter_shared/bin/bigData.R")
source("~/Projects/biter_biter_shared/bin/auxggplot2.R")
#textsize=100
#size=2
#ggplot <- set.ggplot(textsize=textsize)
ggplot <- set.ggplot()


#p1 <- "old_NA"
#p2 <- "young_NA"
#oflet <- "LM"
p1 <- commandArgs(trailingOnly=T)[1]
p2 <- commandArgs(trailingOnly=T)[2]
oflet <- commandArgs(trailingOnly=T)[3]
mafle <- "tfbs.txt.gz" 
mfle <- "featureCounts"; norm <- F
#mafle <- "tmp2.bed"
#mfle <- "TM_allPromoters.tab.gz"
print(oflet)


# methylation rates w methylation coordinates
if (norm) {
m <- readBigTable(mfle, header=T)
head(m,1)
source('~/Projects/biter_biter_shared/bin/auxggplot2.R')
spike <- list(id='ERCC', rm=F) #why is it set?
param.DE <- list(method="edgeRDidNotWork", FDRcut=0.01, contrast=c(1,1,0,0,-1,-1))
caller.profileNormalization.DE(tag="peaks", dta=dta, param.DE=param.DE, spike=spike, species=NA, scale="none")
caller.profileNormalization.DE(tag="peaks", dta=dta, param.DE=param.DE, spike=spike, species=NA, scale="upper")
quit()
}

# choose the best replicate separating one
mfle <- "peaks_scupper_RSRUVr_count.txt.gz"
m <- readBigTable(mfle, header=T)
samples <- grep("_SC_", names(m), value=T)
dta <- m[samples]; rownames(dta) <- m$name;

a <- readBigTable(mafle, header=T); 

ma <- merge(a,m)
summary(ma[samples])
head(ma,1)
ps <- 30 # 1st quadrant
p1n <- grep(p1, names(ma)); ma$p1 <- rowSums(ma[,p1n])/length(p1n)
p2n <- grep(p2, names(ma)); ma$p2 <- rowSums(ma[,p2n])/length(p2n)

ma$acc <- log2(ma$p1+ps) - log2(ma$p2+ps)

plot.new(); pdf(paste(oflet, ".pdf", sep=""), h=8, w=8)
#hist(ma$pm.A4 - ma$pm.Q4)
ggplot(ma, aes(x=log2(p1+ps), y=log2(p2+ps))) + geom_bin2d() + labs(title="Accessibility - post count adjustment")
ggplot(ma, aes(x= log2(ma$p1+ps) - log2(ma$p2+ps))) + geom_density(size=2) + labs(title="Accessibility - post count adjustment", y="Density")

frml <- paste("acc ~", paste(grep(".MA", names(ma), value=T), collapse=" + "), collapse="")
# Check Lasso regression with other factors
#g <- glm(exp ~ met, data=dta1, family="gaussian")
g <- lm(as.formula(frml), data=ma)
k <- options("max.print")[[1]]
options(max.print=100000)
print(summary(g))
options(max.print=k)
#print(cor.test(dta$exp, dta$met))
#print(cor.test(dta1$exp, dta1$met))
writeData(ma, paste(oflet, ".in4pred.gz", sep=""), row.names=F)
dev.off()
quit()


#####################################
#DUMP
library(h2o)
h2o.init()
#h2o.shutdown()
#path <- system.file("extdata", "prostate.csv", package="h2o")
path <- system.file("extdata", "prostate.csv", package="h2o")
h2o_df <- h2o.importFile(path)
h2o_df$CAPSULE <- as.factor(h2o_df$CAPSULE)
#alpha_opts <- list(list(0), list(0.25), list(0.5), list(0.75), list(1))
alpha_opts <- c(0.05, 0.25, 0.5, 0.75, 0.95)
hyper_parameters <- list(alpha = alpha_opts)
grid <- h2o.grid("glm", hyper_params=hyper_parameters, y="CAPSULE", x=c("AGE", "RACE", "PSA", "GLEASON"), training_frame=h2o_df, family="binomial")
grid_models <- lapply(grid@model_ids, function(model_id) model = h2o.getModel(model_id)) 
for (i in 1:length(grid_models)) {
	print(sprintf("regularization:%-50s auc: %f", 
								grid_models[[i]]@model$model_summary$regularization,
								h2o.auc(grid_models[[i]])))
}

data(Insurance)
library(MASS)

class(Insurance$Group) <- "factor"
class(Insurance$Age) <- "factor"
h2o_df <- as.h2o(Insurance)
poisson.fit <- h2o.glm(y="Claims", x=c("District", "Group", "Age"), training_frame = h2o_df, family="poisson")
g <- glm(Claims ~ District + Group + Age + offset(log(Holders)),
		    data = Insurance, family = poisson)
summary(g)

quine.nb1 <- glm.nb(Days ~ Sex/(Age + Eth*Lrn), data = quine)
quine.nbA <- glm.convert(quine.nb1)
quine.nbB <- update(quine.nb1, . ~ . + Sex:Age:Lrn)
anova(quine.nbA, quine.nbB)

quit()

library (RegulatorInference)
regulatorinference.test.installation
regulatorinference.test.installation('sample_data')
dir()
sbs.run.regulator.inference
sbs.run.regulator.inference()
quit()
savehistory("history.R")
load("./sample_data/sample_inputs.rda")
ls()

source("~/Projects/biter_biter_shared/bin/bigData.R")
source("~/Projects/biter_biter_shared/bin/auxggplot2.R")

ggplot <- set.ggplot()

pth <- "./out.exp_ordered_Rle/"
pth <- commandArgs(trailingOnly=T)[1]
oflet <- "all"

plot.new(); pdf(paste(oflet, ".pdf", sep=""), h=8, w=8)

multmerge = function(pth){
	fles  <- list.files(path=pth, ".gz$", full.names=TRUE)
	tagl  <- lapply(fles, function(x){ strsplit(basename(x),'\\.')[[1]][1] })
	dtal  <- lapply(fles, function(x){ readBigTable(x, header=T) })
	dtal <- lapply(1:length(dtal), function(i){ names(dtal[[i]])[3] <- paste(names(dtal[[i]])[3], tagl[i], sep="."); dtal[[i]] })
	Reduce(function(x,y) {merge(x,y, by=c("Id", "exp"))}, dtal)
}

dta <- multmerge(pth); # rownames(dta) <- dta$Id; dta$Id <- NULL
writeData(dta, paste(oflet, ".in4pred.gz", sep=""), row.names=F)
dim(dta)
ggplot(dta, aes(x=exp)) + geom_density()
get.splom(dta[,-1], title="before RLE")

dta <- dta[order(dta$exp),]
frml <- paste("exp ~", paste(grep("^met", names(dta), value=T), collapse=" + "), collapse="")
frml1 <- paste("exp ~", paste(grep("^met", setdiff(names(dta), "met.TM_FOSL1_0"), value=T), collapse=" + "), collapse="")
frmli <- paste("exp ~", paste(grep("^met", names(dta), value=T), collapse=" : "), collapse="")
frmli1 <- paste("exp ~", paste(grep("^met", setdiff(names(dta), "met.TM_FOSL1_0"), value=T), collapse=" : "), collapse="")

g <- lm(as.formula(frml), data=dta, family="gaussian")
g1 <- lm(as.formula(frml1), data=dta, family="gaussian")
gi <- lm(as.formula(frmli), data=dta, family="gaussian")
gi1 <- lm(as.formula(frmli1), data=dta, family="gaussian")
AIC(g,g1,gi,gi1)
summary(g)
#dta <- dta[order(dta$met),] # results are similar when ordered based on met and RLe run
dev.off()

dta1 <- dta["Id"]
dta1$exp <- as.numeric(runmean(Rle(dta$exp), k=20, endrule="constant"))
dta1$met <- as.numeric(runmean(Rle(dta$met), k=20, endrule="constant"))
head(dta1,2)

ggplot(dta1, aes(y=exp, x=met)) + geom_smooth()

# Check Lasso regression with other factors
g <- glm(exp ~ met, data=dta1, family="gaussian")
summary(g,2)
print(cor.test(dta$exp, dta$met))
print(cor.test(dta1$exp, dta1$met))
dev.off()
quit()


library(h2o)
h2o.init()
#h2o.shutdown()
#path <- system.file("extdata", "prostate.csv", package="h2o")
path <- system.file("extdata", "prostate.csv", package="h2o")
h2o_df <- h2o.importFile(path)
h2o_df$CAPSULE <- as.factor(h2o_df$CAPSULE)
#alpha_opts <- list(list(0), list(0.25), list(0.5), list(0.75), list(1))
alpha_opts <- c(0.05, 0.25, 0.5, 0.75, 0.95)
hyper_parameters <- list(alpha = alpha_opts)
grid <- h2o.grid("glm", hyper_params=hyper_parameters, y="CAPSULE", x=c("AGE", "RACE", "PSA", "GLEASON"), training_frame=h2o_df, family="binomial")
grid_models <- lapply(grid@model_ids, function(model_id) model = h2o.getModel(model_id)) 
for (i in 1:length(grid_models)) {
	print(sprintf("regularization:%-50s auc: %f", 
								grid_models[[i]]@model$model_summary$regularization,
								h2o.auc(grid_models[[i]])))
}

data(Insurance)
library(MASS)

class(Insurance$Group) <- "factor"
class(Insurance$Age) <- "factor"
h2o_df <- as.h2o(Insurance)
poisson.fit <- h2o.glm(y="Claims", x=c("District", "Group", "Age"), training_frame = h2o_df, family="poisson")
g <- glm(Claims ~ District + Group + Age + offset(log(Holders)),
		    data = Insurance, family = poisson)
summary(g)

quine.nb1 <- glm.nb(Days ~ Sex/(Age + Eth*Lrn), data = quine)
quine.nbA <- glm.convert(quine.nb1)
quine.nbB <- update(quine.nb1, . ~ . + Sex:Age:Lrn)
anova(quine.nbA, quine.nbB)

quit()

library (RegulatorInference)
regulatorinference.test.installation
regulatorinference.test.installation('sample_data')
dir()
sbs.run.regulator.inference
sbs.run.regulator.inference()
quit()
savehistory("history.R")
load("./sample_data/sample_inputs.rda")
ls()
