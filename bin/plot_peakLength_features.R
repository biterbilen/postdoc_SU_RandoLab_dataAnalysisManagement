source("~/Projects/biter_biter_shared/bin/bigData.R")
source("~/Projects/biter_biter_shared/bin/auxggplot2.R")

rep.sub <- "_(\\d+)$"
# 1.
# All peak lengths
fle <- "length.txt.gz"
dta <- readBigTable(fle, col.names=c("Length", "Tag"))
dta$Tag <- gsub("_peak_\\d*", "", dta$Tag)
dta$Tag <- as.factor(gsub("_young_TMX_SC_NA", "", dta$Tag))
dta %<>% mutate(Name=gsub(rep.sub, "", dta$Tag), Replicate=getRepNumber(dta$Tag, rep.sub))
title1 <- paste0("All Peaks\n",paste.namedList(summary(dta$Tag),sp="=", sp2="\n"))
plot.new(); pdf(gsub(".txt.gz", ".pdf", fle), h=8,w=8)
my_ecdf(dta, aes(x=log10(Length), color=Name), facet="Replicate", long=F, nrow=2, size=rel(2)) +
	labs(title=title1)
#dev.off()

# 2
# CKO and WT specific peak lengths
fle <- "unique_length.txt.gz"
dta <- readBigTable(fle, header=T)
dta %<>% mutate(Name=gsub(rep.sub, "", dta$Tag), Replicate=getRepNumber(dta$Tag, rep.sub))
head(dta)
#plot.new(); pdf(gsub(".txt.gz", ".pdf", fle), h=8,w=8)
title1 <- paste0("Specific Peaks\n",paste.namedList(summary(dta$Tag),sp="=", sp2="\n"))
my_ecdf(dta, aes(x=log10(Length), color=Name), facet="Replicate", long=F, nrow=2, size=rel(2)) +
	labs(title=title1)
#dev.off()

# 3.
# CKO and WT Commmon peak lengths
fle <- "common_length.txt.gz"
dta <- readBigTable(fle, header=T)
head(dta)
#plot.new(); pdf(gsub(".txt.gz", ".pdf", fle), h=8,w=8)
title1 <- paste0("Common Peaks\n",paste.namedList(summary(dta$Tag),sp="=", sp2="\n"))
my_bin(dta %>% sample_frac(0.2) , aes(x=log10(HrCKOlength), y=log10(HrWTlength)), facet="Tag", long=F, nrow=2, size=rel(2), bins=100) +
	labs(title=title1)
dev.off()


