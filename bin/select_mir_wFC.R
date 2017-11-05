source("../../../biter_biter_shared/bin/bigData.R")
source("../../../biter_biter_shared/bin/auxggplot2.R")

otag <- "selected"
dta <- readBigTable("./Data", comment.char="", header=T) 
head(dta)
psa <- min(dta[dta$A>0,"A"])
psbl <- min(dta[dta$BL>0,"BL"])
dta$A <- dta$A + psa
dta$BL <- dta$BL + psbl
dta$logFC <- log(dta$A) - log(dta$BL)
dta$logSum <- log(dta$A) + log(dta$BL)

filtert <- 'logFC > 1.5 | logFC < -1.5'
commentLine <- paste0("A and BL is ps added;A_ps=", psa, ";BL_ps=", psbl, ";ordering:abs(logFC)")
dta %<>% arrange(desc(abs(logFC)))
head(dta)
writeData(dta, paste0("all", ".txt.gz"), row.names=F, commentLine=commentLine)
writeData(dta %>% filter_(filtert), paste0(otag, ".txt.gz"), row.names=F, 
					commentLine=paste0(commentLine,";filtering:", filtert))

plot.new(); pdf(paste0(otag, ".pdf"), h=8, w=8)
my_bin(dta, aes(y=logFC, x=logSum, label=Geneid), long=F, bin=51, slope=0, intercept=0, 
			 textRepelFilter=filtert) +
			 geom_abline(slope=0, intercept=1.5, color="darkgray", linetype=2, size=rel(1)) +
			 geom_abline(slope=0, intercept=-1.5, color="darkgray", linetype=2, size=rel(1)) +
			 labs(x="log(A+ps) + log(BL+ps)", y="log(A+ps) - log(BL+ps)")

dev.off()



