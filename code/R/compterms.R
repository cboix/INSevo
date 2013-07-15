#!/usr/bin/R
library(ggplot2)
TDMK <- read.table('~/Labwork/Rwork/TDMKtab')
GO <- read.delim('~/Labwork/Annotation/termgenes',header=F,sep=" ")
names(GO) <- c('Gene','Term')
mer <- merge(GO,TDMK)

ggplot(mer,aes(Gene,


# using the reduced geneall:
read.delim('~/Labwork/Annotation/red10geneall',header=F,sep=" ") -> ALL
names(ALL) <- c('Gene','Term')
merge(ALL,TDMK) -> MERGEALL

#AVERAGE P.VALUE per term:
with(MERGEALL,by(pvalsout,Term,mean)) -> means
with(MERGEALL,by(Alpha,Term,median)) -> ameans
as.numeric(as.character(means)) -> pvalmeans
as.numeric(as.character(ameans)) -> ameans
data.frame(cbind(levels(MERGEALL$Term),pvalmeans,ameans)) -> dd
dd$pvalmeans <- as.numeric(as.character(dd$pvalmeans))
dd$ameans <- as.numeric(as.character(dd$ameans))
names(dd) <- c('Term','Mean.p.value','Alpha.median')
MERGEALL <- merge(MERGEALL,dd)

ggplot(MERGEALL,aes(reorder(Term,Alpha.median),Alpha,color=Mean.p.value)) + geom_boxplot() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + labs(x='GO Term',y='Alpha')





