#!/usr/bin/R

PFAM <- read.table('~/Labwork/Annotation/PFAM/MERGE.all')
head(PFAM)
EXP <- read.table('EXPtabsim')
head(EXP)

levels(PFAM$Gene) -> genes
count <- c();for(i in genes){sum(PFAM$Gene == i)-> s;count <- c(count,s)}
data.frame(cbind(genes,count)) -> domains
names(domains) <- c('Gene','Num.Domains')
library(reshape)
melt(EXP,measure.vars=c('SY.p.value','NS.p.value')) -> EXP2
EXP2 <- merge(EXP2,domains)
library(ggplot2)
ggplot(EXP2,aes(Num.Domains,value)) + geom_boxplot() + facet_wrap(~variable)


