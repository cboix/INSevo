#!/usr/bin/R
library(ggplot2)
TDMK <- read.table('~/Labwork/Rwork/TDMKtab')


GO <- read.delim('~/Labwork/Annotation/termgenes',header=F,sep=" ")
names(GO) <- c('Gene','Term')
mer <- merge(GO,TDMK)

#GO terms:
TERMS <- read.delim('~/Labwork/Annotation/short.obo',header=F,sep="\t")
names(TERMS) <- c('Term','Name','Namespace')

# using the reduced geneall:
read.delim('~/Labwork/Annotation/red10geneall',header=F,sep=" ") -> ALL
names(ALL) <- c('Gene','Term')
ALL <- merge(ALL,TERMS)

MERGEALL <- merge(ALL,TDMK)

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

# png('~/Labwork/Rwork/Boxplots2.png',height=2000,width=800)
ggplot(MERGEALL,aes(reorder(Name,Alpha.median),Alpha,color=Mean.p.value)) + geom_boxplot() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + labs(x='GO Name',y='Alpha (outliers below -20 removed)') + coord_flip() + facet_wrap(~Namespace) + scale_y_continuous(limits=c(-20,2))


#WITH the annoying outliers which shrink the plot:
# png('~/Labwork/Rwork/Boxplotsfull.png',height=2000,width=1000)
ggplot(MERGEALL,aes(reorder(Name,Alpha.median),Alpha,color=Mean.p.value)) + geom_boxplot() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + labs(x='GO Name',y='Alpha (outliers below -20 removed)') + coord_flip() + facet_wrap(~Namespace)

png('~/Labwork/Rwork/Boxplots3TD.png',height=2000,width=1000)
ggplot(MERGEALL,aes(reorder(Name,TajimaD),TajimaD,color=Mean.p.value)) + geom_boxplot() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + labs(x='GO Name',y='Alpha (outliers below -20 removed)') + coord_flip() + facet_wrap(~Namespace)
