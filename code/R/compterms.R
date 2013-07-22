#!/usr/bin/R
#First, the term data.frame:
TERMS <- read.delim('~/Labwork/Annotation/GO/short.obo',header=F,sep="\t")
names(TERMS) <- c('Term','Name','Namespace')
# using the reduced geneall:
read.delim('~/Labwork/Annotation/GO/red10geneall',header=F,sep=" ") -> ALL
names(ALL) <- c('Gene','Term')
ALL <- merge(ALL,TERMS)

#Then, a function to create the boxplots:
createGO <- function(species,ALL){
    #species <- 'yak'
	library(ggplot2)
	TDMK <- read.table(paste('~/Labwork/Rwork/TDMKtab.',species,sep=""))
	#GO terms:
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
	
    name <- paste('~/Labwork/Rwork/BoxplotsMK.',species,'.png',sep="")
	png(name,height=2000,width=800)
	ggplot(MERGEALL,aes(reorder(Name,Alpha.median),Alpha,color=Mean.p.value)) + geom_boxplot() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + labs(title=paste("GO term v. MK value for ",species,sep=""),x='GO Name',y='Alpha (outliers below -20 removed)') + coord_flip() + facet_wrap(~Namespace) + scale_y_continuous(limits=c(-20,2))
	dev.off()
	
    name1 <- paste('~/Labwork/Rwork/BoxplotsfullMK.',species,'.png',sep="")
	#WITH the annoying outliers which shrink the plot:
	png(name1,height=2000,width=1000)
	ggplot(MERGEALL,aes(reorder(Name,Alpha.median),Alpha,color=Mean.p.value)) + geom_boxplot() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + labs(title=paste("GO term v. MK value for ",species,sep=""),x='GO Name',y='Alpha (no outliers removed)') + coord_flip() + facet_wrap(~Namespace)
	dev.off()
	
    name2 <- paste('~/Labwork/Rwork/BoxplotsTD.',species,'.png',sep="")
	png(name2,height=2000,width=1000)
	ggplot(MERGEALL,aes(reorder(Name,TajimaD),TajimaD,color=Mean.p.value)) + geom_boxplot() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + labs(title=paste("GO term v. Tajima's D for ",species,sep=""),x='GO Name',y='Alpha (outliers below -20 removed)') + coord_flip() + facet_wrap(~Namespace)
	dev.off()
}

createGO('yak',ALL)
createGO('sim',ALL)
