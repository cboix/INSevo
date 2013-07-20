#!/usr/bin/R
#The gene data.frame:
PFAM <- read.delim('~/Labwork/Annotation/pfam',header=F,sep=" ")
names(PFAM) <- c('Gene','ID','Class')

createPFAM <- function(species,PFAM){
    TDMK <- read.table(paste('~/Labwork/Rwork/TDMKtab.',species,sep=" "))
    MERGEPFAM <- merge(PFAM,TDMK)

    with(MERGEALL,by(pvalsout,ID,mean)) -> means
	with(MERGEALL,by(Alpha,ID,median)) -> amedian
	as.numeric(as.character(means)) -> pvalmeans
	as.numeric(as.character(amedian)) -> amedian
	data.frame(cbind(levels(MERGEALL$ID),pvalmeans,amedian))-> dd
	dd$pvalmeans <- as.numeric(as.character(dd$pvalmeans))
	dd$amedian <- as.numeric(as.character(dd$amedian))
	names(dd) <- c('ID','Mean.p.value','Alpha.median')
	MERGEALL <- merge(MERGEALL,dd)

    name <- paste('~/Labwork/Rwork/BoxplotsMK.pfam.',species,'.png',sep="")
    png(name,height=2000,width=1000)
    library(ggplot2)
    ggplot(MERGEALL,aes(reorder(ID,Alpha.median),Alpha,color=Mean.p.value)) + geom_boxplot() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + labs(title=paste("PFAM ID v. MK value for ",species,sep=""),x='PFAM ID',y='Alpha (outliers below -20 removed)') + coord_flip() + facet_wrap(~Class) + scale_y_continuous(limits=c(-20,2))
	dev.off()
	
    name1 <- paste('~/Labwork/Rwork/BoxplotsfullMK.pfam.',species,'.png',sep="")
	#WITH the annoying outliers which shrink the plot:
	png(name1,height=2000,width=1000)
	ggplot(MERGEALL,aes(reorder(ID,Alpha.median),Alpha,color=Mean.p.value)) + geom_boxplot() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + labs(title=paste("PFAM ID v. MK value for ",species,sep=""),x='PFAM ID',y='Alpha (no outliers removed)') + coord_flip() + facet_wrap(~Class)
	dev.off()
	
    name2 <- paste('~/Labwork/Rwork/BoxplotsTD.pfam.',species,'.png',sep="")
	png(name2,height=2000,width=1000)
	ggplot(MERGEALL,aes(reorder(ID,TajimaD),TajimaD,color=Mean.p.value)) + geom_boxplot() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + labs(title=paste("PFAM ID v. Tajima's D for ",species,sep=""),x='PFAM ID',y='Alpha (outliers below -20 removed)') + coord_flip() + facet_wrap(~Class)
	dev.off()

}

createPFAM('yak',PFAM)
createPFAM('yak',PFAM)

# We need to create a heuristic to measure separation between different PFAM groups (and also between different GO term groups). (I could wrap by group if there aren't too many) -- by alpha medians?
