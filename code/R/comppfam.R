#!/usr/bin/R
#The gene data.frame:
PFAM <- read.table('~/Labwork/Annotation/PFAM/MERGE.all')
TD <- read.delim('TajimaD',header=F,sep=" ")
mat <- matrix(unlist(TD),ncol=2,byrow=F)
TD <- data.frame(mat[,1],as.numeric(mat[,2]))
names(TD) <- c('Gene','TajimaD')
MERGETD <- merge(PFAM,TD)
non <- 1*(MERGETD$Clan.Desc == 'None')
MERGETD <- cbind(MERGETD,non)
factor(MERGETD$non) -> MERGETD$non
library(ggplot2)
ggplot(MERGETD,aes(reorder(Clan.Desc,TajimaD),TajimaD)) + geom_boxplot() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + labs(title="PFAM ID vs. Tajima's D value",x='PFAM ID',y="Tajima's D") + coord_flip() 
ggplot(MERGETD,aes(non,TajimaD)) + geom_boxplot()


species <- 'yak'
createPFAM <- function(species,PFAM){

    TDMK <- read.table(paste('~/Labwork/Rwork/TDMKtab.',species,sep=""))
    MERGEPFAM <- merge(PFAM,TDMK)
    non <- 1*(MERGEPFAM$Clan.Desc == 'None')
    MERGEPFAM <- cbind(MERGEPFAM,non)
    factor(MERGEPFAM$non) -> MERGEPFAM$non
#significantly different Dn/Ds
    MERGEALL <- merge(MERGETD,MERGEPFAM)
   
    with(MERGEALL,by(pvalsout,Clan.Desc,mean)) -> means
	with(MERGEALL,by(Alpha,Clan.Desc,median)) -> amedian
	as.numeric(as.character(means)) -> pvalmeans
	as.numeric(as.character(amedian)) -> amedian
	data.frame(cbind(levels(MERGEALL$Clan.Desc),pvalmeans,amedian))-> dd
	dd$pvalmeans <- as.numeric(as.character(dd$pvalmeans))
	dd$amedian <- as.numeric(as.character(dd$amedian))
	names(dd) <- c('Clan.Desc','Mean.p.value','Alpha.median')
	MERGEALL <- merge(MERGEALL,dd)

    # Reduce: (also we need to reduce genes (so points aren't plotted twice over...)
    levels(MERGEALL$Clan.Desc) -> desc
    # For any numbers:
    n <- 5
    dnum <- c()
    for(i in desc){if(sum(MERGEALL$Clan.Desc == i) >= n){dnum <- c(dnum,i)}} 
    redMERGE <- MERGEALL[MERGEALL$Clan.Desc %in% dnum,]
    # make sure that no gene is expressed twice!
    redMERGE2 <- redMERGE[,c(1,2,10:22)]
    redMERGE2 <- unique(redMERGE2)
    redMERGE2 <- redMERGE2[redMERGE2$Alpha > -20,]
    library(reshape)
    redMERGE3 <- melt(redMERGE2,measure.vars=c("TajimaD","Alpha"))
    png(paste('PFAMvTDandMK.',n,'.png',sep=""),width=5400/n,height=1200)
    ggplot(redMERGE3,aes(reorder(Clan.Desc,Alpha.median),value,color=Mean.p.value)) + geom_boxplot() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + labs(title=paste("PFAM ID v. MK and TD value for ",species," for Clans with ",n," or more members",sep=""),x='PFAM ID') + facet_wrap(~variable,ncol=1,scales='free_y')
    dev.off()

    ggplot(redMERGE2,aes(reorder(Clan.Desc,Alpha.median),Alpha,color=Mean.p.value)) + geom_boxplot() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + labs(title=paste("PFAM ID v. MK value for ",species," for Clans with ",n," or more members",sep=""),x='PFAM ID',y='Alpha (outliers below -20 removed)') + coord_flip() + scale_y_continuous(limits=c(-20,2))

    #--------------------------------------------------------------------------------------------------
    # With the exptab:
    EXP <- read.table(paste("~/Labwork/Rwork/EXPtab",species,sep=""))
    ALL <- merge(MERGEALL,EXP)
     # For any numbers:
    n <- 3
    dnum <- c()
    for(i in desc){if(sum(ALL$Clan.Desc == i) >= n){dnum <- c(dnum,i)}} 
    redEXP <- ALL[ALL$Clan.Desc %in% dnum,]
    # make sure that no gene is expressed twice!
    redEXP2 <- redEXP[,c(1,2,10:28)]
    redEXP2 <- unique(redEXP2)
    redEXP2 <- redEXP2[redEXP2$Alpha > -20,]
    with(redEXP2,by(NS.p.value,Clan.Desc,mean)) -> means
	as.numeric(as.character(means)) -> NSpv.mean
	data.frame(cbind(levels(MERGEALL$Clan.Desc),NSpv.mean))-> dd
	dd$NSpv.mean <- as.numeric(as.character(dd$NSpv.mean))
	names(dd) <- c('Clan.Desc','NSpv.mean')
    redEXP2 <- merge(redEXP2,dd)

    ggplot(redEXP2,aes(reorder(Clan.Desc,NSpv.mean),NS.p.value,color=NSpv.mean)) + geom_boxplot() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + labs(title=paste("PFAM ID v. pvalue for GFEXP in ",species," for Clans with ",n," or more members",sep=""),x='PFAM ID',y="GFexp p.value") + coord_flip()

    library(reshape)
    redEXP3 <- melt(redEXP2,measure.vars=c("Alpha","NS.p.value"))
    ggplot(redEXP3,aes(reorder(Clan.Desc,Alpha.median),value,color=Mean.p.value)) + geom_boxplot() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + labs(title=paste("PFAM ID v. MK and GFexp value for ",species," for Clans with ",n," or more members",sep=""),x='PFAM ID') + facet_wrap(~variable,ncol=1,scales='free_y')
   



}

createPFAM('yak',PFAM)
createPFAM('yak',PFAM)

# We need to create a heuristic to measure separation between different PFAM groups (and also between different GO term groups). (I could wrap by group if there aren't too many) -- by alpha medians?
