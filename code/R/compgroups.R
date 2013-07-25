#!/usr/bin/R
TD <- read.delim('~/Labwork/Rwork/TajimaD',header=F,sep=" ")
mat <- matrix(unlist(TD),ncol=2,byrow=F)
dTD <- data.frame(mat[,1],as.numeric(mat[,2]))
names(dTD) <- c('Gene','TajimaD')

#SPECIES TOGETHER:
sim <- read.table('~/Labwork/Rwork/TDMKtab.sim')
yak <- read.table('~/Labwork/Rwork/TDMKtab.yak')
sim[,c(1,5,10)] -> dsim
yak[,c(1,5,10)] -> dyak
names(dyak) <- c('Gene','yAlpha','yakpvals')
names(dsim) <- c('Gene','sAlpha','simpvals')
merge(dyak,dsim) -> BOTH
dAlpha <- BOTH$sAlpha - BOTH$yAlpha
cbind(BOTH,dAlpha) -> BOTH

groups <- read.table('~/Labwork/Annotation/Lists/Gene.Groups',header=F,sep=" ")
data.frame(groups) -> groups
names(groups) <- c('Gene','Group')
factor(groups$Gene) -> groups$Gene
factor(groups$Group) -> groups$Group       
merge(groups,yak) -> yakg
yakg[yakg$Group != 'GPCR',] -> yakg

library(ggplot2)
ggplot(yakg,aes(Group,Alpha)) + geom_boxplot()
merge(dTD,groups) -> TDg
ggplot(yakg,aes(Group,TajimaD)) + geom_boxplot()

#For DoFE analysis:
out <- yakg[,c(1,2,8,10,7,9)]
#Headers:
cat("Phosphatase Genes:\n",file="~/Labwork/Annotation/Lists/outP")
cat("Trans-Membrane Genes:\n",file="~/Labwork/Annotation/Lists/outTM")
cat("Transcription Factors:\n",file="~/Labwork/Annotation/Lists/outTF")
cat("Kinase Genes:\n",file="~/Labwork/Annotation/Lists/outK")
cat("Ubiquitin Genes:\n",file="~/Labwork/Annotation/Lists/outU")
#Files:
write.table(out[out$Group == 'Phosphatase',c(1,3:6)],'~/Labwork/Annotation/Lists/outP',row.names=F,col.names=F,append=TRUE,sep="\t")
write.table(out[out$Group == 'Transmembrane',c(1,3:6)],'~/Labwork/Annotation/Lists/outTM',row.names=F,col.names=F,append=TRUE,sep="\t")
write.table(out[out$Group == 'TranscriptionFactors',c(1,3:6)],'~/Labwork/Annotation/Lists/outTF',row.names=F,col.names=F,append=TRUE,sep="\t")
write.table(out[out$Group == 'Kinase',c(1,3:6)],'~/Labwork/Annotation/Lists/outK',row.names=F,col.names=F,append=TRUE,sep="\t")
write.table(out[out$Group == 'UbiquitinSet',c(1,3:6)],'~/Labwork/Annotation/Lists/outU',row.names=F,col.names=F,append=TRUE,sep="\t")

#Dofe work, after running the program:
dofe <- read.table('~/Labwork/Rwork/dofetab',header=F,sep=" ")
names(dofe) <- c('Group','n','Alpha','SE','C1','C2','mlAlpha','LLV','LLR','p.value','mlC1','mlC2')

library(ggplot2)
ggplot(dofe,aes(Group,Alpha)) + geom_point() + geom_errorbar(aes(ymin=C1,ymax=C2),color='black',width=.1)

dofe[,c(1,3,5,6,7,11,12)] -> DOFE
DOFE <- reshape(DOFE,idvar="Gene",direction="long",times=c('Smith and Eyre-Walker','Bierne and Eyre-Walker'),varying=list(ALPHA=c(2,5),C1=c(3,6),AS=c(4,7)),v.names=c("Alpha","C1","C2"))
rownames(DOFE) <- 1:nrow(DOFE)

#png('DoFEanalysis.png',width=700,height=500)
ggplot(DOFE,aes(Group,Alpha)) + geom_point(pch=16) + geom_errorbar(aes(ymin=C1,ymax=C2),color='black',width=.1) + facet_wrap(~time)
#dev.off()
