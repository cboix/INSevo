#!/usr/bin/R
#For this, we will look at yakuba data ONLY and see if we can segregate by statistic and see any type of GO enrichment/PFAM domains (we use Merge(all=TRUE) for this.
#MOTIF ENRICHMENT:
read.delim('IN',header=F,sep=" ") -> IN
names(IN) <- c('Count','Gene','Set')
sums <- with(IN,by(Count,Set,sum))
sums <- data.frame(matrix(c('Control','Insulin',as.numeric(as.character(sums))),ncol=2))
names(sums) <- c('Set','Total')


#First create the huge annotated table:
SEMK <- read.table('~/Labwork/Rwork/SEMKproc.yakfulldiv')
TD <- read.delim('~/Labwork/Rwork/TajimaD',header=F,sep=" ")
mat <- matrix(unlist(TD),ncol=2,byrow=F)
dTD <- data.frame(mat[,1],as.numeric(mat[,2]),'Insulin')
names(dTD) <- c('Gene','TajimaD','Set')
TDMK <- merge(SEMK,dTD,all=TRUE)

tar <- read.delim('~/Labwork/Rwork/targets1',header=F)
names(tar) <- 'Gene'
tdmk <- merge(TDMK,tar)
tdmk$Set <- 'Canonical'
for (gene in tar$Gene){TDMK <- subset(TDMK,Gene!=gene)}
TDMK <- rbind(TDMK,tdmk)

#Fold in the simulans data: 
SEMK <- read.table('~/Labwork/Rwork/SEMKproc.simfulldiv')
SEMK <- SEMK[,c(1,5:7,10)]
names(SEMK) <- c('Gene','Alpha.sim','Ds.sim','Dn.sim','pvalsout.sim')
TDMK <- merge(TDMK,SEMK,all=TRUE)

#Add ratios: 
TDMK$DnDs <- TDMK$Dn/TDMK$Ds
TDMK$DnDs.sim <- TDMK$Dn.sim/TDMK$Ds.sim
TDMK$PnPs <- TDMK$Pn/TDMK$Ps

conTDMK <- TDMK

gene <- data.frame(TDMK$Gene)
names(gene) <- 'Gene'

#Groupings:
groups <- read.table('~/Labwork/Annotation/Lists/Gene.Groups',header=F,sep=" ")
groups <- data.frame(groups)
names(groups) <- c('Gene','Group')
group <- merge(gene,groups)
# Note number of rows increases here because some genes belong to TWO sets.
TDMK <- merge(TDMK,group,all=TRUE)

# GO terms:
TERMS <- read.delim('~/Labwork/Annotation/GO/short.obo',header=F,sep="\t")
names(TERMS) <- c('Term','Name','Namespace')
ALL <- read.delim('~/Labwork/Annotation/GO/geneall',header=F,sep=" ")
names(ALL) <- c('Gene','Term')
ALL <- merge(ALL,TERMS)
ALL <- merge(gene,ALL)
# Up to 7401 rows!
TDMK <- merge(TDMK,ALL,all=TRUE)

#PFAM terms:


# A treemodel using only certain terms of conservation:
library(rpart)
treemodel <- rpart(Set ~ Alpha + DnDs + DnDs.sim + PnPs + TajimaD,data=conTDMK)
plot(treemodel) 
text(treemodel, use.n=T)




