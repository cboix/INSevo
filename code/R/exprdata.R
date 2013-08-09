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

expr <- read.delim('~/Labwork/Expdata/targets.exp',header=F,sep=" ")
names(expr) <- c('FBgn','Gene','T0','T1','T2','T3','T4','T5','T6','T7','T8','T9')
slopes <- sapply(1:nrow(expr),data=expr,function(i,data){a <- as.numeric(as.character(data[i,3:12])); b <- 1:10; dtf <- data.frame(cbind(a,b)); lm(a~b,dtf) -> model; return(model$coefficients[2]);})
expr <- cbind(expr,slopes)
mEXPR <- melt(expr,id.vars=c('Gene','FBgn','slopes'),measure.vars=c('T0','T1','T2','T3','T4','T5','T6','T7','T8','T9'))
eyak <- merge(expr,yak)
nyak <- yak[!(yak$Gene %in% eyak$Gene),]
eTD <- merge(expr,dTD)

library(ggplot2)
g1 <- ggplot(eTD,aes(TajimaD,slopes)) + geom_point(alpha=.6)
g2 <- ggplot(eyak,aes(Alpha,slopes)) + geom_point(alpha=.6)
ggsave('SlopeTD.png',g1,width=4,height=4,dpi=1200)
ggsave('SlopeAlpha.png',g2,width=4,height=4,dpi=1200)

eyak$slopes < 0 -> down
cbind(eyak,down) -> ed
names(ed)[ncol(ed)] <- "Dir"
data.frame(matrix(c(TRUE,FALSE,"Down","Up"),nrow=2,ncol=2)) -> dummy
names(dummy) <- c('Dir','Direction')
merge(ed,dummy) -> ed
ggplot(ed,aes(TajimaD,Alpha,color=slopes)) + geom_point() + facet_wrap(~Direction)
up <- ed[ed$Direction == 'Up',]
down <- ed[ed$Direction == 'Down',]
with(down,1-mean(Ds)/mean(Dn)*mean(Pn/(Ps +1)))
# [1] -0.7877475
with(up,1-mean(Ds)/mean(Dn)*mean(Pn/(Ps +1)))
# [1] -0.7475928
with(nyak,1-mean(Ds)/mean(Dn)*mean(Pn/(Ps +1)))
# [1] -0.6911052
mean(down$TajimaD)
# [1] -1.685584
mean(up$TajimaD)
# [1] -1.632076
mean(nyak$TajimaD)
# [1] -1.622571

#   To test for enrichment in either direction:
#   SEG <- data.frame(Gene=ed$Gene,Dir=ed$Direction)
#   read.table('~/Labwork/StarrSeq/InsulinMotifs.tab') -> IM
   merge(IM,SEG,all=TRUE) -> IMS
   IMS[IMS$Dir == 'Down',] -> down
   IMS[IMS$Dir == 'Up',] -> up
   colSums(na.omit(up[,2:105]))/nrow(up) -> upsums
   colSums(na.omit(down[,2:105]))/nrow(down) -> dsums
   rbind(upsums,dsums) -> sums
   data.frame(sums,Dir=c('Up','Down')) -> sums
   melt(sums,id.vars='Dir') -> sm



#----------------heatmap and binning -------------------------------
library(pheatmap)
library(gridExtra)
datamatrix=as.matrix(expr[,3:12])
m1=apply(datamatrix,1, mean)
datamatrix=datamatrix-m1

cor.ex=cor(t(expr[,3:12]))
#cor.ex[1:10, 1:30]
d.ex=1-cor.ex
h.ex=hclust(as.dist(d.ex))
plot(h.ex)

datamatrix=as.matrix(expr[,4:12]-expr[,3])
#datamatrix=as.matrix(timepart2.sig[,3:11]/timepart2.sig[,2])

p <-pheatmap(datamatrix, clustering_distance_rows='correlation',cluster_cols=F)
png('ExprHeatmap.png',width=6,height=6,res=1200,units='in')

ord <- p$tree_row$order
width <- 5
bins <- data.frame(cbind(as.character(reorder(expr$Gene,ord)),head(sort(rep.int(1:(length(ord)/width+width),width)),length(ord))))
names(bins) <- c('Gene','bin')

#-----------------Plots by binning ----------------------------------
byak <- merge(yak,bins)
factor(byak$bin) -> byak$bin
ggplot(byak,aes(bin,Alpha)) + geom_boxplot() + coord_flip()





