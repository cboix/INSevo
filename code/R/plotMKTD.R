#!/usr/bin/R
# -------------------- Tajima's D Test--------------------
TD <- read.delim('TajimaD',header=F,sep=" ")
mat <- matrix(unlist(TD),ncol=2,byrow=F)
dTD <- data.frame(mat[,1],as.numeric(mat[,2]))
names(dTD) <- c('Gene','TajimaD')

png('~/Labwork/Rwork/TajimaD.png',height=1200,width=1200)
gTD <- ggplot(dTD,aes(Gene,TajimaD))
gTD + geom_point(alpha=.8) + geom_hline(yintercept=c(0,-1.539,-1.765,-2.132,-2.462),color=c('black','red','red','red','red'),linetype='dashed') + theme_bw() + theme(axis.text.x = element_text(angle = 90, hjust = 1))+labs(x="Gene",y="Tajima's D")
dev.off()

#SPECIES TOGETHER:
sim <- read.table('TDMKtab.sim')
yak <- read.table('TDMKtab.yak')
sim[,c(1,5,10)] -> dsim
yak[,c(1,5,10)] -> dyak
names(dyak) <- c('Gene','yAlpha','yakpvals')
names(dsim) <- c('Gene','sAlpha','simpvals')
merge(dyak,dsim) -> BOTH
dAlpha <- BOTH$sAlpha - BOTH$yAlpha
cbind(BOTH,dAlpha) -> BOTH

#Plotting:
library(ggplot2)
png('YakvSimAlpha.png',height=1000,width=1000)
ggplot(BOTH,aes(sAlpha,yAlpha,color=yakpvals)) + geom_point(alpha=.8) + labs(x="Alpha of Simulans Divergence",y="Alpha of Yakuba Divergence") + scale_x_continuous(limits=c(-10,2)) + scale_y_continuous(limits=c(-10,2))
dev.off()

png('AlphabySpecies2.png',width=1500,heigh=1500)
ggplot(BOTH,aes(sAlpha,yAlpha,color=yakpvals)) + geom_point(alpha=.8) + labs(x="Alpha of Simulans Divergence",y="Alpha of Yakuba Divergence") + scale_x_continuous(limits=c(-10,2)) + scale_y_continuous(limits=c(-10,2)) + geom_text(data=BOTH,aes(sAlpha,yAlpha,label=Gene),size=3,vjust=1.5,hjust=.5)
dev.off()

png('AlphaDifferenceDensity.png')
ggplot(BOTH,aes(dAlpha)) + geom_density() + labs(x="Difference in Alpha value")
dev.off()

# FROM now on, by species:
# ------------------ McDonald Test ---------------------
species <- 'yak'
readname <- paste('TDMKtab.',species,sep="")
dtf <- read.table(readname)
aDTF <- dtf
dtf2 <- dtf[dtf$pvalsout != 1,]

#Different tables for plotting:
idx <- (aDTF$pvalsout != 1)*(aDTF$TajimaD < -2.132)*1:nrow(aDTF)
labDTF <- aDTF[idx,]
idx2 <- (aDTF$pvalsout !=1)
lab2DTF <- aDTF[idx2,]
idx3 <-(aDTF$pvalsout != 1)*(aDTF$TajimaD > -2.132)*1:nrow(aDTF)
lab3DTF <- aDTF[idx3,]
idx4 <-(aDTF$pvalsout != 1)*(aDTF$TajimaD > -2.132)*(aDTF$Alpha > 0)*1:nrow(aDTF)
lab4DTF <- aDTF[idx4,]


#Plotting the McDonald Kreitman values ------------------------------------------------------
library(ggplot2)
#png('MKall.png')
g1 <- ggplot(dtf,aes(x=reorder(Gene,Alpha),Alpha,color=p.value)) 
g1 + geom_point() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + geom_hline(yintercept=AvA,color='red',linetype='dashed')
#dev.off()

g1 <- ggplot(dtf,aes(Gene,Alpha,color=p.value)) 
g1 + geom_point() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + geom_hline(yintercept=AvA,color='red',linetype='dashed') 

#png('MKbypval.png')
g2 <- ggplot(dtf2,aes(x=reorder(Gene,Alpha),Alpha,color=pvalsout)) 
g2 + geom_point() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + labs(x="Gene") + geom_hline(yintercept=AvA,color='red',linetype='dashed')
#dev.off()

#Plotting the TD and MK values together -----------------------------------------------------

dev.new()
glab1 <- ggplot(aDTF,aes(Alpha,TajimaD,color=pvalsout))
glab1 + geom_point(alpha=0.8) + labs(x="MK test - Alpha",y="Tajima's D value") +  geom_hline(yintercept=c(0,-2.132),color=c('black','red'),linetype='dashed') + geom_vline(xintercept=AvA,color='green',linetype='dashed') + geom_text(data=labDTF,aes(Alpha,TajimaD,label= Gene),vjust=1.5,hjust=.1,size=3)

dev.new()
glab2 <- ggplot(aDTF,aes(Alpha,TajimaD,color=pvalsout))
glab2 + geom_point(alpha=0.8) + labs(x="MK test - Alpha",y="Tajima's D value") +  geom_hline(yintercept=c(0,-2.132),color=c('black','red'),linetype='dashed') + geom_vline(xintercept=AvA,color='green',linetype='dashed') + geom_text(data=lab2DTF,aes(Alpha,TajimaD,label= Gene),vjust=1.5,hjust=.1,size=3)

#Significant TD AND Alpha
dev.new()
gsmall <- ggplot(labDTF,aes(Alpha,TajimaD,color=pvalsout))
gsmall + geom_point(alpha=0.8) + labs(x="MK test - Alpha",y="Tajima's D value") +  geom_hline(yintercept=c(-2.132),color=c('red'),linetype='dashed') + geom_vline(xintercept=AvA,color='green',linetype='dashed') + geom_text(data=labDTF,aes(Alpha,TajimaD,label= Gene),vjust=1.5,hjust=.1,size=3)

#SIDE OF PLOT (high alpha, nonsignificant TD)
dev.new()
gside2 <- ggplot(lab4DTF,aes(Alpha,TajimaD,color=pvalsout))
gside2+ geom_point(alpha=0.8) + labs(x="MK test - Alpha",y="Tajima's D value") +  geom_hline(yintercept=c(0,-2.132),color=c('black','red'),linetype='dashed') + geom_vline(xintercept=AvA,color='green',linetype='dashed') + geom_text(data=lab4DTF,aes(Alpha,TajimaD,label= Gene),vjust=1.5,hjust=.1,size=3)
