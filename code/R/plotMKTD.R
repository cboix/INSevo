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

# FROM now on, by species:
# ------------------ McDonald Test ---------------------
MKplot <- function(species){
readname <- paste('SEMKproc.',species,sep="")
dtf <- read.table(readname)
MK <- read.table('MKdivsim')
SE <- read.table('SEMKdivsim')
tab <- cbind(MK,SE[,-1])
names(tab) <- c("Gene","Gstat","Chisq","p.value","Alpha","Ds","Dn","Ps","Pn")
tabfull <- tab
#REMOVE OUTLIERS by 3 Standard deviations out of NORMAL assumptions! (set the other variable fixed, otherwise we will have a RATIO distribution)

sdd <- with(tabfull,sd(Ds/Ps))
sdp <- with(tabfull,sd(Ps/Ds))
idxout <- c(with(tabfull,Ds/Ps > mean(Ds/Ps) + sdd*2)*1:nrow(tabfull),1:nrow(tabfull)*with(tabfull,Ps/Ds > mean(Ps/Ds) + sdp*2))
tab <- tabfull[-idxout,]

#Calculate average alpha value:
AvA <- 1 - mean(tab$Ds)/mean(tab$Dn)*mean(tab$Pn/(tab$Ps+1))
pvals <- sort(tab$p.value)
#APPLY CORRECTION TO PVALS (Benjamini-Hochberg)
pvalsout <- sapply(1:length(pvals),data=pvals,function(i,data){if(pvals[i] <= i/length(pvals)*.05){return(pvals[i])}else{return(1)}})
order(tab$p.value) -> ord
tab <- tab[ord,]
cbind(tab,pvalsout) -> dtf

dtf2 <- dtf[dtf$pvalsout != 1,]





#Plotting: (after we have put it into a dataframe by reading the data)
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


#---------------------BOTH TOGETHER:---------------------------------------------------
merge(dtf,dTD,by="Gene") -> aDTF
write.table(aDTF,'TDMKtab')
idx <- (aDTF$pvalsout != 1)*(aDTF$TajimaD < -2.132)*1:nrow(aDTF)
labDTF <- aDTF[idx,]
idx2 <- (aDTF$pvalsout !=1)
lab2DTF <- aDTF[idx2,]
idx3 <-(aDTF$pvalsout != 1)*(aDTF$TajimaD > -2.132)*1:nrow(aDTF)
lab3DTF <- aDTF[idx3,]

idx4 <-(aDTF$pvalsout != 1)*(aDTF$TajimaD > -2.132)*(aDTF$Alpha > 0)*1:nrow(aDTF)
lab4DTF <- aDTF[idx4,]

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

