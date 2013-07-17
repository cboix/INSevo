#!/usr/bin/R
# Most of the code is co-opted from various sources.
# ------------------ McDonald Test ---------------------
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

#WRITE the table:
write.table(dtf,'SEMKprocessed')

# -------------------- Tajima's D Test---------------------------
TD <- read.delim('TajimaD',header=F,sep=" ")
mat <- matrix(unlist(TD),ncol=2,byrow=F)
dTD <- data.frame(mat[,1],as.numeric(mat[,2]))
names(dTD) <- c('Gene','TajimaD')

#---------------------BOTH TOGETHER:-----------------------------
merge(dtf,dTD,by="Gene") -> aDTF
write.table(aDTF,'TDMKtab')
