#!/usr/bin/R
EXP <- read.delim('GFExp',header=F,sep=" ")
EXP <- EXP[!EXP$V2 == 'EMPTY',]
EXP <- EXP[!EXP$V4 == 'NaN',]
EXP <- EXP[!EXP$V2 == 'NaN',]
EXP$V2 <- as.numeric(as.character(EXP$V2))
EXP$V3 <- as.numeric(as.character(EXP$V3))
EXP$V4 <- as.numeric(as.character(EXP$V4))
EXP$V5 <- as.numeric(as.character(EXP$V5))


E <- data.frame(EXP)
names(E) <- c('Gene','SY.p.value','SY.df','NS.p.value','NS.df')

#ORDER BY PVALUE AND APPLY BENJAMINI CORRECTION:
	pvals <- sort(E$SY.p.value)
	SYpvalsout <- sapply(1:length(pvals),data=pvals,function(i,data){if(pvals[i] <= i/length(pvals)*.05){return(pvals[i])}else{return(1)}})
	ord <- order(E$SY.p.value)
	E <- E[ord,]
	E <- cbind(E,SYpvalsout)
	
	pvals <- sort(E$NS.p.value)
	NSpvalsout <- sapply(1:length(pvals),data=pvals,function(i,data){if(pvals[i] <= i/length(pvals)*.05){return(pvals[i])}else{return(1)}})
	ord <- order(E$NS.p.value)
	E <- E[ord,]
	E <- cbind(E,NSpvalsout)

write.table(E,'EXPtab')
SYisR <- E[E$SYpvalsout == 1,]

DTF <- read.table('TDMKtab')
merge(E,DTF) -> AL
AL[AL$NSpvalsout == 1,] -> out
AL[!AL$NSpvalsout == 1,] -> sig
AvO <- with(out,1 - mean(Ds)/mean(Dn)*mean(Pn/(Ps+1)))
AvS <- with(sig,1 - mean(Ds)/mean(Dn)*mean(Pn/(Ps+1)))

1*(AL$NSpvalsout ==  1) -> issig
AL <- cbind(AL,issig)
dummy1 <- data.frame(issig=c(1,0),line=c(AvO,AvS))

library(ggplot2)
#AlphaVdfbySig.png
 ggplot(AL,aes(NS.df,Alpha,color=NS.p.value)) + geom_point(alpha=.8) + facet_grid(~issig) + geom_hline(data=dummy1,aes(yintercept=line),col='red',linetype='dashed')

#AlphavTDbySig.png
ggplot(AL,aes(Alpha,TajimaD,color=NS.p.value)) + geom_point(alpha=.8) + facet_grid(~issig)  + geom_vline(data=dummy1,aes(xintercept=line),col='red',linetype='dashed')

# FIXME IS THIS RESULT BECAUSE I based this on DIV data?

#ggplot(E,aes(Gene,SY.p.value)) + geom_point(alpha=.5) + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + labs(x="Gene")



