#!/usr/bin/R
MK <- read.table('MKdivsim')
SE <- read.table('SEMKdivsim')

tab <- cbind(MK,SE[,-1])
names(tab) <- c("Gene","Gstat","Chisq","p.value","Alpha","Ds","Dn","Ps","Pn")

pvals <- sort(tab$p.value)
#APPLY CORRECTION TO PVALS (Benjamini-Hochberg)
pvalsout <- sapply(1:length(pvals),data=pvals,function(i,data){if(pvals[i] <= i/length(pvals)*.05){return(pvals[i])}else{return(1)}})
order(tab$p.value) -> ord
tab <- tab[ord,]
cbind(tab,pvalsout) -> dtf

dtf2 <- dtf[dtf$pvalsout != 1]

#Plotting: (after we have put it into a dataframe by reading the data)
library(ggplot2)
png('MKall.png')
g1 <- ggplot(dtf,aes(x=reorder(Gene,Alpha),Alpha,color=p.value)) 
g1 + geom_point() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) 
dev.off()
g1 <- ggplot(dtf,aes(Gene,Alpha,color=p.value)) 
g1 + geom_point() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) 

png('MKbypval.png')
g2 <- ggplot(dtf2,aes(x=reorder(Gene,Alpha),Alpha,color=pvalsout)) 
g2 + geom_point() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + labs(x="Gene")
dev.off()

#g1 <- ggplot(dtf,aes(x=X1,y=X2))
#g1 + geom_point(pch=1) + geom_hline(yintercept=c(0,-1.539,-1.765,-2.132,-2.462),color=c('black','red','red','red','red'),linetype='dotted') + theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + labs(x="Gene",y="Tajima's D")  

