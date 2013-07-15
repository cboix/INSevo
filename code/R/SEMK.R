#!/usr/bin/R

SEMK <- function(x,data){
    v <- c(unlist(data[x,2:5]))
    Ds <- v[1]
    Dn <- v[2]
    Ps <- v[3]
    Pn <- v[4] 
    alpha <- 1 - (Ds/Dn)*(Pn/Ps) 
    return(c(alpha,Ds,Dn,Ps,Pn))
}

# Read and Remove bad data points:
read.delim('~/project/divergence/fulldiv',header=FALSE,sep=" ") -> tab
#read.delim('~/Labwork/Rwork/Divall',header=FALSE,sep=" ") -> tab
for(i in 2:5){tab <- tab[!tab[,i] == 0,]}

#Run MK function and write file
out <- sapply(1:nrow(tab),data=tab,SEMK)
out <- data.frame(t(out))
out <- cbind(tab$V1,out)
names(out) <- c('Gene','alpha','Ds','Dn','Ps','Pn')
Average <- with(out, mean(Ds)/mean(Dn)*mean(Pn/(Ps + 1)))

#Write the out data (genes/alpha values)
write.table(out,file="~/project/SEMKdivsim",sep=" ")

#Create figure and send to png
library(ggplot2)
g1 <- ggplot(out,aes(Gene,alpha)) + geom_point(pch=1) + geom_hline(yintercept=Average,color='red',linetype='dotted') + theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1))

png('~/project/output/SEMK.png')
g1
dev.off()
