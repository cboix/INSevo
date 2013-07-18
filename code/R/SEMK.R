#!/usr/bin/R
PREFIX <- '/net/home/carlesba'
SEMK <- function(x,data){
    v <- c(unlist(data[x,2:5]))
    Ds <- v[1]
    Dn <- v[2]
    Ps <- v[3]
    Pn <- v[4] 
    alpha <- 1 - (Ds/Dn)*(Pn/Ps) 
    return(c(alpha,Ds,Dn,Ps,Pn))
}

semkwrite <- function(species,prefix){

	# Read and Remove bad data points:
    readname <- paste(prefix,'/project/divergence/',species,'fulldiv',sep="")
	read.delim(readname,header=FALSE,sep=" ") -> tab

	for(i in 2:5){tab <- tab[!tab[,i] == 0,]}
	#Run MK function and write file
	out <- sapply(1:nrow(tab),data=tab,SEMK)
	out <- data.frame(t(out))
	out <- cbind(tab$V1,out)
	names(out) <- c('Gene','alpha','Ds','Dn','Ps','Pn')
	Average <- with(out, mean(Ds)/mean(Dn)*mean(Pn/(Ps + 1)))

	#Write the out data (genes/alpha values)
    tabname <- paste(prefix,'/project/SEMKdiv',species,sep="")
	write.table(out,file=tabname,sep=" ")
	
	#Create figure and send to png
	library(ggplot2)
    pngname <- paste(prefix,'/project/output/SEMK',species,'.png',sep="")
	png(pngname)
	ggplot(out,aes(Gene,alpha)) + geom_point(pch=1) + geom_hline(yintercept=Average,color='red',linetype='dotted') + theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
	dev.off()
}

semkwrite('yak',PREFIX)
semkwrite('sim',PREFIX)
