#!/usr/bin/R
#Open file:
#fileConn<-file("output.txt")
#writeLines(c("#Tajima's D test results:"),fileConn)

#Take arguments
args <- commandArgs(TRUE)
#Arguments should be 1 argument, the name of the gene:
filename <- paste("~/db/DGRP/fullgenes/",args[1],"_SNP",sep="")

if (file.info(filename)$size == 0){

out <- paste(args[1],"EMPTY",sep=" ")
cat(out)
cat("\n")

} else {
tab <- read.delim(filename,header=FALSE,sep=" ")

pairwise <- sapply(seq(10,416,2),data=tab,function(x,data){out=0;for(i in seq(x+2,418,2)){out=out+sum(abs(data[,x]-data[,i]))}; return(out)})
#Variables to be used:
    n <- length(seq(10,418,2))
    k <- sum(pairwise)/(n*(n-1))
    S <- nrow(tab)
    a1 <- sum(sapply(1:(n-1),function(x){return(1/x)}))
    a2 <- sum(sapply(1:(n-1),function(x){return(1/(x^2))}))
    b1 <- (n + 1)/(3*(n-1))
    b2 <- 2*(n^2 + n + 3)/(9*n*(n-1))
    c1 <- b1 - 1/a1
    c2 <- b2 - (n+2)/(a1*n) + a2/(a1^2)
    e1 <- c1/a1
    e2 <- c2/(a1^2 + a2)

#The test statistics:
    d <- k - S/a1
    Vd <- e1*S + e2*S*(S-1)

    D <- d/sqrt(Vd)

out <- paste(args[1],D,sep=" ")
cat(out)
cat("\n")
}

# NOTE if Tajima's D is above 0, population may have suffered a recent bottle neck or there might be OVERdominant selection.
# If Tajima's D is below 0, population size may be increasing (due to high #differences but not in large quantities so high S/a1 but not k) OR there is purifying selection at the locus

#Plotting: (after we have put it into a dataframe by reading the data)
#library(ggplot2)
#g1 <- ggplot(dtf,aes(x=X1,y=X2))
#g1 + geom_point(pch=1) + geom_hline(yintercept=c(0,-1.539,-1.765,-2.132,-2.462),color=c('black','red','red','red','red'),linetype='dotted') + theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + labs(x="Gene",y="Tajima's D")  

