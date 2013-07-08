#!/usr/bin/R

SEMK <- function(x,data){

    v <- c(unlist(data[x,2:5]))
    Ds <- v[1]
    Dn <- v[2]
    Ps <- v[3]
    Pn <- v[4] 
    scalefactor <- (Pn + Ps)/(Dn + Ds)
    Dn <- scalefactor * Dn
    Ds <- scalefactor * Ds
    alpha <- 1 - (Ds/Dn)*(Pn/Ps) 
    return(alpha)
}

# Read and Remove bad data points:
read.delim('~/project/divergence/fulldiv',header=FALSE,sep=" ") -> tab
for(i in 2:5){tab <- tab[!tab[,i] == 0,]}

#Run MK function and write file
out <- sapply(1:nrow(tab),data=tab,SEMK)
out <- data.frame(out)
out <- cbind(tab$V1,out)

write.table(out,file="~/project/SEMKdivsim",sep=" ")

