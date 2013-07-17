#!/usr/bin/R

MK <- function(x,data){

    mat <- matrix(unlist(data[x,2:5]),ncol=2,byrow=FALSE)
    chmat <- chisq.test(mat)
    gstat <- sum(chmat$observed*log(chmat$observed/chmat$expected))*2
    fi <- fisher.test(mat)
    return(c(gstat,chmat$p.value,fi$p.value))

}

# Read and Remove bad data points:
read.delim('~/project/divergence/fulldiv',header=FALSE,sep=" ") -> tab

for(i in 2:5){tab <- tab[!tab[,i] == 0,]}

#Run MK function and write file
out <- sapply(1:nrow(tab),data=tab,MK)
out <- data.frame(t(out))
out <- cbind(tab$V1,out)

write.table(out,file="~/project/MKdivsim",sep=" ")

