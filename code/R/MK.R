#!/usr/bin/R
PREFIX <- '/net/home/carlesba'
MK <- function(x,data){

    mat <- matrix(unlist(data[x,2:5]),ncol=2,byrow=FALSE)
    chmat <- chisq.test(mat)
    gstat <- sum(chmat$observed*log(chmat$observed/chmat$expected))*2
    fi <- fisher.test(mat)
    return(c(gstat,chmat$p.value,fi$p.value))

}

mkwrite <- function(species,prefix){

	# Read and Remove bad data points:
	readname <- paste(prefix,'/project/divergence/',species,'fulldiv',sep="")
	read.delim(readname,header=FALSE,sep=" ") -> tab
	
	for(i in 2:5){tab <- tab[!tab[,i] == 0,]}
	#Run MK function and write file
	out <- sapply(1:nrow(tab),data=tab,MK)
	out <- data.frame(t(out))
	out <- cbind(tab$V1,out)
	tabname<-paste(prefix,'/project/MKdiv',species,sep="")
	write.table(out,file=tabname,sep=" ")
	# --- for YAKUBA:
}

mkwrite('yak',PREFIX)
mkwrite('sim',PREFIX)
