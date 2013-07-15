#!/usr/bin/R
# Read in length 
args <- commandArgs(TRUE)
len <- args[1] 

#TODO read in file
filenameSY <- paste(args[2],"/project/divergence/",args[3],".sy.dist",sep="")
filenameNS <- paste(args[2],"/project/divergence/",args[3],".ns.dist",sep="")

# We require at least the SY file.... otherwise it might not be very useful.
if (file.info(filenameSY)$size == 0){

out <- paste(args[1],"EMPTY",sep=" ")
cat(out)
cat("\n")

} else {
    # Read in synonymous and non-synonymous separately:
    SY <- read.delim(filenameSY,header=FALSE,sep=" ")
    NS <- read.delim(filenameNS,header=FALSE,sep=" ")

# Make function only if we actually have non-empty files.
GFExp <- function(dist,len){
    # Make a model distribution:
    ld <- length(dist)
    h1 <- hist(dist,breaks = round(ld/4),plot=FALSE)
    obs <-h1$counts

    # Create model distribution:
    p <- pexp(h1$breaks,ld/len)
    exp <- sapply(2:length(p),data=p,function(x,data){return(data[x] -data[x-1])})*ld

    #Assess the model:
    stat <- sum((obs-exp)^2/exp)
    degf <- length(exp) - 1
    p.value <- 1 - pchisq(stat,degf)
    return(p.value)
}
# Run for each:
SYp.value <- GFExp(SY[[1]],len)
NSp.value <- GFExp(NS[[1]],len)

# Print out results:
out <- paste(args[3],SYp.value,NSp.value)
cat(out)
cat("\n")
}
