#!/usr/bin/R
# Read in length 
args <- commandArgs(TRUE)
len <- as.numeric(args[1])
prefix <- args[2]
gene <- args[3]

# Read in file
filenameSY <- paste(prefix,"/project/divergence/",gene,".sy.dist",sep="")
filenameNS <- paste(prefix,"/project/divergence/",gene,".ns.dist",sep="")

# We require at least the SY file.... otherwise it might not be very useful.
if (file.info(filenameSY)$size == 0){

out <- paste(gene,"EMPTY",sep=" ")
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
    return(c(degf,p.value))
}
# Run for each:
SYstat <- GFExp(SY[[1]],len)
NSstat <- GFExp(NS[[1]],len)
SYdf <- SYstat[1] 
SYp.value <- SYstat[2]
NSdf <- NSstat[1]
NSp.value <- NSstat[2]

#METHOD II - The area:
patharea <- function(i,x,p){
    diff <- p[i] - x[i]
    area <- 0.5*((x[i] - x[i-1])^2 + (diff)*abs(diff))
    return(abs(area))
}

path <- function(path){
    # Normalize and create a similar, model, vector:
    l <- length(path)
    step <- 1/l
    path <- path/sum(path)
    #TODO START HALFWAY!
    steps <- rep.int(step,l)
    model <- c(step/2,rep.int(step,l-1),step/2)

    x <- c(0,cumsum(path))
    m <- c(cumsum(model))
    p <- c(0,cumsum(steps))
    parea <- sum(sapply(2:(l+1),x=x,p=p,patharea))
    marea <- sum(sapply(2:(l+1),x=m,p=p,patharea))

    return(c(parea,marea,parea-marea))
}

    
# Print out results:
out <- paste(gene,SYp.value,SYdf,NSp.value,NSdf)
cat(out)
cat("\n")
}

#SIMULATION for path method
if (FALSE){
    N <- 5000
	BernouilliProcess <- function(k,lambda){
	    n <- c()
	    while(sum(n) < k){
	        if(lambda > runif(1)){
	            n <- c(n,1)
	        } else { n <- c(n,0) }
	    }
	    return(n)
	}
    # Over 1500 nucleotides (similar numbers to Akt1):
    simpath <- function(length,occ){
        #TODO SIMULATE BY POISSON PROCESS to keep length
	    #gene <- rpois(length,occ/length)
        gene <- BernouilliProcess(occ,occ/length)
	    idx <- (1:length(gene))*(gene != 0)
	    idx <- idx[(idx != 0)]
	    dist <- c(idx[1],diff(idx))
	    return(path(dist))
    }

    sim <- sapply(rep.int(1500,N),occ=15,simpath)
    dtf <- data.frame(cbind(t(sim),15))
    names(dtf) <- c('Path.area','Model.area','Test.Statistic','Occurrences')

for(o in c(30,60,120,240,480)){
    sim <- sapply(rep.int(1500,N),occ=o,simpath)
    sim <- data.frame(cbind(t(sim),o))
    names(sim) <- c('Path.area','Model.area','Test.Statistic','Occurrences')
    dtf <- rbind(dtf,sim)

}

    library(ggplot2)
    ggplot(dtf,aes(Test.Statistic)) + geom_histogram(binwidth=.02) + facet_wrap(~Occurrences)

}




# Return position of k events w/ poisson process w/ fixed lambda:
PoisProcess <- function(k,lambda){
    n <- 0
    t < p
    while(n < k){
        

    }

}

