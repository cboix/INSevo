#!/usr/bin/R
# Read in length 
args <- commandArgs(TRUE)
len <- as.numeric(args[1])
prefix <- args[2]
gene <- args[3]
species <- as.character(args[4])

# Read in file
if (species == 'sim'){
    filenameSY <- paste(prefix,"/project/divergence/",gene,".sy.dist",sep="")
    filenameNS <- paste(prefix,"/project/divergence/",gene,".ns.dist",sep="")
}
if (species == 'yak'){
    filenameSY <- paste(prefix,"/db/DGRP/yakgenes/",gene,".sy.dist",sep="")
    filenameNS <- paste(prefix,"/db/DGRP/yakgenes/",gene,".ns.dist",sep="")
}

# We require at least the SY file.... otherwise it might not be very useful.
if (file.info(filenameSY)$size == 0){

out <- paste(gene,"EMPTY",sep=" ")
cat(out)
cat("\n")

} else {
    # Read in synonymous and non-synonymous separately:
    SY <- read.delim(filenameSY,header=FALSE,sep=" ")
    NS <- read.delim(filenameNS,header=FALSE,sep=" ")
    SY <- unique(SY)
    NS <- unique(NS)
# Make function only if we actually have non-empty files.
GFExp <- function(dist,len){
    # Make a model distribution:
    ld <- length(dist)
    bk <- round(ld/4)
    if(bk < 2){bk <- 2}
    h1 <- hist(dist,breaks = round(bk),plot=FALSE)
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
    y <- c(0,cumsum(steps))
    m <- c(cumsum(model))

    parea <- sum(sapply(2:(l+1),x=x,p=y,patharea))
    marea <- sum(sapply(2:(l+1),x=m,p=y,patharea))

    #return(c(parea,marea,parea-marea))
    return(c(parea-marea))
}
SYa.stat <- path(SY[[1]])
NSa.stat <- path(NS[[1]])

curvature <- function(path){
    # Normalize and create a similar, model, vector:
    l <- length(path)
    step <- 1/l
    path <- path/sum(path)
    #TODO START HALFWAY!
    steps <- rep.int(step,l)
    model <- c(step/2,rep.int(step,l-1),step/2)

    x <- c(0,cumsum(path))
    y <- c(0,cumsum(steps))

    # Extraneous: 
    fder <- diff(y)/diff(x)
    # Second der is (fder[i] -fder[i-1])/(x[i] - x[i-1])
    sder <- sapply(2:length(fder),x=x,fder=fder,function(i,fder,x){return((fder[i] -fder[i-1])/(x[i] - x[i-1]))})
    return(mean(abs(sder)))
}

firstder <- function(path){
    # Normalize and create a similar, model, vector:
    l <- length(path)
    step <- 1/l
    path <- path/sum(path)
    #TODO START HALFWAY!
    steps <- rep.int(step,l)
    model <- c(step/2,rep.int(step,l-1),step/2)

    x <- c(0,cumsum(path))
    y <- c(0,cumsum(steps))

    # Extraneous: 
    fder <- diff(y)/diff(x)
    out <- sapply(fder,function(x){if(x < 1 && x > 0){return(1/x)}else{return(x)}})
    return(mean(out))
}

SYc.stat <- curvature(SY[[1]])
NSc.stat <- curvature(NS[[1]])
    
# Print out results:
out <- paste(gene,SYp.value,SYdf,NSp.value,NSdf,SYa.stat,NSa.stat,SYc.stat,NSc.stat,length(SY[[1]]),length(NS[[1]]))
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
    simpath <- function(length,occ,type){
        #TODO SIMULATE BY POISSON PROCESS to keep length
	    #gene <- rpois(length,occ/length)
        gene <- BernouilliProcess(occ,occ/length)
	    idx <- (1:length(gene))*(gene != 0)
	    idx <- idx[(idx != 0)]
	    dist <- c(idx[1],diff(idx))
        #curvature, or path:
        if(type == 'c'){return(curvature(dist))}
        if(type == 'a'){return(path(dist))}
        if(type == 'f'){return(firstder(dist))}
    }

    sim <- sapply(rep.int(1500,N),occ=15,type='a',simpath)
    dtf <- data.frame(cbind(t(sim),15))
    names(dtf) <- c('Path.area','Model.area','Test.Statistic','Occurrences')

for(o in c(30,60,120,240,480)){
    sim <- sapply(rep.int(1500,N),occ=o,type='a',simpath)
    sim <- data.frame(cbind(t(sim),o))
    names(sim) <- c('Path.area','Model.area','Test.Statistic','Occurrences')
    dtf <- rbind(dtf,sim)
}

#FOR curvature:
    N <- 2000
    sim <- sapply(rep.int(1500,N),occ=15,type='f',simpath)
    dtf <- data.frame(cbind(sim,15))
    names(dtf) <- c('Test.Statistic','Occurrences')

for(o in seq(8,24,2)){
    sim <- sapply(rep.int(1500,N),type='f',occ=o,simpath)
    sim <- data.frame(cbind(sim,o))
    names(sim) <- c('Test.Statistic','Occurrences')
    dtf <- rbind(dtf,sim)
}

    log(dtf$Test.Statistic) -> lTS
    dtf <- cbind(dtf,lTS)
    dtf$Occurrences <- factor(dtf$Occurrences)

    library(ggplot2)
    ggplot(dtf,aes(Occurrences,lTS)) + geom_violin()
    ggplot(dtf,aes(Occurrences,lTS)) + geom_boxplot()

    ggplot(dtf,aes(lTS,color=Occurrences)) + geom_density()

}




# Return position of k events w/ poisson process w/ fixed lambda:
PoisProcess <- function(k,lambda){
    n <- 0
    t < p
    while(n < k){
        

    }

}

