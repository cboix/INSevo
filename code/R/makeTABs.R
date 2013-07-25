#!/usr/bin/R
#First line of processing, creating the tables.
# ------------------ McDonald Test ---------------------
MKprocess <- function(species,dTD){
    readnameMK <- paste('MKdiv',species,sep="")
    readnameSE <- paste('SEMKdiv',species,sep="")
    MK <- read.table(readnameMK)
    SE <- read.table(readnameSE)
    tab <- cbind(MK,SE[,-1])
    names(tab) <- c("Gene","Gstat","Chisq","p.value","Alpha","Ds","Dn","Ps","Pn")
    tabfull <- tab
    #REMOVE OUTLIERS by 3 Standard deviations out of NORMAL assumptions! (set the other variable fixed, otherwise we will have a RATIO distribution)
    sdd <- with(tabfull,sd(Ds/Ps))
    sdp <- with(tabfull,sd(Ps/Ds))
    idxout <- c(with(tabfull,Ds/Ps > mean(Ds/Ps) + sdd*2)*1:nrow(tabfull),1:nrow(tabfull)*with(tabfull,Ps/Ds > mean(Ps/Ds) + sdp*2))
    tab <- tabfull[-idxout,]
    
    #Calculate average alpha value:
    AvA <- 1 - mean(tab$Ds)/mean(tab$Dn)*mean(tab$Pn/(tab$Ps+1))
    pvals <- sort(tab$p.value)
    #APPLY CORRECTION TO PVALS (Benjamini-Hochberg)
    pvalsout <- sapply(1:length(pvals),data=pvals,function(i,data){if(pvals[i] <= i/length(pvals)*.05){return(pvals[i])}else{return(1)}})
    order(tab$p.value) -> ord
    tab <- tab[ord,]
    cbind(tab,pvalsout) -> dtf
    dtf2 <- dtf[dtf$pvalsout != 1,]
    #WRITE the table:
    tabname <- paste('SEMKproc.',species,sep="")
    write.table(dtf,tabname)

    #WRITE merged table:
    aDTF <- merge(dtf,dTD,by='Gene')
    tabname2 <- paste('TDMKtab.',species,sep="")
    write.table(aDTF,tabname2)
}

# -------------------- Tajima's D Test---------------------------
TD <- read.delim('TajimaD',header=F,sep=" ")
mat <- matrix(unlist(TD),ncol=2,byrow=F)
dTD <- data.frame(mat[,1],as.numeric(mat[,2]))
names(dTD) <- c('Gene','TajimaD')

MKprocess('yak',dTD)
MKprocess('sim',dTD)
