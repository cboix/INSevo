#!/usr/bin/R
PFAM <- read.table('~/Labwork/Annotation/PFAM/MERGE.all')
species <- 'yak'
EXP <- read.delim(paste('GFExp',species,sep=""),header=F,sep=" ")
EXP <- EXP[!EXP$V2 == 'EMPTY',]
EXP <- EXP[!EXP$V4 == 'NaN',]
EXP <- EXP[!EXP$V2 == 'NaN',]
EXP$V2 <- as.numeric(as.character(EXP$V2))
EXP$V3 <- as.numeric(as.character(EXP$V3))
EXP$V4 <- as.numeric(as.character(EXP$V4))
EXP$V5 <- as.numeric(as.character(EXP$V5))
EXP <- EXP[EXP$V3 > 1,]
EXP <- EXP[EXP$V5 > 1,]

E <- data.frame(EXP)
names(E) <- c('Gene','SY.p.value','SY.df','NS.p.value','NS.df','SYa.stat','NSa.stat','SYc.stat','NSc.stat')
library(reshape)
#Create a condensed matrix in this way:
E3 <- reshape(E,idvar="Gene",direction="long",varying=list(PV=c(2,4),DF=c(3,5),AS=c(6,7),CS=c(8,9),LS=c(10,11)),v.names=c("p.value","df","Area.Statistic","Curv.Statistic","Length"))
d <- c(); for(i in 1:nrow(E3)){if(E3$time[i] == 1){d <- c(d,"SY")}else{d <- c(d,"NS")}}
E3$time <- d
rownames(E3) <- 1:nrow(E3)
names(E3)[2] <- "Type"
E3$Type <- factor(E3$Type)

levels(PFAM$Gene) -> genes
count <- c();for(i in genes){sum(PFAM$Gene == i)-> s;count <- c(count,s)}
data.frame(cbind(genes,count)) -> domains
names(domains) <- c('Gene','Num.Domains')
E3 <- merge(E3,domains)


#ORDER BY PVALUE AND APPLY BENJAMINI CORRECTION:
pvals <- sort(E$SY.p.value)
SYpvalsout <- sapply(1:length(pvals),data=pvals,function(i,data){if(pvals[i] <= i/length(pvals)*.05){return(pvals[i])}else{return(1)}})
ord <- order(E$SY.p.value)
E <- E[ord,]
E <- cbind(E,SYpvalsout)

pvals <- sort(E$NS.p.value)
NSpvalsout <- sapply(1:length(pvals),data=pvals,function(i,data){if(pvals[i] <= i/length(pvals)*.05){return(pvals[i])}else{return(1)}})
ord <- order(E$NS.p.value)
E <- E[ord,]
E <- cbind(E,NSpvalsout)

write.table(E,paste('EXPtab',species,sep=""))
write.table(E3,paste('EXPtabmelt',species,sep=""))
SYisR <- E[E$SYpvalsout == 1,]

#Data from TDMK:
DTF <- read.table(paste('TDMKtab.',species,sep=""))
merge(E,DTF) -> AL
AL[AL$NSpvalsout == 1,] -> out
AL[!AL$NSpvalsout == 1,] -> sig
AvO <- with(out,1 - mean(Ds)/mean(Dn)*mean(Pn/(Ps+1)))
AvS <- with(sig,1 - mean(Ds)/mean(Dn)*mean(Pn/(Ps+1)))

1*(AL$NSpvalsout ==  1) -> issig
AL <- cbind(AL,issig)
dummy1 <- data.frame(issig=c(1,0),line=c(AvO,AvS))

library(ggplot2)
#AlphaVdfbySig.png
 ggplot(AL,aes(NS.df,Alpha,color=NS.p.value)) + geom_point(alpha=.8) + facet_grid(~issig) + geom_hline(data=dummy1,aes(yintercept=line),col='red',linetype='dashed')

#AlphavTDbySig.png
ggplot(AL,aes(Alpha,TajimaD,color=NS.p.value)) + geom_point(alpha=.8) + facet_grid(~issig)  + geom_vline(data=dummy1,aes(xintercept=line),col='red',linetype='dashed')

# FIXME IS THIS RESULT BECAUSE I based this on DIV data?

#ggplot(E,aes(Gene,SY.p.value)) + geom_point(alpha=.5) + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + labs(x="Gene")

#PLOTS OF PFAM VS. EXP (Important for validity)

EPFAM <- merge(E3,PFAM)
non<- 1*(EPFAM$Clan.Desc=='None')
cbind(EPFAM,non) -> EPFAM
factor(EPFAM$non) -> EPFAM$non

#Reduction in size:
    levels(EPFAM$Clan.Desc) -> desc
    n <- 5
    dnum <- c()
    for(i in desc){if(sum(EPFAM$Clan.Desc == i) >= n){dnum <- c(dnum,i)}} 
    redEP <- EPFAM[EPFAM$Clan.Desc %in% dnum,]
    # make sure that no gene is expressed twice!
    redEP2 <- redEP[,c(1:7,15:16)]
    redEP2 <- unique(redEP2)
    #Add a ratio score:
    redEP2[redEP2$Type == 'SY',]$Curv.Statistic -> c
    redEP2[redEP2$Type == 'NS',]$Curv.Statistic -> d
    Ratio.CS <- d/c
    redEP2[redEP2$Type == 'SY',]$Gene -> g
    dtf2 <- data.frame(cbind(as.character(g),Ratio.CS))
    names(dtf2) <- c('Gene','Ratio.CS')
    redEP2 <- merge(redEP2,dtf2)
    as.numeric(as.character(redEP2$Ratio.CS)) -> redEP2$Ratio.CS
    # Plots for this smaller set:

    with(redEP2,by(Ratio.CS,Clan.Desc,median)) -> ameans
    as.numeric(as.character(ameans)) -> ameans
    data.frame(cbind(levels(redEP2$Clan.Desc),ameans)) -> dd
    dd$ameans <- as.numeric(as.character(dd$ameans))
    names(dd) <- c('Clan.Desc','R.CS.median')
    redEP2 <- merge(redEP2,dd)

    redEP2[redEP2$Type == 'SY',]$Area.Statistic -> c
    redEP2[redEP2$Type == 'NS',]$Area.Statistic -> d
    Ratio.AS <- d/c
    redEP2[redEP2$Type == 'SY',]$Gene -> g
    dtf2 <- data.frame(cbind(as.character(g),Ratio.AS))
    names(dtf2) <- c('Gene','Ratio.AS')
    redEP2 <- merge(redEP2,dtf2)
    as.numeric(as.character(redEP2$Ratio.AS)) -> redEP2$Ratio.AS
    


#PREDICTIVE POWER OF STATISTICS:
#Curv statistic seems the most significant in distinguishing both SY and NS and No domain v. Domain!
    ggplot(EPFAM,aes(non,log(Curv.Statistic))) + geom_boxplot() + facet_wrap(~Type) + labs(x="Has no labeled domains?",y="log of d'' Statistic")
    ggplot(EPFAM,aes(non,log(Curv.Statistic))) + geom_boxplot() + facet_wrap(~df + Type)
    
    ggplot(EPFAM,aes(Num.Domains,log(Curv.Statistic))) + geom_boxplot() + facet_wrap(~Type) + labs(y="log of curv.stat",x="Number of Domains")

    #Interesting plots for domain conservation!
    ggplot(redEP2,aes(reorder(Clan.Desc,Curv.Statistic),log(Curv.Statistic))) + geom_boxplot() + facet_wrap(~Type) + coord_flip()

    ggplot(redEP2,aes(Type,log(Curv.Statistic))) + geom_boxplot() + facet_wrap(~Clan.Desc) + coord_flip() + labs(title=paste("Curve Statistic for SY against NS for PFAM clans with ",n," members or more in the gene set",sep=""),y="Log of the Curvature Statistic",x="")

    ggplot(redEP2,aes(reorder(Clan.Desc,R.CS.median),log(Ratio.CS),color=reorder(Clan.Desc,R.CS.median))) + geom_boxplot() + coord_flip() + labs(x="Clan Description",y="Log of Ratio of Curve Statistic",title=paste("Curve Statistic Ratio for",species,"for Clans with",n,"members or more")) + theme(legend.position="none") + geom_hline(yintercept=0,color="black",linetype="dashed")


    ggplot(redEP2,aes(reorder(Clan.Desc,Ratio.AS),log(Ratio.AS),color=reorder(Clan.Desc,Ratio.AS))) + geom_boxplot() + coord_flip() + labs(x="Clan Description",y="Log of Ratio of Area Statistic",title=paste("Area Statistic Ratio for",species,"for Clans with",n,"members or more")) + theme(legend.position="none") + geom_hline(yintercept=0,color="black",linetype="dashed")


    # Now we can add the TDMK data! YAY :)


