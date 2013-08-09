#!/usr/bin/R
# Two folders, ~/Labwork/siggenes hAandlTD/ and lAandlTD/
num <- 2
folders <- c('~/Labwork/siggenes/yak/hAandlTD/','~/Labwork/siggenes/yak/lAandlTD/')
folders <- c('~/Labwork/siggenes/yak/Cf2/','~/Labwork/siggenes/yak/Foxo/')
folders2<- c('~/Labwork/siggenes/DGRP/Cf2/','~/Labwork/siggenes/DGRP/Foxo/')

folders2 <- c('~/Labwork/siggenes/DGRP/hATD/','~/Labwork/siggenes/DGRP/lATD/')
pngnames <- c("highAlpha.DnDs.yak","lowAlpha.DnDs.yak")
titles <- c("Genes with high alpha and significant Tajima's D in Yakuba","Genes with low Alpha and significant Tajima's D in Yakuba")

library(ggplot2)
for(num in 1:2) {
    folder <- folders[num]
    folder2 <- folders2[num]
    names <- list.files(folder)
    names2 <- list.files(folder2)
    genes <- sub('.out','',names)
    
    a <- c();for(j in 1:length(names2)){a <- c(a,list(read.delim(paste(folder2,names2[j],sep=""),header=F,sep=" ")))}
    sums <- c();for(j in 1:length(a)){b <- data.frame(a[[j]]);sums <- c(sums,sum(b$V6=='SY'),sum(b$V6=='NS'))}
    sums <- matrix(sums,ncol=2,byrow=T)

    DGRP <- data.frame(cbind(genes,sums))
    names(DGRP) <- c('Gene','SY','NS')

    a <- c();for(j in 1:length(names)){a <- c(a,list(read.delim(paste(folder,names[j],sep=""),header=F,sep=" ")))}
    dtf <- c();for(i in 1:length(a)){tmp <- cbind(a[[i]],genes[i]); dtf <- rbind(dtf,tmp)};names(dtf) <- c('Melanogaster','Yakuba','Mutations','Codon','Position','Type','Gene');
    out <- merge(dtf,DGRP)
    
    png(paste(pngnames[num],".png",sep=""),width=1200,height=1200)
    ggplot(dtf,aes(Codon,color=Type)) + geom_histogram(binwidth=50,alpha=.5,aes(fill=Type),position='stack') + labs(x='Nucleotide Number',y='Count',title=titles[num]) + facet_wrap(~Gene,scale='free_x')
    dev.off()

    png(paste(pngnames[num],".lab.png",sep=""),width=1200,height=1200)
    ggplot(out,aes(Codon,color=Type)) + geom_histogram(binwidth=50,alpha=.5,aes(fill=Type),position='stack') + facet_wrap(~Gene,scales='free_x') + geom_text(aes(0,25,label=paste("P_n = ",NS,"; P_s = ",SY,";",sep="")),color='black',size=4,hjust=-1)
    dev.off()

}

