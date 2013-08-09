#FISHERS TESTS pairing TWO TYPES:
IM <- read.table('FullIM.tab')
library(reshape)
library(ggplot2)

# Change these to: Canonical, RNAseq, Top50, Top200, Top500, or Control.
type1 = 'Canonical'
type2 = 'Control'

IM[IM$Set == type1,] -> canon
IM[IM$Set == type2,] -> rnaseq
can2 <- canon[,2:(ncol(canon)-2)]
rna2 <- rnaseq[,2:(ncol(rnaseq)-2)]

cfish <- rbind(colSums(na.omit(can2)),nrow(canon))
rfish <- rbind(colSums(na.omit(rna2)),nrow(rnaseq))

scfish <- sapply(1:ncol(rna2),function(i){return(sd(rna2[,i]))})
srfish <- sapply(1:ncol(can2),function(i){return(sd(can2[,i]))})
sdtot <- data.frame(Name=names(rna2),sd=sqrt(scfish^2 + srfish^2))
sdtot$Motif <- substring(sdtot$Name,1,nchar(as.character(sdtot$Name))-2)
sdtot$Type <- substring(sdtot$Name,nchar(as.character(sdtot$Name))-1,nchar(as.character(sdtot$Name)))
sdSY <- sdtot[sdtot$Type == 'SY',]

fish <- rbind(cfish,rfish)
fish <- data.frame(fish,Set=c(type1,'C',type2,'R'))

dtf <- data.frame(Motif="",P.value=0);for(i in 1:(ncol(fish)-1)){ft <- fisher.test(round(matrix(fish[,i],ncol=2)));dtf <- rbind(dtf,data.frame(Motif=names(fish)[i],P.value=ft$p.value))}

mots <- unique(substring(names(fish[-ncol(fish)]),1,nchar(names(fish[-ncol(fish)]))-2))
fish2 <- data.frame(None=c(0,0,0,0)); for(name in mots){ns <- fish[paste(name,'NS',sep="")];sy <- fish[paste(name,'SY',sep="")];fish2 <- cbind(fish2,c(ns[1,],sy[1,],ns[3,],sy[3,])); names(fish2)[ncol(fish2)] <- name;}

dtf2 <- data.frame(Motif="",P.value=0);for(i in 1:(ncol(fish2)-1)){ft <- fisher.test(round(matrix(fish2[,i],ncol=2)));dtf2 <- rbind(dtf2,data.frame(Motif=names(fish2)[i],P.value=ft$p.value))}

fish3 <- data.frame(None=c(0,0,0,0)); for(name in mots){ns <- fish[paste(name,'NS',sep="")];sy <- fish[paste(name,'SY',sep="")];fish3 <- cbind(fish3,c(ns[1,]/ns[2,],sy[1,]/sy[2,],ns[3,]/ns[4,],sy[3,]/sy[4,])); names(fish3)[ncol(fish3)] <- name;}
fish3 <- data.frame(fish3,Type=c('NS','SY','NS','SY'),Set=c(type1,type1,type2,type2))

fish4 <- melt(fish3[,-1],idvars=c('Type','Set'))
# CHANGE THIS VALUE HERE FROM .5 if you want to see more motifs in plot:
f4var <- fish4[fish4$value > .5,]$variable
fish4 <- fish4[fish4$variable %in% f4var,]

ggplot(fish4,aes(reorder(variable,value),value,fill=Set)) + geom_bar(position='dodge') + coord_flip() + facet_wrap(~Type) + labs(x="Transcription Factor",y="Average Occurrence per Gene",title=paste(type1,"vs.",type2))

#fish4[fish4$Type == 'SY',] -> fish5
#dev.new()
#ggplot(fish5,aes(reorder(variable,value),value,fill=Set)) + geom_bar(position='dodge') + coord_flip() + labs(x="Transcription Factor",y="Average Occurrence per Gene")


# RATIOS:
fish3[,c(-(ncol(fish3)-1),-ncol(fish3))] -> fish6
nam <- names(fish3)[c(fish3[2,c(-(ncol(fish3)-1),-ncol(fish3))] > .3)]

gg_color_hue <- function(n) {
    hues = seq(15, 375, length=n+1)
    hcl(h=hues, l=65, c=100)[1:n]
}
cols = gg_color_hue(4)

f <- fish6[2,-1]/fish6[4,-1]
f <- data.frame(Motif=names(f),Value=t(f[1,]))
f[f$Motif %in% nam,] -> f2
f2[f2$X2 > 1,] -> f3
ggplot(f3,aes(reorder(Motif,X2),X2)) + geom_bar(fill=cols[3]) + coord_flip() + labs(x='Transcription Factor',y=expression(paste(frac(Canonical,Control),'  Occurrence per gene')))
dat <- f3
dat$p <- 1
# BASED ON FISHERS, not very exact:
dat$p[dat$X2 > 1.7] <- .04
dat$p[dat$X2 > 2] <- .04
dat$p[dat$X2 > 2.5] <- .0055
dat$star <- ""
dat$star[dat$p <= .1] <- "+"
dat$star[dat$p <= .05]  <- "*"
dat$star[dat$p <= .01]  <- "**"
dat$star[dat$p <= .001] <- "***"

dat <- merge(dat,sdSY)
dat2 <- dat[dat$sd < 1.7,]

ggplot(dat2,aes(reorder(Motif,X2),X2,ymin=X2-sd,ymax=X2+sd,width=.5)) + geom_bar(fill=cols[3]) + geom_errorbar(stat="identity",alpha=.5) + geom_text(aes(Motif,X2 + .1,label=star), colour="black", vjust=.5,hjust=-.7, size=10) + coord_flip() + labs(x='Transcription Factor',y=expression(paste(frac(Canonical,Control),'  Occurrence per gene')))

ggplot(dat2,aes(reorder(Motif,X2),X2,ymin=X2-sd,ymax=X2+sd,width=.5)) + geom_bar(fill='black',color='black',alpha=.5) + geom_errorbar(stat="identity",alpha=.5) + geom_text(aes(Motif,X2 + .1,label=star), colour="black", vjust=.5,hjust=-.7, size=10) + coord_flip() + labs(x='Transcription Factor',y=expression(paste(frac(Canonical,Control),'  Occurrence per gene')))




fish3[2,]/fish3[4,] -> plt
melt(plt[,-1],id.vars=c('Type','Set')) -> nan
dev.new()
ggplot(nan,aes(reorder(variable,value),value-1)) + geom_bar() + coord_flip() + labs(x="Transcription Factor",y=paste("Ratio of",type1,"over",type2,"Occurrence per Gene - 1"))

