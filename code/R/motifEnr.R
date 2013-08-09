#!/usr/bin/R
# Create at a level of 0.1:
#   awk '$3 > 0.3 {print $4,$5,"Insulin"}' InsulinMotifs | sort | uniq -c | sort -rn | awk '{print $1,$2,$3,$4}'> IM
#   awk '$3 > 0.3 {print $4,$5,"Control"}' InsulinMotifsControl | sort | uniq -c | sort -rn | awk '{print $1,$2,$3,$4}' >> IM

IM <- read.delim('IM',header=F,sep=" ")
names(IM) <- c('Count','Gene','Disrupt','Set')
sums <- with(IM,by(Count,Set,sum))
sums <- data.frame(matrix(c('Control','Insulin',as.numeric(as.character(sums))),ncol=2))
names(sums) <- c('Set','Total')
IM <- merge(IM,sums)
IM$Total <- as.numeric(as.character(IM$Total))
IM$Frac <- IM$Count/IM$Total

data <- data.frame(Gene='Something',X1=0,X2=0,X3=0,X4=0);for(name in levels(IM$Gene)){b <- IM[IM$Gene == name,]$Frac;a <- data.frame(gene=name,b[1],b[2],b[3],b[4]);if(length(a) == 5){names(a) <- c('Gene','X1','X2','X3','X4');data<- rbind(data,a)}};data<- data[-1,]
data$fold <- data$X1/data$X3
data <- na.omit(data)

library(ggplot2)
ggplot(data,aes(Gene,fold,fill=log(X1))) + geom_bar() + coord_flip()

a <- data$fold > 1.5
a <- na.omit(a *1:length(a))
data[a,]


IM <- read.table('InsulinMotifs.tab')
tar <- read.delim('~/Labwork/Rwork/targets1',header=F)
names(tar) <- 'Gene'
tar$Set <- 'Canonical'
tIM <- merge(tar,IM)
for (gene in tar$Gene){IM <- subset(IM,Gene !=gene)}
IM$Set <- 'RNAseq'
IM <- rbind(IM,tIM)
IM[IM$Set == 'Canonical',] -> canon
IM[IM$Set == 'RNAseq',] -> rnaseq
csums <- colSums(na.omit(canon[,2:105]))/nrow(canon)
rsums <- colSums(na.omit(rnaseq[,2:105]))/nrow(rnaseq)
sums <- rbind(csums,rsums)
sums <- data.frame(sums,Set=c('Canonical','RNA-seq'))
library(reshape)
sm <- melt(sums,id.vars='Set')

library(ggplot2)
ggplot(sm,aes(variable,value,fill=Set)) + geom_bar(position='dodge') + coord_flip()

#FISHERS TESTS pairing TWO TYPES:
type1 = 'Canonical'
type2 = 'RNAseq'
IM[IM$Set == type1,] -> canon
IM[IM$Set == type2,] -> rnaseq
cfish <- rbind(colSums(na.omit(canon[,2:(ncol(canon)-2)])),nrow(canon))
rfish <- rbind(colSums(na.omit(rnaseq[,2:(ncol(rnaseq)-2)])),nrow(rnaseq))
fish <- rbind(cfish,rfish)
fish <- data.frame(fish,Set=c(type1,'C',type2,'R'))

dtf <- data.frame(Motif="",P.value=0);for(i in 1:(ncol(fish)-1)){ft <- fisher.test(round(matrix(fish[,i],ncol=2)));dtf <- rbind(dtf,data.frame(Motif=names(fish)[i],P.value=ft$p.value))}

mots <- unique(substring(names(fish[-ncol(fish)]),1,nchar(names(fish[-ncol(fish)]))-2))
fish2 <- data.frame(None=c(0,0,0,0)); for(name in mots){ns <- fish[paste(name,'NS',sep="")];sy <- fish[paste(name,'SY',sep="")];fish2 <- cbind(fish2,c(ns[1,],sy[1,],ns[3,],sy[3,])); names(fish2)[ncol(fish2)] <- name;}

dtf2 <- data.frame(Motif="",P.value=0);for(i in 1:(ncol(fish2)-1)){ft <- fisher.test(round(matrix(fish2[,i],ncol=2)));dtf2 <- rbind(dtf2,data.frame(Motif=names(fish2)[i],P.value=ft$p.value))}

fish3 <- data.frame(None=c(0,0,0,0)); for(name in mots){ns <- fish[paste(name,'NS',sep="")];sy <- fish[paste(name,'SY',sep="")];fish3 <- cbind(fish3,c(ns[1,]/ns[2,],sy[1,]/sy[2,],ns[3,]/ns[4,],sy[3,]/sy[4,])); names(fish3)[ncol(fish3)] <- name;}
fish3 <- data.frame(fish3,Type=c('NS','SY','NS','SY'),Set=c(type1,type1,type2,type2))

fish4 <- melt(fish3[,-1],idvars=c('Type','Set'))
library(ggplot2)
f4var <- fish4[fish4$value > .5,]$variable
fish4 <- fish4[fish4$variable %in% f4var,]

ggplot(fish4,aes(reorder(variable,value),value,fill=Set)) + geom_bar(position='dodge') + coord_flip() + facet_wrap(~Type) + labs(x="Transcription Factor",y="Average Occurrence per Gene",title=paste(type1,"vs.",type2))
#ggsave('EnMotifsCanonRNAseq.png',width=8,height=6)

fish3[2,]/fish3[4,] -> plt
melt(plt[,-1],id.vars=c('Type','Set')) -> nan
dev.new()
ggplot(nan,aes(reorder(variable,value),value-1)) + geom_bar() + coord_flip() + labs(x="Transcription Factor",y=paste("Ratio of",type1,"over",type2,"Occurrence per Gene - 1"))

# ----------------------In a totally new R session---------------------------------------------
# Integrate Control into pipeline:
IM <- read.table('InsulinMotifs.tab')
IC <- read.table('InsulinMotifsControl.tab')
tar <- read.delim('~/Labwork/Rwork/targets1',header=F)
names(tar) <- 'Gene'
tar$Set <- 'Canonical'
tIM <- merge(tar,IM)
for (gene in tar$Gene){IM <- subset(IM,Gene !=gene)}
IM$Set <- 'RNAseq'
IM <- rbind(IM,tIM)
IC <- data.frame(IC,Set='Control')

# Common colums -- in this case IM is a strict super set of IC
cnames <- names(IM)[names(IM) %in% names(IC)]
IC <- subset(IC,select=cnames)
IM <- subset(IM,select=cnames)
IM <- rbind(IM,IC)

# Look at sums over different datasets:
IM[IM$Set == 'Canonical',] -> canon
IM[IM$Set == 'RNAseq',] -> rnaseq
IM[IM$Set == 'Control',] -> cont 
len <- ncol(IM) - 2
csums <- colSums(na.omit(canon[,2:len]))/nrow(canon)
rsums <- colSums(na.omit(rnaseq[,2:len]))/nrow(rnaseq)
tsums <- colSums(na.omit(cont[,2:len]))/nrow(cont)
sums <- rbind(csums,rsums,tsums)
sums <- data.frame(sums,Set=c('Canonical','RNA-seq','Control'))

library(reshape)
sm <- melt(sums,id.vars='Set')

library(ggplot2)
sm2var <- sm[sm$value > .5,]$variable
sm2 <- sm[sm$variable %in% sm2var,]
ggplot(sm2,aes(reorder(variable,value),value,fill=Set)) + geom_bar(position='dodge') + coord_flip()
IMstor <- IM


# Take genes that have high expression:
t800 <- read.delim('~/Labwork/project/smalltar',header=F,sep=" ")
names(t800) <- c('Gene','Rank')
t800 <- t800[order(t800$Rank),]
t50 <- head(t800,50)
t200 <- head(t800,200)
t500 <- head(t800,500)
t0 <- merge(IM,data.frame(Gene=t50$Gene))
t2 <- merge(IM,data.frame(Gene=t200$Gene))
t5 <- merge(IM,data.frame(Gene=t500$Gene))
t0$Set <- 'Top50'
t2$Set <- 'Top200'
t5$Set <- 'Top500'
IM <- rbind(IM,t0,t2,t5)
IM$Set <- factor(IM$Set)

# Look at sums over different datasets:
IM[IM$Set == 'Top200',] -> c2
IM[IM$Set == 'Top500',] -> c5
len <- ncol(IM) - 2
sums2 <- colSums(na.omit(c2[,2:len]))/nrow(c2)
sums5 <- colSums(na.omit(c5[,2:len]))/nrow(c5)
SUMS <- rbind(csums,rsums,tsums,sums2,sums5)
SUMS <- data.frame(SUMS,Set=c('Canonical','RNA-seq','Control','Top200','Top500'))

SM <- melt(SUMS,id.vars='Set')

library(ggplot2)
SM2var <- SM[SM$value > .5,]$variable
SM2 <- SM[SM$variable %in% SM2var,]
ggplot(SM2,aes(reorder(variable,value),value,fill=Set)) + geom_bar(position='dodge') + coord_flip() + labs(x='Transcription Factors',y='Average Occurrence per Gene')

# CAN WE CLUSTER?
SM2var <- unique(SM2var)
nam <- names(IM)[names(IM) %in% SM2var]
nam <- c('Set',nam)
TOCL <- subset(IM,select=nam)
TOCL$Set <- factor(TOCL$Set)

top <- TOCL[TOCL$Set == 'Top200',]
can <- TOCL[TOCL$Set == 'Canonical',]
MOD <- rbind(top,can)
as.character(MOD$Set) -> MOD$Set
factor(MOD$Set) -> MOD$Set

# In case: TOCL <- IMstor[,c(-1,-(ncol(IMstor)-1)]

# Look predicting ONLY these two:
TOCL[TOCL$Set == 'RNAseq',] -> rnaseq
TOCL[TOCL$Set == 'Control',] -> cont 
TOCR <- rbind(rnaseq,cont)
TOCR$Set <- as.character(TOCR$Set)
TOCR$Set <- factor(TOCR$Set)
randomForest(Set ~.,data=TOCR,ntree=2000) -> rfmodel
#        Control RNAseq class.error
#Control     260    420   0.6176471
#RNAseq      272    558   0.3277108

# GBM:
newcol = data.frame(isRNAseq = (TOCR$Set == 'RNAseq'))
TOCR <- cbind(TOCR[-ncol(TOCR)],newcol)
gmodel <- gbm(isRNAseq ~ ., data=TOCR,n.trees=1000,interaction.depth=2,distribution="bernoulli")
summary(gmodel)

pca <- prcomp(TOCL[,-1])
summary(pca)
plot(pca)
biplot(pca)
library(e1071)
pcaout <- predict(pca)
pcaIM <- data.frame(pcaout[,1:14]);pcaIM <- cbind(Set=as.character(TOCL$Set),pcaIM)
cl <- cmeans(pcaIM[,-1],centers=3)
cl <- kmeans(pcaIM[,-1],centers=3)
library(cluster)
clusplot(pcaIM[,-1],cl$cluster, color=TRUE, shade=TRUE, labels=2, lines=0)

# along similar lines, working off of a matrix of about 200 different components, we do PCA:

