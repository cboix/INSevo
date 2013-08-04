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


