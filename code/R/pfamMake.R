#!/usr/bin/R
clans <- read.delim("~/Labwork/Annotation/PFAM/clans",sep="\t",header=F);
names(clans) <- c("Clan.Num","Clan.ID","Clan.Name","Clan.Desc");
genes <- read.delim("mclans",sep="\t",header=F);
names(genes) <- c("Clan.Num","Fam.Num");
pfam <- read.delim("genes.ann",sep=" ",header=F);
names(pfam) <- c("Fam.Name","Gene","Fam.Type","Fam.Num","Fam.ID");

PFAM <- merge(genes,clans);
as.character(PFAM$Clan.ID) -> PFAM$Clan.ID
as.character(PFAM$Clan.Name) -> PFAM$Clan.Name
as.character(PFAM$Clan.Desc) -> PFAM$Clan.Desc
a <- c("None","None","None","None","None")
PFAM <- rbind(PFAM,a);
write.table(PFAM."~/Labwork/Annotation/PFAM/PFAM.both")
MERGE <- data.frame(merge(pfam,PFAM))
write.table(MERGE,"~/Labwork/Annotation/PFAM/MERGE.all")
