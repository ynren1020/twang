##############2019-08-05###############
##prad normal sample filter condition##
#######################################

library("edgeR")
library("sva")
#library("R.utils")
library(tidyverse)
library(purrr)
library(limma)

input2<-"PRAD.normals.counts.txt"
normal<-read.delim(input2, stringsAsFactors = FALSE, check.names=FALSE);

#group<-as.factor(rep(2,52))
#row.names(normal)<-normal$Geneid
data <- normal[,-c(1,2)];

#y <- DGEList(counts=data, group=group);
tmm <- cpm(data);
#tmm$geneid<-row.names(tmm)
cpm.ig<-as.vector(tmm[59563,]) #IGF1RAS1 #mean 0.02033955
#sum(cpm.ig==0) 29

keep <- rowSums(cpm(normal[,3:ncol(normal)])==0) <=40 #39032
keep <- rowSums(cpm(normal[,3:ncol(normal)])==0) <=50 #50092 (use this one,loose)
keep <- rowSums(cpm(normal[,3:ncol(normal)])==0) <=52 #61484
table(keep);
normal.sub<-normal[keep,]

write.table(normal.sub,"PRAD.normals.counts.filtered.loose.txt",quote = FALSE,col.names = TRUE,row.names = FALSE,sep = "\t")
