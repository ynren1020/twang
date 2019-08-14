################################2019-08-14#################################
################tcga CHOL normal rpkm calculation##########################
###########################################################################

library(tidyverse)

input1<-"CHOL.normal.txt"
input2<-"chol.txt"
input3<-"mapping_names.dedup.corrected.txt"

pb69<-read.delim(input1,header = TRUE,stringsAsFactors = FALSE,check.names = FALSE)
genes<-read.delim(input2,header = TRUE,stringsAsFactors = FALSE,check.names = FALSE)
mapping<-read.delim(input3,header = FALSE,stringsAsFactors = FALSE,check.names = FALSE)

##subset genes by protein coding##
mapping.pc<-mapping[mapping$V4=="protein_coding",]
genes.pc<-genes[genes$Geneid%in%mapping.pc$V1,]

##combine pb69 and protein coding##
pb69$Geneid<-"PB.69"
pb69$Length<-3460
df4rpkm<-bind_rows(pb69,genes.pc)
df4rpkm<-df4rpkm[,c(10,11,1:9)]
##calculate rpkm##
calrpkm<-function(x,y){x*10^9/sum(x)/y} #x is the column of read counts, y is the length,sum(x) is total counts

for (i in 3:ncol(df4rpkm)){
  df4rpkm[,i]<-calrpkm(df4rpkm[,i],df4rpkm[,2])
}

##choose PB.69##
df4rpkm.sub<-df4rpkm[df4rpkm$Geneid=="PB.69",]

##transpose##
rownames(df4rpkm.sub)<-df4rpkm.sub$Geneid
df4rpkm.sub$Geneid<-NULL
df4rpkm.sub$Length<-NULL
df4rpkm.subT<-as.data.frame(t(as.matrix(df4rpkm.sub)))

##add cohort and groups##
cohort<-strsplit(input1,"[.]")[[1]][1]
groups<-strsplit(input1,"[.]")[[1]][2]
df4rpkm.subT$cohort<-cohort
df4rpkm.subT$group<-groups

##output##
write.table(df4rpkm.subT,"CHOL.wide2long.normal.txt",quote = FALSE,col.names = FALSE,row.names = FALSE,sep = "\t")


