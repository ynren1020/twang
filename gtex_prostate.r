####################2019-09-06######################
##select GTEX normal prostate norm_count samples####
##unit log(norm_count+1)  ##########################
####################################################

#library(tidyverse)

args<-commandArgs(TRUE)
gtex<-read.delim(args[1],header = TRUE,stringsAsFactors = FALSE,check.names = FALSE)
pheno<-read.delim(args[2],header = FALSE,stringsAsFactors = FALSE)

#gtex<-read.delim("gtex_RSEM_Hugo_norm_count_sub.txt",header = TRUE,stringsAsFactors = FALSE,check.names = FALSE)
#pheno<-read.delim("GTEX_phenotype_prostate",header = FALSE,stringsAsFactors = FALSE)
pheno<-pheno[,1]
#gtex.prostate<-dplyr::select(gtex,one_of(c("sample",normal.join$all)))
#normal<-data.frame(all=names(gtex)[-1])
#normal.join<-full_join(normal,pheno,by=c("all"="V1"))
#normal.join<-na.omit(normal.join)

rownames(gtex)<-gtex$sample
gtex$sample<-NULL

##Transpose##
gtexT<-as.data.frame(t(gtex))
gtexT$sample<-rownames(gtexT)

##select##
gtexT.sub<-gtexT[gtexT$sample%in%pheno,]
gtexT.sub$sample<-NULL

##Transpose back to gene as rows##
gtexT.subT<-as.data.frame(t(gtexT.sub))
gtexT.subT$sample<-rownames(gtexT.subT)
gtexT.subT<-gtexT.subT[,c(ncol(gtexT.subT),1:ncol(gtexT.subT)-1)]

write.table(gtexT.subT,"gtex.prostate.RSEM_norm_count.txt",quote = FALSE,row.names = FALSE,col.names = TRUE,sep = "\t")




