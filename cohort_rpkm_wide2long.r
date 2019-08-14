#########################2019-08-13#######################
##rpkm wide to long, df prep for boxplot #################
##########################################################
#library(dplyr)
#library(tidyr)

##Rstudio##
#input<-"SKCM_normal.rpkm.txt"
#df<-read.delim(input,stringsAsFactors = FALSE,check.names = FALSE,header = TRUE)

##Rscript##
args <- commandArgs(TRUE)
df <- read.table(args[1],header = TRUE,check.names = FALSE,stringsAsFactors = FALSE)


##choose PB.69##
df.sub<-df[df$GeneID=="PB.69",]
rownames(df.sub)<-df.sub$GeneID
df.sub$GeneID<-NULL
df.sub$Length<-NULL

##transpose##
df.subT<-as.data.frame(t(as.matrix(df.sub)))

##add cohort and group##
cohort<-strsplit(input,"_")[[1]][1]
groups<-strsplit(strsplit(input,"_")[[1]][2],"[.]")[[1]][1]
df.subT$cohort<-cohort
df.subT$group<-groups
  
write.table(df.subT,paste0(cohort,".wide2long.",groups,".txt"),quote = FALSE,row.names = FALSE,col.names = FALSE,sep="\t")


  