#######################2019-08-12###############################
##align cohort feature count results' colname with ty's result##
################################################################

library(dplyr)

input<-"TCGA.metainfo.txt"
meta<-read.delim(input,header = TRUE,stringsAsFactors = FALSE,sep = "\t")
meta$fileid<-paste0(meta$FILE_ID,".bam")

args <- commandArgs(TRUE)
df <- read.table(args[1],header = TRUE,check.names = FALSE,stringsAsFactors = FALSE)

#input2<-"UVM.tumor.txt"
#df <-read.table(input2,header = TRUE,check.names = FALSE,stringsAsFactors = FALSE)


##rename column names##
for (i in 1:ncol(df)){
n<-length(strsplit(names(df)[i],"/")[[1]])

names(df)[i]<-strsplit(names(df)[i],"/")[[1]][n]

}

dfT<-as.data.frame(t(as.matrix(df)))
dfT$fileid<-rownames(dfT)

##join dfT and meta by fileid##
dfT.join<-left_join(dfT,meta,by="fileid")
##col3 and col1##
dfT.join.sub<-dfT.join[,c(3,1)]
names(dfT.join.sub)<-NULL
dfT.join.sub.mat<-t(as.matrix(dfT.join.sub))
##first row as column names##
dfT.join.sub.mat.out<-as.data.frame(dfT.join.sub.mat)
names(dfT.join.sub.mat.out) <- as.character(unlist(dfT.join.sub.mat.out[1,]))
dfT.join.sub.mat.out<-dfT.join.sub.mat.out[-1,]
##factor to integer##
dfT.join.sub.mat.out[] <- as.integer(as.character((unlist(dfT.join.sub.mat.out, use.names = FALSE))))
dfT.join.sub.mat.out$GeneID<-"PB.69"
dfT.join.sub.mat.out$Length <-3460
##reorder column##
dfT.join.sub.mat.out<-dfT.join.sub.mat.out[,c(ncol(dfT.join.sub.mat.out)-1,ncol(dfT.join.sub.mat.out),1:(ncol(dfT.join.sub.mat.out)-2))]
#names(dfT.join.sub.mat.out)

##rbind with ty output##
#input3<-"UVM.txt"
#df3<-read.table(input3,header = TRUE,check.names = FALSE,stringsAsFactors = FALSE)
df3<-read.table(args[2],header = TRUE,check.names = FALSE,stringsAsFactors = FALSE)
#names(df3)
#if(dim(df3)[2]==dim(dfT.join.sub.mat.out)[2]) 
df4<-bind_rows(dfT.join.sub.mat.out,df3)

#calculate rpkm##
calrpkm<-function(x,y){x*10^9/sum(x)/y} #x is the column of read counts, y is the length,sum(x) is total counts

for (i in 3:ncol(df4)){
  df4[,i]<-calrpkm(df4[,i],df4[,2])
}

#write.table(df4,paste0(args[1],"rpkm.txt"),sep="\t",quote = FALSE,row.names = FALSE,col.names = TRUE)
write.table(df4,paste0(strsplit(args[1],'[.]')[[1]][1],"_",strsplit(args[1],'[.]')[[1]][2],".rpkm.txt"),sep="\t",quote = FALSE,row.names = FALSE,col.names = TRUE)
#write.table(df4,"UVM.tumor.rpkm.txt",sep="\t",quote = FALSE,row.names = FALSE,col.names = TRUE)
