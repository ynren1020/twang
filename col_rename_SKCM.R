#######################2019-08-12######################################################################################
##align cohort feature count results' colname with ty's result##
##Tricks learned from this script: When do data frame transpose, 
##the important thing to remember is that we should make sure the cells (or each element in the df) 
##are numeric values, usually, in data frame, we also have some other variables as factors or characters,
##remove it or transform it to rownames, or keep it separately. In our case here, we transform the character 
##column as rownames, and then remove that column, after this procedure, all cells are numeric, and then transform 
##the df to matrix by as.matrix, and then transpose it by t(),if we want to transform matrix back to df, 
##we can do as.data.frame(), and now the previous rownames of old df becomes the colnames of new df.
###############################################################################################################################

library(dplyr)

input<-"TCGA.metainfo.txt"
meta<-read.delim(input,header = TRUE,stringsAsFactors = FALSE,sep = "\t")
meta$fileid<-paste0(meta$FILE_ID,".bam")

##Rscript in linux##
args <- commandArgs(TRUE)
df <- read.table(args[1],header = TRUE,check.names = FALSE,stringsAsFactors = FALSE)
##Rstudio##
input2<-"SKCM.normal.txt"
df <-read.table(input2,header = TRUE,check.names = FALSE,stringsAsFactors = FALSE)
##BRCA first and second column are the same, remove one duplicate##
#df<-df[,-1]

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

##IMPORTANT LESSONS HERE:DATA FRAME TRANSPOSE SHOULD BE DONE WITH NUMERIC CELLS, OTHER CELLS REMOVE OR USE AS ROWNAMES, LIKE HERE##
##CHOOSE COLUMN TO MAKE IT AS ROWNAMES, MAKE NUMERIC DF AS MATRIX AND THEN DO TRANSPOSE, NOW THE COLNAME IS PREVIOUS ROWNAME
rownames(dfT.join.sub)<-dfT.join.sub$ALIQUOT_BARCODE
dfT.join.sub$ALIQUOT_BARCODE<-NULL
dfT.join.sub.mat<-t(as.matrix(dfT.join.sub))
##first row as column names##
dfT.join.sub.mat.out<-as.data.frame(dfT.join.sub.mat)
dfT.join.sub.mat.out$GeneID<-"PB.69"
dfT.join.sub.mat.out$Length <-3460
##reorder column##
dfT.join.sub.mat.out<-dfT.join.sub.mat.out[,c(ncol(dfT.join.sub.mat.out)-1,ncol(dfT.join.sub.mat.out),1:(ncol(dfT.join.sub.mat.out)-2))]
##rbind with ty output##
##Rstudio##
input3<-"SKCM.txt"
df3<-read.table(input3,header = TRUE,check.names = FALSE,stringsAsFactors = FALSE)
##Rscript on Linux##
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
write.table(df4,"BRCA_tumor.rpkm.txt",sep="\t",quote = FALSE,row.names = FALSE,col.names = TRUE)
write.table(df4,"SKCM_normal.rpkm.txt",sep="\t",quote = FALSE,row.names = FALSE,col.names = TRUE)
