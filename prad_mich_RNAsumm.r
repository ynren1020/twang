##############################2019-09-13#############################
##Arul nature paper's data GSE35988 (RNA microarray data)  ##########
##summarize probe value into gene expression value  #################
#####################################################################

library(GEOquery)
library(tidyverse)
library(knitr)

# If you have network access, the more typical way to do this
# would be to use this:
 #GSE35988
gse_s <- getGEO("GSE35988",GSEMatrix=FALSE,destdir = "/Users/yren/Projects/twang/paper") 
head(Meta(gse_s))

# names of all the GSM objects contained in the GSE
names(GSMList(gse_s))

# and get the first GSM object on the list
GSMList(gse_s)[[1]]

# and the names of the GPLs represented
names(GPLList(gse_s))

# Converting to BioConductor ExpressionSets and limma MALists
# GEO datasets are (unlike some of the other GEO entities), quite similar to the limma data structure MAList and to the Biobase data structure ExpressionSet. 
# Therefore, there are two functions, GDS2MA and GDS2eSet that accomplish that task.

# Getting GSE Series Matrix files as an ExpressionSet
# GEO Series are collections of related experiments. In addition to being available as SOFT format files, which are quite large, NCBI GEO has prepared a simpler format file based on tab-delimited text. The getGEO function can handle this format and will parse very large GSEs quite quickly. The data structure returned from this parsing is a list of ExpressionSets. As an example, we download and parse GSE2553.

# Note that GSEMatrix=TRUE is the default
gse <- getGEO("GSE35988",GSEMatrix=TRUE,destdir = "/Users/yren/Projects/twang/paper") 
show(gse)

GPL6480<-gse[[1]]
GPL6848<-gse[[2]]

show(pData(phenoData(GPL6480))[1:5,c(1,6,8)])
#show(pData(phenoData(GPL6480))[1:5,])

#Converting GSE to an ExpressionSet
gsmlist_6480 = Filter(function(gsm) {Meta(gsm)$platform_id=='GPL6480'},GSMList(gse_s))
length(gsmlist_6480)

gsmlist_6848 = Filter(function(gsm) {Meta(gsm)$platform_id=='GPL6848'},GSMList(gse_s))
length(gsmlist_6848)

Table(gsmlist_6480[[1]])[1:5,]
Table(gsmlist_6848[[1]])[1:5,]

# and get the column descriptions
Columns(gsmlist_6480[[1]])[1:5,]

# get the probeset ordering
probesets_6480 <- Table(GPLList(gse_s)[[1]])$ID
probesets_6848 <- Table(GPLList(gse_s)[[2]])$ID
# make the data matrix from the VALUE columns from each GSM
# being careful to match the order of the probesets in the platform
# with those in the GSMs
##FUNCTION to make probe-gene set expression##
probe_gene.exp<-function(gsmlist,probesets,n){
  data.matrix <- do.call('cbind',lapply(gsmlist,function(x) 
  {tab <- Table(x)
  mymatch <- match(probesets,tab$ID_REF)
  return(tab$VALUE[mymatch])
  }))
  data.matrix <- apply(data.matrix,2,function(x) {as.numeric(as.character(x))})
  #data.matrix <- log2(data.matrix)
  #data.matrix_6480[1:5,]
  rownames(data.matrix) <- probesets
  colnames(data.matrix) <- names(gsmlist)
  
  ##match gene and probe##
  genesets<-Table(GPLList(gse_s)[[n]])[,c("ID","GENE_SYMBOL")]
  df<-as.data.frame(data.matrix)
  df$probe<-rownames(df)
  df.join<-full_join(df,genesets,by=c("probe"="ID"))
  df.join<-select(df.join,c("probe","GENE_SYMBOL",everything()))
  ##output##
 # write.table(df_6480.join,"./prad_mich/GPL_6480_probe_gene.exp.txt",quote = FALSE,sep = "\t",col.names = TRUE,row.names = FALSE)
 return(df.join)
}

#apply function to 6840##
df_6848<-probe_gene.exp(gsmlist_6848,probesets_6848,2)
df_6480<-probe_gene.exp(gsmlist_6480,probesets_6480,1)
write.table(df_6848,"./prad_mich/GPL_6848_probe_gene.exp.txt",quote = FALSE,sep = "\t",col.names = TRUE,row.names = FALSE)
write.table(df_6480,"./prad_mich/GPL_6480_probe_gene.exp.txt",quote = FALSE,sep = "\t",col.names = TRUE,row.names = FALSE)

##summarize by gene##
df_6480<-df_6480[,-1]
df_6480<-filter(df_6480,GENE_SYMBOL!="")
df_6480<-df_6480%>%group_by(GENE_SYMBOL)%>%
  summarise_each(funs(mean(., na.rm = TRUE))) #remove NA

expmean<-function(df){
  df<-df[,-1]
  df<-filter(df,GENE_SYMBOL!="")
  df<-df%>%group_by(GENE_SYMBOL)%>%
    summarise_each(funs(mean(., na.rm = TRUE)))
  return(df)
}

df_6848.sum<-expmean(df_6848)

write.table(df_6848.sum,"./prad_mich/GPL_6848_probe_gene.exp_mean.txt",quote = FALSE,sep = "\t",col.names = TRUE,row.names = FALSE)
write.table(df_6480,"./prad_mich/GPL_6480_probe_gene.exp_mean.txt",quote = FALSE,sep = "\t",col.names = TRUE,row.names = FALSE)
