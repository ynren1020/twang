###########################2019-09-22##############################
##Find genes similar to NEFH, tumor (low), normal(high),  #########
###################################################################

library(dplyr)
#args<-commandArgs(TRUE)

##TALE ##
input1<-"MSKCC_PCa_Clinical_Annotation.txt"
input2<-"MSKCC_PCa_mRNA_data.txt"
##survival data
tale.surv<-read.delim(input1,header = TRUE,stringsAsFactors = FALSE,check.names = FALSE)
tumor<-c("PRIMARY","MET")
tale.surv.sub<-tale.surv[tale.surv$Type%in%tumor,]
tale.surv.sub.sub<-tale.surv.sub[,1:2]

##microarray data##
tale.expr<-read.delim(input2,header = TRUE,stringsAsFactors = FALSE,check.names = FALSE) #26447 187
tale.expr$GeneID<-NULL
tale.expr.sub<-tale.expr[!duplicated(tale.expr$GeneSymbol),]
tale.expr.sub<-na.omit(tale.expr.sub) #26446
##col to rowname##
rownames(tale.expr.sub)<-tale.expr.sub$GeneSymbol
tale.expr.sub$GeneSymbol<-NULL

##transpose##
tale.expr.subT<-as.data.frame(t(tale.expr.sub))
tale.expr.subT$sample<-rownames(tale.expr.subT)
##20 genes##
#tale.expr.subT.sub<-tale.expr.subT[,1:20]
#tale.expr.subT.sub$sample<-rownames(tale.expr.subT.sub)
##boxplot##
tale.join<-full_join(tale.expr.subT,tale.surv.sub.sub,by=c("sample"="Sample ID"))
tale.join$Type[is.na(tale.join$Type)]<-"Normal"
#tale.join<-na.omit(tale.join) #185
#tale.join.box<-tale.join[1:179,1:4]

#for (i in 1:nrow(tale.join.box)){
#  tale.join.box$group[i]<-ifelse(tale.join.box$Type[i]=="PRIMARY","Localized",ifelse(tale.join.box$Type[i]=="MET","MET","Normal"))
#}
#tale.join.box$group[is.na(tale.join.box$group)]<-"Normal"
tale.join$Type<-factor(tale.join$Type,levels = c("Normal","PRIMARY","MET"))
#localized 131
#MET 19
#Normal 35

##t test or anova##
n<-ncol(tale.join)-2
#n<-2
primary.met.test<-list()
normal.met.test<-list()
normal.primary.test<-list()
df.out<-vector(mode = "list", length = n)

for (i in 1:n){
  primary.met.test[[i]]<-t.test(tale.join[tale.join$Type=="PRIMARY",i],tale.join[tale.join$Type=="MET",i],na.action=na.omit)
  normal.met.test[[i]]<-t.test(tale.join[tale.join$Type=="Normal",i],tale.join[tale.join$Type=="MET",i],na.action=na.omit)
  normal.primary.test[[i]]<-t.test(tale.join[tale.join$Type=="Normal",i],tale.join[tale.join$Type=="PRIMARY",i],na.action=na.omit)
  df.out[[i]]<-data.frame("primary.vs.met"=primary.met.test[[i]]$p.value,"normal.vs.met"=normal.met.test[[i]]$p.value,"normal.vs.primary"=normal.primary.test[[i]]$p.value)
  df.out[[i]]$gene<-names(tale.join)[i]
  #df.out<-rbind(df.out[i])
  #return(df.out)
}

df <- do.call("rbind", df.out)

df$primary.vs.met.p.adj<-p.adjust(df$primary.vs.met, method = "fdr", n = n) #BH == fdr
df$normal.vs.met.p.adj<-p.adjust(df$normal.vs.met,method = "fdr", n = n)
df$normal.vs.primary.p.adj<-p.adjust(df$normal.vs.primary,method = "fdr", n = n)

write.table(df,"TALE.gene.ttest.pvalue.txt",quote = FALSE,sep = "\t",col.names = TRUE,row.names = FALSE)

df.sub<-df[df$primary.vs.met.p.adj<=0.05,]
write.table(df.sub,"TALE.gene.ttest.pvalue.sig.txt",quote = FALSE,sep = "\t",col.names = TRUE,row.names = FALSE)
