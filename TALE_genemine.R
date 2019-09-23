###########################2019-09-22##############################
##Find genes similar to NEFH, tumor (low), normal(high),  #########
###################################################################

library(tidyverse)

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

##20 genes##
tale.expr.subT.sub<-tale.expr.subT[,1:20]
tale.expr.subT.sub$sample<-rownames(tale.expr.subT.sub)
##boxplot##
tale.join<-full_join(tale.expr.subT.sub,tale.surv.sub.sub,by=c("sample"="Sample ID"))
tale.join$Type[is.na(tale.join$Type)]<-"Normal"
tale.join<-na.omit(tale.join) #185
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
n<-2
primary.met.test<-list()
normal.met.test<-list()
normal.primary.test<-list()
df.out<-vector(mode = "list", length = 2)

for (i in 1:2){
primary.met.test[[i]]<-t.test(tale.join[tale.join$Type=="PRIMARY",i],tale.join[tale.join$Type=="MET",i])
normal.met.test[[i]]<-t.test(tale.join[tale.join$Type=="Normal",i],tale.join[tale.join$Type=="MET",i])
normal.primary.test[[i]]<-t.test(tale.join[tale.join$Type=="Normal",i],tale.join[tale.join$Type=="PRIMARY",i])
df.out[[i]]<-data.frame("primary.vs.met"=primary.met.test[[i]]$p.value,"normal.vs.met"=normal.met.test[[i]]$p.value,"normal.vs.primary"=normal.primary.test[[i]]$p.value)
df.out[[i]]$gene<-names(tale.join)[i]
#df.out<-rbind(df.out[i])
#return(df.out)
}

df <- do.call("rbind", df.out)

df$primary.vs.met.p.adj<-p.adjust(df$primary.vs.met, method = "fdr", n = 2) #BH == fdr
df$normal.vs.met.p.adj<-p.adjust(df$normal.vs.met,method = "fdr", n = 2)
df$normal.vs.primary.p.adj<-p.adjust(df$normal.vs.primary,method = "fdr", n = 2)

write.table(df,"TALE.gene.ttest.pvalue.txt",quote = FALSE,sep = "\t",col.names = TRUE,row.names = FALSE)

df.sub<-df[df$primary.vs.met.p.adj<=0.05,]
write.table(df.sub,"TALE.gene.ttest.pvalue.sig.txt",quote = FALSE,sep = "\t",col.names = TRUE,row.names = FALSE)

##df for pvalue output##

#names(primary.met.test)
#[1] "statistic"   "parameter"   "p.value"     "conf.int"    "estimate"    "null.value"  "alternative" "method"     
#[9] "data.name" 

p <- ggboxplot(tale.join.box, x = "group", y = "NEFH",
               color = "group", palette =c("#00AFBB", "#E7B800","#FC4E07"),
               add = "jitter", 
               xlab=FALSE,
               ylab="NEFH expressoin",
               width = 0.3)
p


# Add p-values comparing groups
my_comparisons <- list( c("MET","Localized"), c("MET", "Normal"), c("Localized", "Normal") )
p1<-p + stat_compare_means(comparisons = my_comparisons,size=6)+ # Add pairwise comparisons p-value
  stat_compare_means(label.y = 15,size=6) +
  rotate_x_text(90)

##change tick mark labels##
p2<-p1+scale_x_discrete(labels=c("Normal" = "Normal\n(N=29)", "Localized" = "Localized\n(N=131)",
                                 "MET" = "Metastatis\n(N=19)")) +
  theme(text = element_text(size=20)) 
p2



ggsave("TALE_NEFH.boxplot.pval.updated.pdf",width = 8,height = 8, dpi = 300)







tale.expr.sub<-tale.expr[tale.expr$GeneSymbol=="NEFH",3:ncol(tale.expr)]
tale.expr.sub<-na.omit(tale.expr.sub)
tale.expr.subT<-as.data.frame(t(tale.expr.sub))
tale.expr.subT$sample<-rownames(tale.expr.subT)
tale.expr.subT$NEFH<-tale.expr.subT$`17159`
tale.expr.subT$`17159`<-NULL

##median group##
for (i in 1:nrow(tale.expr.subT)){
  tale.expr.subT$group[i]<-ifelse(tale.expr.subT$NEFH[i]>=median(tale.expr.subT$NEFH),"high","low")
}

##join##
tale.surv.sub<-tale.surv.sub[,1:39]
tale.join<-full_join(tale.expr.subT,tale.surv.sub,by=c("sample"="Sample ID"))
tale.join<-tale.join[!is.na(tale.join$Type),]
for (i in 1:nrow(tale.join)){
  tale.join$status[i]<-ifelse(tale.join$Event[i]=="NO",0,1)
}
##BCR##
for (i in 1:nrow(tale.join)){
  tale.join$status2[i]<-ifelse(tale.join$BCR_Event[i]=="NO",0,1)
}

#tale.join<-tale.join[,c("SurvTime","status","group")]
#tale.join<-na.omit(tale.join)
##only primary tumor##
tale.join<-tale.join[tale.join$Type=="PRIMARY",]
#surv##
fit1<-survfit(Surv(BCR_FreeTime,status2)~group,data=tale.join)
res<-ggsurvplot(fit1,data=tale.join,
                xlab = "Days",
                ylab = "BCR_Free Survival Probability (%)",
                conf.int=TRUE,
                pval = TRUE,
                fun="pct",
                risk.table = TRUE,
                size=1,
                linetype = "strata"
)
ggsave("tale.NEFH.median.BCR_Free.PRIMARY.pdf", plot = print(res), width = 8, height = 8, dpi = 500)

