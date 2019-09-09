############################2019-09-09#############################
####GPL6480, GPL6848, expression and annotation file ##############
###################################################################

library(GEOquery)
library(tidyverse)
library(survival)
library(survminer)

##manually download GPL6848 data##
##microarray##
expmat1<-read.delim("GSE35988-GPL6848_series_matrix.txt",header = TRUE,stringsAsFactors = FALSE,check.names = FALSE)
annot1<-read.delim("GPL6848.annot.txt",header = TRUE,stringsAsFactors = FALSE,check.names = FALSE)


##subset annot and match expmat##
annot1.sub<-annot1[,c(13,3)]
expmat1.join<-full_join(expmat1,annot1.sub,by=c("ID_REF"="Platform_SPOTID"))

##choose NEFH##
expmat1.join.sub<-expmat1.join[expmat1.join$`Gene symbol`=="NEFH",]
expmat1.join.sub<-na.omit(expmat1.join.sub)
expmat1.join.sub$ID_REF<-NULL
expmat1.join.sub$gene_name<-expmat1.join.sub$`Gene symbol`
expmat1.join.sub$`Gene symbol`<-NULL
rownames(expmat1.join.sub)<-expmat1.join.sub$gene_name
expmat1.join.sub$gene_name<-NULL

##transpose##
expmat1.join.subT<-as.data.frame(t(expmat1.join.sub))
expmat1.join.subT$sample<-rownames(expmat1.join.subT)

##assign groups##
expmat1.join.subT$group<-NULL
for (i in 1:nrow(expmat1.join.subT)){
  expmat1.join.subT$group[i]<-ifelse(str_detect(expmat1.join.subT$sample[i],"N"),"Normal",ifelse(str_detect(expmat1.join.subT$sample[i],"T"),"Localized","MET"))
}

##join expmat1.join.subT and expmat.join.subT together##
expmat1.join.subT.2p<-bind_rows(expmat1.join.subT,expmat.join.subT)

expmat1.join.subT.2p$group<-factor(expmat1.join.subT.2p$group,levels = c("Normal","Localized","MET"))
#localized 131
#MET 19
#Normal 29
p <- ggboxplot(expmat1.join.subT.2p, x = "group", y = "NEFH",
               color = "group", palette =c("#00AFBB", "#E7B800","#FC4E07"),
               add = "jitter", 
               xlab=FALSE,
               ylab="NEFH expressoin",
               width = 0.3)
p


# Add p-values comparing groups
my_comparisons <- list( c("MET","Localized"), c("MET", "Normal"), c("Localized", "Normal") )
p1<-p + stat_compare_means(comparisons = my_comparisons,size=4)+ # Add pairwise comparisons p-value
  stat_compare_means(label.y = 6.3,size=4) +
  rotate_x_text(90)

##change tick mark labels##
p2<-p1+scale_x_discrete(labels=c("Normal" = "Normal\n(N=28)", "Localized" = "Localized\n(N=59)",
                                 "MET" = "Metastatis\n(N=35)")) +
  theme(text = element_text(size=20),legend.position = "none") 
p2
ggsave("PRAD_MICH_NEFH.boxplot.pval.GPL6480_6848.pdf",width = 8,height = 8, dpi = 300)


##clinical sub##
##median group##
expmat1.join.subT.2p<-expmat1.join.subT.2p[expmat1.join.subT.2p$group=="MET",]
for (i in 1:nrow(expmat1.join.subT.2p)){
  expmat1.join.subT.2p$strata[i]<-ifelse(expmat1.join.subT.2p$NEFH[i]>=median(expmat1.join.subT.2p$NEFH),"high","low")
}


##join##
clin.sub<-clin[,c(1,9,12)]
tale.join<-full_join(expmat1.join.subT.2p,clin.sub,by=c("sample"="PATIENT_ID"))
tale.join<-na.omit(tale.join)
#tale.join<-tale.join[!is.na(tale.join$Type),]
for (i in 1:nrow(tale.join)){
  tale.join$status[i]<-ifelse(tale.join$OS_STATUS[i]=="DECEASED",1,0)
}


#tale.join<-tale.join[,c("SurvTime","status","group")]
#tale.join<-na.omit(tale.join)
##only primary tumor##
#tale.join<-tale.join[tale.join$Type=="PRIMARY",]
#surv##
fit1<-survfit(Surv(OS_MONTHS,status)~strata,data=tale.join)
res<-ggsurvplot(fit1,data=tale.join,
                xlab = "Months",
                ylab = "Overall Survival Probability (%)",
                #conf.int=TRUE,
                pval = TRUE,
                fun="pct",
                risk.table = TRUE,
                size=1,
                linetype = "strata"
)
ggsave("PRAD_UMICH.NEFH.median.Overall.met.GPL6480_6848.pdf", plot = print(res), width = 8, height = 8, dpi = 500)



