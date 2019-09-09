####################2019-09-09####################
##NEFH expression ~ survival######################
################prad_mich#########################
##UCSF:mmc5,2018_04_15_matrix_rna_tpm.txt#########
##PROMOTE:tier.info.txt, /home/tywang/Projects/lncRNAs/PROMOTE/known_counts/results/V1/rpkm.count.txt
##TALE:MSKCC_****,xls,MSKCC**.txt#################

##install package and download data from GEO##
BiocManager::install("GEOquery")
library(GEOquery)
library(tidyverse)
library(survival)
library(survminer)

# If you have network access, the more typical way to do this
# would be to use this:
gse <- getGEO("GSE35988",GSEMatrix=FALSE)  #GSE35988
head(Meta(gse))

##manually download GPL6480 data##
##microarray##
expmat<-read.delim("GSE35988-GPL6480_series_matrix.txt",header = TRUE,stringsAsFactors = FALSE,check.names = FALSE)
annot<-read.delim("GPL6480.annot.txt",header = TRUE,stringsAsFactors = FALSE,check.names = FALSE)
clin<-read.delim("data_clinical_patient_pradmich.txt",header = TRUE,stringsAsFactors = FALSE,check.names = FALSE)

##subset annot and match expmat##
annot.sub<-annot[,c(1,3)]
expmat.join<-full_join(expmat,annot.sub,by=c("ID_REF"="ID"))

##choose NEFH##
expmat.join.sub<-expmat.join[expmat.join$`Gene symbol`=="NEFH",]
expmat.join.sub<-na.omit(expmat.join.sub)
expmat.join.sub$ID_REF<-NULL
expmat.join.sub$gene_name<-expmat.join.sub$`Gene symbol`
expmat.join.sub$`Gene symbol`<-NULL
rownames(expmat.join.sub)<-expmat.join.sub$gene_name
expmat.join.sub$gene_name<-NULL

##transpose##
expmat.join.subT<-as.data.frame(t(expmat.join.sub))
expmat.join.subT$sample<-rownames(expmat.join.subT)

##assign groups##
expmat.join.subT$group<-NULL
for (i in 1:nrow(expmat.join.subT)){
  expmat.join.subT$group[i]<-ifelse(str_detect(expmat.join.subT$sample[i],"N"),"Normal",ifelse(str_detect(expmat.join.subT$sample[i],"T"),"Localized","MET"))
  
}

##Normal 12, primary 49, Mets 27##

##boxplot##
#tale.join<-full_join(tale.expr.subT,tale.surv.sub,by=c("sample"="Sample ID"))
#tale.join.box<-tale.join[1:179,1:4]

#for (i in 1:nrow(tale.join.box)){
#  tale.join.box$group[i]<-ifelse(tale.join.box$Type[i]=="PRIMARY","Localized",ifelse(tale.join.box$Type[i]=="MET","MET","Normal"))
#}
#tale.join.box$group[is.na(tale.join.box$group)]<-"Normal"
expmat.join.subT$group<-factor(expmat.join.subT$group,levels = c("Normal","Localized","MET"))
#localized 131
#MET 19
#Normal 29
p <- ggboxplot(expmat.join.subT, x = "group", y = "NEFH",
               color = "group", palette =c("#00AFBB", "#E7B800","#FC4E07"),
               add = "jitter", 
               xlab=FALSE,
               ylab="NEFH expressoin",
               width = 0.3)
p


# Add p-values comparing groups
my_comparisons <- list( c("MET","Localized"), c("MET", "Normal"), c("Localized", "Normal") )
p1<-p + stat_compare_means(comparisons = my_comparisons,size=4,face = "bold")+ # Add pairwise comparisons p-value
  stat_compare_means(label.y = 6.3,size=4,face = "bold") +
  rotate_x_text(90)

##change tick mark labels##
p2<-p1+scale_x_discrete(labels=c("Normal" = "Normal\n(N=12)", "Localized" = "Localized\n(N=49)",
                                 "MET" = "Metastatis\n(N=27)")) +
  theme(text = element_text(size=20),legend.position = "none") 
p2



ggsave("PRAD_MICH_NEFH.boxplot.pval.pdf",width = 8,height = 8, dpi = 300)

##clinical sub##
##median group##
expmat.join.subT.sub<-expmat.join.subT[expmat.join.subT$group=="MET",]
for (i in 1:nrow(expmat.join.subT.sub)){
  expmat.join.subT.sub$strata[i]<-ifelse(expmat.join.subT.sub$NEFH[i]>=median(expmat.join.subT.sub$NEFH),"high","low")
}


##join##
clin.sub<-clin[,c(1,9,12)]
tale.join<-full_join(expmat.join.subT.sub,clin.sub,by=c("sample"="PATIENT_ID"))
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
ggsave("PRAD_UMICH.NEFH.median.Overall.met.pdf", plot = print(res), width = 8, height = 8, dpi = 500)
















input1<-"mmc5.txt"
input2<-"2018_04_15_matrix_rna_tpm.txt"
ucsf.surv<-read.delim(input1,header = TRUE,stringsAsFactors = FALSE,check.names = FALSE)
ucsf.expr<-read.delim(input2,header = TRUE,stringsAsFactors = FALSE,check.names = FALSE)

##NEFH expression##
ucsf.expr.sub<-ucsf.expr[ucsf.expr$IDENTIFIER=="NEFH",2:ncol(ucsf.expr)]
ucsf.expr.subT<-as.data.frame(t(ucsf.expr.sub))
ucsf.expr.subT$Patient.ID<-rownames(ucsf.expr.subT)
ucsf.expr.subT<-ucsf.expr.subT%>%rename("NEFH"="9534")
ucsf.expr.subT$sample<-str_replace_all(ucsf.expr.subT$Patient.ID,c("-BL"),"")
ucsf.expr.subT$sample<-str_replace_all(ucsf.expr.subT$sample,"-PRO","")
ucsf.expr.subT$sample<-str_replace_all(ucsf.expr.subT$sample,"-PRO2","")
##median group##
for (i in 1:nrow(ucsf.expr.subT)){
  ucsf.expr.subT$group[i]<-ifelse(ucsf.expr.subT$NEFH[i]>=median(ucsf.expr.subT$NEFH),"high","low")
  }
##<10 group##
for (i in 1:nrow(ucsf.expr.subT)){
  ucsf.expr.subT$group2[i]<-ifelse(ucsf.expr.subT$NEFH[i]>=10,"high","low")
}

##top25 and bottom25 group##
for (i in 1:nrow(ucsf.expr.subT)){
  ucsf.expr.subT$group3[i]<-ifelse(ucsf.expr.subT$NEFH[i]>=quantile(ucsf.expr.subT$NEFH,probs = 0.75),"high",ifelse(ucsf.expr.subT$NEFH[i]<=quantile(ucsf.expr.subT$NEFH,probs = 0.25),"low",NA))
}



##join clinical and expr together##
ucsf.join<-full_join(ucsf.expr.subT,ucsf.surv,by=c("sample"="Patient.ID"))
ucsf.join<-na.omit(ucsf.join)
##survival analysis##
fit1<-survfit(Surv(OS.mCRPC,Event)~group3,data=ucsf.join)
res<-ggsurvplot(fit1,data=ucsf.join,
                xlab = "Days",
                ylab = "Overall Survival Probability (%)",
                conf.int=TRUE,
                pval = TRUE,
                fun="pct",
                risk.table = TRUE,
                size=1,
                linetype = "strata"
)
ggsave("ucsf.NEFH.top25bot25.OS.pdf", plot = print(res), width = 8, height = 8, dpi = 500)




##TALE ##
input1<-"MSKCC_PCa_Clinical_Annotation.txt"
input2<-"MSKCC_PCa_mRNA_data.txt"
tale.surv<-read.delim(input1,header = TRUE,stringsAsFactors = FALSE,check.names = FALSE)
tumor<-c("PRIMARY","MET")
tale.surv.sub<-tale.surv[tale.surv$Type%in%tumor,]

tale.expr<-read.delim(input2,header = TRUE,stringsAsFactors = FALSE,check.names = FALSE)
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




##promote##
input1<-"tier2.info.txt"
input2<-"rpkm.count.txt"
promote.surv<-read.delim(input1,header = TRUE,stringsAsFactors = FALSE,check.names = FALSE)
promote.expr<-read.delim(input2,header = TRUE,stringsAsFactors = FALSE,check.names = FALSE)
promote.expr<-promote.expr[promote.expr$GeneID=="NEFH",2:ncol(promote.expr)]
promote.exprT<-as.data.frame(t(promote.expr))
promote.exprT$sample<-rownames(promote.exprT)
promote.exprT<-promote.exprT%>%rename("NEFH"="53407")

for (i in 1:nrow(promote.exprT)){
  promote.exprT$group[i]<-ifelse(promote.exprT$NEFH[i]>=median(promote.exprT$NEFH),"high","low")
}

##join#3
promote.join<-full_join(promote.exprT,promote.surv,by=c("sample"="sample_id"))
promote.join<-promote.join[,c(3,13,14)]

#surv##
fit1<-survfit(Surv(weeks_to_death,death)~group,data=promote.join)
res<-ggsurvplot(fit1,data=promote.join,
                xlab = "Days",
                ylab = "Overall Survival Probability (%)",
                conf.int=TRUE,
                pval = TRUE,
                fun="pct",
                risk.table = TRUE,
                size=1,
                linetype = "strata"
)
ggsave("Promote.NEFH.median.OS.pdf", plot = print(res), width = 8, height = 8, dpi = 500)


##GSE46691##
input<-"GSE46691_family.soft.txt"
dat<-read.delim(input,header=TRUE,stringsAsFactors = FALSE)

#NEFH<-str_which(dat$X.11,"NEFH")
#dat<-dat[str_which(dat$X.11,"NEFH"),]
datsub<-dat[dat$seqname=="chr22",]


##GSE62116##
input<-"NEFH_GPL5188-122.txt"
dat<-read.delim(input,header=FALSE,stringsAsFactors = FALSE)
ID_REF.GSE62116<-as.data.frame(dat$V1)

write.table(ID_REF.GSE62116,"NEFH.ID_REF.GSE62116.txt",sep="\t",col.names = FALSE,row.names = FALSE,quote = FALSE)



