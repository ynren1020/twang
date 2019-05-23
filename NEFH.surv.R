####################2019-05-23####################
##NEFH expression ~ survival######################
##UCSF, promote, and TALE#########################
##UCSF:mmc5,2018_04_15_matrix_rna_tpm.txt#########
##PROMOTE:tier.info.txt, /home/tywang/Projects/lncRNAs/PROMOTE/known_counts/results/V1/rpkm.count.txt
##TALE:MSKCC_****,xls,MSKCC**.txt#################
library(tidyverse)
library(survival)
library(survminer)


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
tale.surv<-tale.surv[1:232,1:39]
tale.join<-full_join(tale.expr.subT,tale.surv,by=c("sample"="Sample ID"))
for (i in 1:nrow(tale.join)){
  tale.join$status[i]<-ifelse(tale.join$Event[i]=="NO",0,1)
}
##BCR##
for (i in 1:nrow(tale.join)){
  tale.join$status2[i]<-ifelse(tale.join$BCR_Event[i]=="NO",0,1)
}

tale.join<-tale.join[,c("SurvTime","status","group")]
tale.join<-na.omit(tale.join)
#surv##
fit1<-survfit(Surv(BCR_FreeTime,status2)~group,data=tale.join)
res<-ggsurvplot(fit1,data=tale.join,
                xlab = "Days",
                ylab = "Overall Survival Probability (%)",
                conf.int=TRUE,
                pval = TRUE,
                fun="pct",
                risk.table = TRUE,
                size=1,
                linetype = "strata"
)
ggsave("tale.NEFH.median.BCR_FreeTime.pdf", plot = print(res), width = 8, height = 8, dpi = 500)






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


