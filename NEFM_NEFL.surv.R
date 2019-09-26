##############################2019-09-26################################
##check NEFM NEFL expression trend and survival pattern in tale ########
########################################################################

library(tidyverse)
library(survival)
library(survminer)


##TALE ##
input1<-"MSKCC_PCa_Clinical_Annotation.txt"
input2<-"MSKCC_PCa_mRNA_data.txt"
tale.surv<-read.delim(input1,header = TRUE,stringsAsFactors = FALSE,check.names = FALSE)
tumor<-c("PRIMARY","MET")
tale.surv.sub<-tale.surv[tale.surv$Type%in%tumor,]

tale.expr<-read.delim(input2,header = TRUE,stringsAsFactors = FALSE,check.names = FALSE)
tale.expr.sub<-tale.expr[tale.expr$GeneSymbol=="NEFM",]
tale.expr.sub<-na.omit(tale.expr.sub)
rownames(tale.expr.sub)<-tale.expr.sub$GeneSymbol
tale.expr.sub$GeneID<-tale.expr.sub$GeneSymbol<-NULL

tale.expr.subT<-as.data.frame(t(tale.expr.sub))
tale.expr.subT$sample<-rownames(tale.expr.subT)


##median group##
for (i in 1:nrow(tale.expr.subT)){
  tale.expr.subT$group[i]<-ifelse(tale.expr.subT$NEFM[i]>=median(tale.expr.subT$NEFM),"high","low")
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

tale.join<-tale.join[!is.na(tale.join$NEFM),]
#tale.join<-tale.join[,c("SurvTime","status","group")]
#tale.join<-na.omit(tale.join)
##only primary tumor##
#tale.join<-tale.join[tale.join$Type=="PRIMARY",]
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
ggsave("tale.NEFM.median.BCR_Free.pdf", plot = print(res), width = 8, height = 8, dpi = 500)

fit2<-survfit(Surv(SurvTime,status)~group,data=tale.join)
res<-ggsurvplot(fit2,data=tale.join,
                xlab = "Days",
                ylab = "Overall Survival Probability (%)",
                conf.int=TRUE,
                pval = TRUE,
                fun="pct",
                risk.table = TRUE,
                size=1,
                linetype = "strata"
)
ggsave("tale.NEFM.median.OS.pdf", plot = print(res), width = 8, height = 8, dpi = 500)



##boxplot##
tale.join<-full_join(tale.expr.subT,tale.surv.sub,by=c("sample"="Sample ID"))
tale.join.box<-tale.join[1:179,1:4]

for (i in 1:nrow(tale.join.box)){
  tale.join.box$group[i]<-ifelse(tale.join.box$Type[i]=="PRIMARY","Localized",ifelse(tale.join.box$Type[i]=="MET","MET","Normal"))
}
tale.join.box$group[is.na(tale.join.box$group)]<-"Normal"
tale.join.box$group<-factor(tale.join.box$group,levels = c("Normal","Localized","MET"))
#localized 131
#MET 19
#Normal 29
p <- ggboxplot(tale.join.box, x = "group", y = "NEFM",
               color = "group", palette =c("#00AFBB", "#E7B800","#FC4E07"),
               add = "jitter",
               xlab=FALSE,
               ylab="NEFM expressoin",
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
  theme(text = element_text(size=20),legend.position = "none")
p2


ggsave("TALE_NEFM.boxplot.pval.updated.pdf",width = 8,height = 8, dpi = 300)