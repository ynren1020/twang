##########################2019-09-18###############################
##TCGA prad normal tumor and su2c meta normalized voom count ######
##"prad.normal.tumor.meta.counts.voom_4boxplot.txt" ###############
###################################################################

library(tidyverse)
library(ggpubr)

dat<-read.delim("prad.normal.tumor.meta.counts.voom_4boxplot.txt",header = TRUE,stringsAsFactors = FALSE,check.names = FALSE)

##choose TMPRSS2##
datsub<-dat[rownames(dat)=="TMPRSS2",]

##transpose##
datsubT<-as.data.frame(t(datsub))
##group##
datsubT$group<-c(rep("Normal",52),rep("Localized",500),rep("Mets",101))
datsubT$group<-factor(datsubT$group,levels=c("Normal","Localized","Mets"))

##plot##
p <- ggboxplot(datsubT, x = "group", y = "TMPRSS2",xlab = FALSE,ylab="TMPRSS2 expression",
               color = "group", palette =c("#00AFBB", "#E7B800", "#FC4E07"),
               add = "jitter")

p
# Change the plot orientation: horizontal
#ggpar(p, orientation = "horiz",xlab="Metastatic/n (N=59)")
# Add p-values comparing groups
# Specify the comparisons you want
my_comparisons <- list( c("Normal", "Localized"), c("Localized", "Mets"), c("Mets", "Normal") )
p1<-p + stat_compare_means(comparisons = my_comparisons)+ # Add pairwise comparisons p-value
  stat_compare_means(label.y = 18)+                   # Add global p-value
  scale_x_discrete(labels=c("Normal" = "Normal\n(N=52)", "Localized" = "Localized\n(N=500)","Mets" = "Mets\n(N=101)"))+
  rotate_x_text(90)+theme(legend.position = "none")
p1

ggsave("TMPRSS2_TCGASU2C.boxplot.pdf",width = 8,height = 8)
