##################2019-07-29#######################
##boxplot for exitron paper########################
###################################################

library(tidyverse)
library(ggpubr)

colony<-c(27,23,19,4,4,8,18,13,21,10,9,5)
lncap<-c("Empty Vector","Empty Vector","Empty Vector","NEFH cDNA","NEFH cDNA","NEFH cDNA","Empty Vector","Empty Vector","Empty Vector","NEFH cDNA","NEFH cDNA","NEFH cDNA")

df<-data.frame(LNCaP=lncap,colony=colony)



##plot logrpkm##
p <- ggbarplot(df, x = "LNCaP", y = "colony",
               color = "LNCaP",palette =c("#00AFBB", "#FC4E07"),
               #fill = "LNCaP",
               add = c("mean_se","jitter"),
               width = 0.3,
               xlab="LNCaP Cells",
               ylab="Colony Number"
               )

p1<-p + stat_compare_means(method="t.test",
                           ref.group = "Empty Vector", 
                           label = "p.format",
                           method.args = list(alternative = "less"),
                           label.y = 20,
                           label.x.npc = "center",
                           size = 8)

# Add p-values comparing groups
# Specify the comparisons you want
#my_comparisons <- list( c("AR_V7 negtive", "AR_V7 positive (high)"), c("AR_V7 positive (high)", "AR_V7 positive (low)"), c("AR_V7 negtive", "AR_V7 positive (low)") )
#p1<-p + stat_compare_means(label.x=1.5,size=8,label = "p.format",tip.length = 0.3) # Add pairwise comparisons p-value

##change tick mark labels##
p2<-p1+theme(axis.text.x = element_text(color = "grey20", size = 17, angle = 0, hjust = .5, vjust = .5, face = "bold"),
             axis.text.y = element_text(color = "grey20", size = 17, angle = 0, hjust = 1, vjust = 0, face = "bold"),  
             axis.title.x = element_text(color = "grey20", size = 22, angle = 0, hjust = .5, vjust = 0, face = "bold"),
             axis.title.y = element_text(color = "grey20", size = 22, angle = 90, hjust = .5, vjust = .5, face = "bold"))

ggsave("LNCaP.NEFH.colony.barplot.dot.pdf",width=8,height=8,dpi = 300)                             
##Ct##

NEFH.ctrl<-c(31.898,32.140,31.989)
NEFH.exp<-c(28.113,28.120,28.207)
B2M.ctrl<-c(15.046,15.085,15.002)
B2M.exp<-c(15.032,14.850,15.307)

del.exp<-NEFH.exp-B2M.exp
del.ctrl<-NEFH.ctrl-B2M.ctrl

del.del.expvsctrl<-del.exp-del.ctrl

ctrl.fc<-2^abs((del.ctrl-mean(del.ctrl)))
exp.fc<-2^abs((del.exp-mean(del.ctrl)))

lncap<-rep(c("Empty Vector","NEFH cDNA"),each =3)
fc<-c(ctrl.fc,exp.fc)
df<-data.frame(lncap=lncap,fc=fc)

p <- ggbarplot(df, x = "lncap", y = "fc",
               color = "lncap", palette =c("#00AFBB", "#FC4E07"),
               fill = "lncap",
               add = c("mean_se"),
               width = 0.3,
               xlab="LNCaP cells",
               ylab="Relative NEFH expression\n over control")
p1<-p + stat_compare_means(method="t.test",
                           ref.group = "Empty Vector", 
                           label = "p.format",
                           method.args = list(alternative = "great"),
                           label.y = 20,
                           label.x.npc = "center",
                           size = 8)

p2<-p1+theme(axis.text.x = element_text(color = "grey20", size = 17, angle = 0, hjust = .5, vjust = .5, face = "bold"),
             axis.text.y = element_text(color = "grey20", size = 17, angle = 0, hjust = 1, vjust = 0, face = "bold"),  
             axis.title.x = element_text(color = "grey20", size = 22, angle = 0, hjust = .5, vjust = 0, face = "bold"),
             axis.title.y = element_text(color = "grey20", size = 22, angle = 90, hjust = .5, vjust = .5, face = "bold"))

ggsave("LNCaP.NEFH.expression.bar.pdf",width=8,height=8,dpi = 300) 




# Add p-values comparing groups
# Specify the comparisons you want
#my_comparisons <- list( c("AR_V7 negtive", "AR_V7 positive (high)"), c("AR_V7 positive (high)", "AR_V7 positive (low)"), c("AR_V7 negtive", "AR_V7 positive (low)") )
p1<-p + stat_compare_means(
                           method.args = list(alternative = "great"),
                           ref.group = "Normal"
) # Add pairwise comparisons p-value

##change tick mark labels##
p2<-p1+scale_x_discrete(labels=c("Tumor" = "Tumor\n(N=422)", "Normal" = "Normal\n(N=114)"))

ggpar(p2,legend.title = "category")



