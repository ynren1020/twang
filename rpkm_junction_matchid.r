########################2019-08-16######################
##choose samples which have characteristic junctions####
################TCGA  ##################################
library(tidyverse)
library(ggpubr)

input<-"TCGA.metainfo.txt"
meta<-read.delim(input,header = TRUE,stringsAsFactors = FALSE,sep = "\t")

input1<-"tcga.normal_tumor.rpkm.id.txt"
rpkm<-read.delim(input1,header = FALSE,stringsAsFactors = FALSE,sep = "\t")

input2<-"tcga.junction.txt"
junction<-read.delim(input2,header = TRUE,stringsAsFactors = FALSE,sep = "\t")

##BRCA . to -##
rpkm.BRCA<-rpkm[rpkm$V3=="BRCA",]
rpkm.BRCA$V1<-str_replace_all(rpkm.BRCA$V1,"[.]","-")

##combine##
rpkm<-rbind(rpkm[rpkm$V3!="BRCA",],rpkm.BRCA)

##meta junction join##
junction.meta<-left_join(junction,meta,by=c("sample"="FILE_ID"))

##join rpkm##
junction.meta.rpkm<-left_join(junction.meta,rpkm,by=c("ALIQUOT_BARCODE"="V1"))
junction.meta.rpkm<-na.omit(junction.meta.rpkm)
write.table(junction.meta.rpkm,"TCGA.rpkm.junction.txt",quote=FALSE,col.names = TRUE,row.names = FALSE,sep = "\t")

df<-junction.meta.rpkm[order(junction.meta.rpkm$V3),]
df$V3<-factor(df$V3,levels = c("BLCA","BRCA","CESC", "CHOL", "COAD", "ESCA", "GBM","HNSC", "KICH", "KIRP", "LIHC", "LUAD", "LUSC", "PAAD", "PRAD", "READ", "SARC", "SKCM", "STAD", "UCEC"))

##ggdotchart##
##,'#51acc0','#b25c35','#a7cbda', '#8f7b30', '#516590','#c5d994', '#9d5356',
'#569973', ', '#ccbeaa','#886b82' ,'#dd8999', '#4b7267','#daa476','#d0aed7'
p <- ggdotchart(df, x = "V3", y = "V2", order = c("BLCA","BRCA","CESC", "CHOL", "COAD", "ESCA", "GBM","HNSC", "KICH", "KIRP", "LIHC", "LUAD", "LUSC", "PAAD", "PRAD", "READ", "SARC", "SKCM", "STAD", "UCEC"),
                color = "V3", palette =c('#7b6a4a', '#8049d8',  '#5dd74f', '#cf45cc','#dde73f','#646acd', '#a0d64b','#d3458d', '#67d88e','#dd4529', '#7fdfcd', '#d34058', '#568734',
                                         '#cd83d9','#bebe49', '#964891', '#dfb23f', '#729bd9','#db842d','#51acc0'),
                ylab = "IGF1R-AS1 RPKM",
                facet.by = "V4"
                
)
p1<-p+theme(legend.position = "none",axis.title.x=element_blank())

ggsave("TCGA.rpkm.tumor.junctionmatched.pdf",width=8,height=6,dpi = 300)

##boxplot##
df<-df%>%group_by(V4,V3)%>%
  count(ncount=n())
df$label<-paste0(df$V3,", ",df$V4, "\n(n=",df$ncount,")")

df1.join<-full_join(junction.meta.rpkm,df,by=c("V3"="V3","V4"="V4"))
df1.join<-df1.join[c(1:306,308,307,309:nrow(df1.join)),]
df1.join<-df1.join[df1.join$V3!="KICH"&df1.join$V3!="KIRP",]
df1.join<-rename(df1.join,"Type"="V4")


##color by cohort##
p2<-ggboxplot(df1.join, "label", "V2", color = "V3",
          palette =c('#7b6a4a', '#8049d8',  '#5dd74f', '#cf45cc','#dde73f','#646acd', '#a0d64b','#d3458d', '#67d88e','#dd4529', '#7fdfcd', '#d34058', '#568734',
                               '#cd83d9','#bebe49', '#964891', '#dfb23f', '#729bd9'),
          
          linetype = 1,  #solid line is 1
          size = 1.5,  #change point and box outline size
          add = c("jitter"),
          alpha = 0.2,
          #add.params = list(color = "red"),
          #xlab = "Dose (mg)",
          ylab = "IGF1R-AS1 expression (FPKM)"
        
)+font("ylab",size=18)+font("y.text",size=18)

#axis.ticks.x=element_blank(),axis.text.x=element_blank()

p2+theme(legend.position = "none",axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank())

ggsave("TCGA.FPKM.boxplot.colorcohort.withoux.pdf",width=10,height=6,dpi = 300)


##color by tumor or normal,redo junction of SKCM and PRAD##
#'#7b6a4a', '#8049d8' '#5dd74f', '#cf45cc','#dde73f','#646acd', '#a0d64b','#d3458d', '#67d88e','#dd4529','#7fdfcd' , '#d34058', '#568734',
#'#cd83d9','#bebe49', '#964891', '#dfb23f', '#729bd9' legend.position = "none",
df1.join$label2<-paste0(df1.join$V3," (n=",df$ncount,")")
p2<-ggboxplot(df1.join, "label", "V2", color = "Type",
              palette =c('#646acd','#dd4529'),
              linetype = 1,  #solid line is 1
              size = 1.5,  #change point and box outline size
              #          notch = TRUE, #default is FALSE,Notches are used to compare groups; if the notches of two boxes do not overlap, this suggests that the medians are significantly different.
              #order = c("2", "1", "0.5"), #change order of items
              add = c("jitter"),
              alpha = 0.2,
              #add.params = list(color = "red"),
              #xlab = "Dose (mg)",
              ylab = "IGF1R-AS1 expression (FPKM)"
              
)

p2+theme(axis.title.x=element_blank(),axis.text.x = element_text(angle = 90),axis.ticks.x=element_blank())
 

ggsave("TCGA.FPKM.tumor.junctionredo.boxplot.color.pdf",width=10,height=6,dpi = 300)


##normal,tumor,meta prostate cancer boxplot and test##
meta.fpkm<-read.delim("meta.fpkm.txt",header = FALSE,stringsAsFactors = FALSE)
df1.join.sub<-df1.join[df1.join$V3=="PRAD",c("V2","V3","Type")]
meta.fpkm$V2<-meta.fpkm$V1
meta.fpkm$V3<-"PRAD"
meta.fpkm$Type<-"Mets"
meta.fpkm$V1<-NULL

prad<-bind_rows(df1.join.sub,meta.fpkm)
prad$Type<-factor(prad$Type,levels=c("normal","tumor","Mets"))
##boxplot##
p <- ggboxplot(prad, x = "Type", y = "V2",
               color = "Type", palette =c("#00AFBB", "#E7B800", "#FC4E07"),
               add = "jitter",
               xlab=FALSE,
               ylab="IGF1R-AS1 expression (FPKM)")
p
# Add p-values comparing groups
# Specify the comparisons you want
#my_comparisons <- list( c("AR_V7 negtive", "AR_V7 positive (high)"), c("AR_V7 positive (high)", "AR_V7 positive (low)"), c("AR_V7 negtive", "AR_V7 positive (low)") )
p1<-p + stat_compare_means(
  method.args = list(alternative = "greater"),
  ref.group = "normal"
) # Add pairwise comparisons p-value

##change tick mark labels##
p2<-p1+scale_x_discrete(labels=c("normal" = "Benign\n(n=2)", "tumor" = "PCa\n(n=90)","Mets" = "Mets\n(n=71)"))


###########Gleason score PCa boxplot#########
score<-read.delim("data_bcr_clinical_data_patient.txt",header = TRUE,stringsAsFactors = FALSE,check.names = FALSE)
score.join<-full_join(score,df1.join,by="PATIENT_ID")
score.join<-na.omit(score.join)
score.join.sub<-score.join[,c("GLEASON_SCORE","V2","Type")]
score.join.sub<-score.join.sub[score.join.sub$Type!="normal",]
score.join.sub$GLEASON_SCORE<-factor(score.join.sub$GLEASON_SCORE,levels=c("6","7","8","9","10"))

table(score.join.sub$GLEASON_SCORE)
#6  7  8  9 10 
#6 41 16 26  1 

##boxplot##
p <- ggboxplot(score.join.sub, x = "GLEASON_SCORE", y = "V2",
               color = "GLEASON_SCORE", palette = get_palette("Dark2", 5),
               add = "jitter",
               xlab=FALSE,
               ylab="IGF1R-AS1 expression (FPKM)")
p
# Add p-values comparing groups
my_comparisons <- list( c("6", "7"), c("6", "8"), c("6", "9"),c("6","10"),c("8","9"),c("7","9") )
p1<-p+stat_compare_means(comparisons = my_comparisons)

# Specify the comparisons you want
#my_comparisons <- list( c("AR_V7 negtive", "AR_V7 positive (high)"), c("AR_V7 positive (high)", "AR_V7 positive (low)"), c("AR_V7 negtive", "AR_V7 positive (low)") )
#p1<-p + stat_compare_means(
#  method.args = list(alternative = "greater"),
#  ref.group = "normal"
#) # Add pairwise comparisons p-value

##change tick mark labels##
#p2<-p1+scale_x_discrete(labels=c("normal" = "Benign\n(n=2)", "tumor" = "PCa\n(n=90)","Mets" = "Mets\n(n=71)"))



