######################2019-08-14######################
############boxplot for tcga rpkm#####################
######################################################
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)

input<-"tcga.tumor_normal.rpkm.txt"
df<-read.delim(input, header = FALSE,stringsAsFactors = FALSE)
df<-na.omit(df)
df$logV1<-log(df$V1)
df.normal<-df[df$V3=="normal",]
df.tumor<-df[df$V3=="tumor",]
df.tumor.sub<-df.tumor[df.tumor$V2=="BRCA"&df.tumor$V1<5.9,]
df.tumor.sub2<-df.tumor[df.tumor$V2!="BRCA",]
df.tumor.subC<-rbind(df.tumor.sub,df.tumor.sub2)

##remove SKCM outlier##
df.tumor.subC1<-df.tumor.subC[df.tumor.subC$V2=="SKCM"&df.tumor.subC$V1<10,]
df.tumor.subC2<-df.tumor.subC[df.tumor.subC$V2!="SKCM",]
df.tumor.subCC<-rbind(df.tumor.subC1,df.tumor.subC2)

##LAML 5018
df.tumor.subCC1<-df.tumor.subCC[df.tumor.subCC$V2=="LAML"&df.tumor.subCC$V1<5.0,]
df.tumor.subCC2<-df.tumor.subCC[df.tumor.subCC$V2!="LAML",]
df.tumor.subCCC<-rbind(df.tumor.subCC1,df.tumor.subCC2)

p <- ggdotchart(df.normal, x = "V2", y = "logV1",
               color = "V2", palette =c('#7b6a4a', '#8049d8',  '#5dd74f', '#cf45cc','#dde73f','#646acd', '#a0d64b','#d3458d', '#67d88e','#dd4529', '#7fdfcd', '#d34058', '#568734',
                                        '#cd83d9','#bebe49', '#964891', '#dfb23f', '#729bd9','#db842d','#51acc0','#b25c35','#a7cbda', '#8f7b30', '#516590','#c5d994', '#9d5356',
                                        '#569973', '#dd8999', '#4b7267','#daa476','#d0aed7', '#ccbeaa','#886b82'),
               ylab = "IGF1R-AS1 log(RPKM)"
               
               )
p1<-p+theme(legend.position = "none",axis.title.x=element_blank())

ggsave("TCGA.rpkm.normal.logscale.pdf",width=8,height=6,dpi = 300)

##tumor##
p <- ggdotchart(df.tumor.subCCC[order(df.tumor.subCCC$V2),], x = "V2", y = "V1",
                color = "V2", palette =c('#7b6a4a', '#8049d8',  '#5dd74f', '#cf45cc','#dde73f','#646acd', '#a0d64b','#d3458d', '#67d88e','#dd4529', '#7fdfcd', '#d34058', '#568734',
                                         '#cd83d9','#bebe49', '#964891', '#dfb23f', '#729bd9','#db842d','#51acc0','#b25c35','#a7cbda', '#8f7b30', '#516590','#c5d994', '#9d5356',
                                         '#569973', '#dd8999', '#4b7267','#daa476','#d0aed7', '#ccbeaa','#886b82'),
                ylab = "IGF1R-AS1 RPKM"
                
                
)
p1<-p+theme(legend.position = "none",axis.title.x=element_blank())

ggsave("TCGA.rpkm.tumor.outlierreomove.pdf",width=8,height=6,dpi = 300)



#############################BRCA tumor data recheck#######################################
input2<-"BRCA_tumor.rpkm.txt"
brca<-read.delim(input2,header = TRUE,stringsAsFactors = FALSE)
brca.sub<-brca[brca$GeneID=="PB.69",]
rownames(brca.sub)<-brca.sub$GeneID
brca.sub$GeneID<-NULL
brca.sub$Length<-NULL
brca.subT<-as.data.frame(t(as.matrix(brca.sub)))
brca.subT<-na.omit(brca.subT)

tops<-top_n(brca.subT,5,PB.69)
PB.69
1  18.271675 TCGA-AC-A3QQ-01B-06R-A22O-07
2   3.815868
3  39.281845 TCGA-A7-A13G-01B-04R-A22O-07

4 162.511133
5   5.906650 TCGA-AC-A2QH-01B-04R-A22O-07

input3<-"BRCA.tumor.txt"
brca.r<-read.delim(input3,header = FALSE,stringsAsFactors = FALSE)

##remove 5.9 above patients##

##recheck PRAD tumor##
input2<-"PRAD_tumor.rpkm.txt"
PRAD<-read.delim(input2,header = TRUE,stringsAsFactors = FALSE)
PRAD.sub<-PRAD[PRAD$GeneID=="PB.69",]
rownames(PRAD.sub)<-PRAD.sub$GeneID
PRAD.sub$GeneID<-NULL
PRAD.sub$Length<-NULL
PRAD.subT<-as.data.frame(t(as.matrix(PRAD.sub)))
PRAD.subT<-na.omit(PRAD.subT)

tops<-top_n(PRAD.subT,5,PB.69)
PB.69
1 4.997459
2 3.701073
3 3.717243
4 3.719253
5 9.183603 TCGA-WW-A8ZI-01A-11R-A37L-07 (reliable, check bam file in igv, can detect junctions)

##recheck SKCM tumor ##
input2<-"SKCM_tumor.rpkm.txt"
SKCM<-read.delim(input2,header = TRUE,stringsAsFactors = FALSE)
SKCM.sub<-SKCM[SKCM$GeneID=="PB.69",]
rownames(SKCM.sub)<-SKCM.sub$GeneID
SKCM.sub$GeneID<-NULL
SKCM.sub$Length<-NULL
SKCM.subT<-as.data.frame(t(as.matrix(SKCM.sub)))
SKCM.subT<-na.omit(SKCM.subT)

tops<-top_n(SKCM.subT,5,PB.69)
PB.69
1  1.856620
2  2.292124
3  2.081816
4  1.876222
5 10.342661 TCGA-ER-A3PL-06A-11R-A239-07 (this one does not have the specific junction, should be removed)

##recheck LAML tumor ##
input2<-"LAML_tumor.rpkm.txt"
LAML<-read.delim(input2,header = TRUE,stringsAsFactors = FALSE)
LAML.sub<-LAML[LAML$GeneID=="PB.69",]
rownames(LAML.sub)<-LAML.sub$GeneID
LAML.sub$GeneID<-NULL
LAML.sub$Length<-NULL
LAML.subT<-as.data.frame(t(as.matrix(LAML.sub)))
LAML.subT<-na.omit(LAML.subT)

tops<-top_n(LAML.subT,5,PB.69)
PB.69
1 2.635927
2 2.841299
3 3.360241
4 2.658916
5 5.123639 TCGA-AB-2817-03A-01T-0736-13 896b1ee3-adaf-42d4-ab3a-282714e821ad.bam, remove this one
  