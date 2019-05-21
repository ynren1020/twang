########################2019-05-21###########################
##barplot for NEFH expression in prostate cancer cell line###
#############################################################
library(tidyverse)
library(ggpubr)

input<-"mRNA expression (RNAseq)_ NEFH.txt"
dat<-read.delim(input,header=TRUE,stringsAsFactors = FALSE,check.names = FALSE)
names(dat)
prostate<-str_which(names(dat),"PROSTATE")
datsub<-dat[,c(1,prostate)]
datsub<-datsub[,-c(1,10:11)]
datsubT<-as.data.frame(t(datsub))
datsubT<-datsubT%>%rename("NEFH"="V1")
datsubT$cell<-rownames(datsubT)
datsubT$NEFH.round<-round(datsubT$NEFH,2)

##barplot##
p<-ggbarplot(datsubT, x = "cell", y = "NEFH.round",
          fill = "steelblue",
          color = "steelblue",
          xlab = FALSE,
          ylab = "NEFH expression",
          label = TRUE, 
          lab.pos = "out", 
          lab.col = "black",
          #color = "white",            # Set bar border colors to white
          #palette = "jco",            # jco journal color palett. see ?ggpar
          sort.val = "desc",          # Sort the value in dscending order
          sort.by.groups = FALSE,     # Don't sort inside each group
          x.text.angle = 90           # Rotate vertically x axis texts
)
