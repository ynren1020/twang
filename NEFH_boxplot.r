##########################2019-05-24####################
##SU2C,PRAD(normal,primary) NEFH expression boxplot#####
########################################################

library(tidyverse)
library(ggpubr)

input1<-"rpkm_prad.n.join.txt"
input2<-"rpkm_prad.t.join.txt"
input3<-"rpkm_su2c.join.txt"
prad.n<-read.delim(input1,header = TRUE,stringsAsFactors = FALSE,check.names = FALSE)
prad.t<-read.delim(input2,header = TRUE,stringsAsFactors = FALSE,check.names = FALSE)
su2c<-read.delim(input3,header = TRUE,stringsAsFactors = FALSE,check.names = FALSE)

##NEFH
prad.n.sub<-prad.n[prad.n$V2=="NEFH",]
prad.t.sub<-prad.t[prad.t$V2=="NEFH",]
su2c.sub<-su2c[su2c$V2=="NEFH",]

prad.n.sub<-na.omit(prad.n.sub)
prad.t.sub<-na.omit(prad.t.sub)
su2c.sub<-na.omit(su2c.sub)

##transpose##
prad.n.sub<-prad.n.sub[,2:(ncol(prad.n.sub)-6)]
prad.t.sub<-prad.t.sub[,2:(ncol(prad.t.sub)-6)]
su2c.sub<-su2c.sub[,2:(ncol(su2c.sub)-5)]

prad.n.subT<-as.data.frame(t(prad.n.sub))
prad.t.subT<-as.data.frame(t(prad.t.sub))
su2c.subT<-as.data.frame(t(su2c.sub))

prad.n.subT<-prad.n.subT%>%rename("NEFH"="2246")
prad.n.subT$group<-"Normal (N=52)"
prad.t.subT<-prad.t.subT%>%rename("NEFH"="2246")
prad.t.subT$group<-"Localized (N=491)"
su2c.subT<-su2c.subT%>%rename("NEFH"="57173")
su2c.subT$group<-"Metastasis (N=58)"

##rbind##
all.join<-rbind(prad.n.subT,prad.t.subT,su2c.subT)
all.join$group<-factor(all.join$group,levels = c("Normal (N=52)","Localized (N=491)","Metastasis (N=58)"))
all.join$logNEFH<-log10(all.join$NEFH)
##plot##
p <- ggboxplot(all.join, x = "group", y = "logNEFH",xlab = FALSE,ylab="NEFH expression",
               color = "group", palette =c("#00AFBB", "#E7B800", "#FC4E07"),
               add = "jitter", shape = "group")
  
p
# Change the plot orientation: horizontal
#ggpar(p, orientation = "horiz",xlab="Metastatic/n (N=59)")
# Add p-values comparing groups
# Specify the comparisons you want
my_comparisons <- list( c("Normal (N=52)", "Localized (N=491)"), c("Localized (N=491)", "Metastasis (N=58)"), c("Metastasis (N=58)", "Normal (N=52)") )
p1<-p + stat_compare_means(comparisons = my_comparisons)+ # Add pairwise comparisons p-value
  stat_compare_means(label.y = 6)                   # Add global p-value
p1

ggsave("NEFH_TCGASU2C.boxplot.pdf",width = 8,height = 8)




