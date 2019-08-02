####################2019-08-01#######################################################
##create https://www.cell.com/cancer-cell/fulltext/S1535-6108(18)30306-4 Figure 4C###
##33 cohort with 33 defined color in len_boxplot.py##################################
#####################################################################################
devtools::install_github("thomasp85/patchwork")
library(dplyr)
library(tidyr)
library(ggplot2)
library(gridExtra)
library(patchwork)

exitron<-read.delim("TCGA.tumor-specific.exitron.txt",header=TRUE,stringsAsFactors = FALSE)
exitron.iffs<-read.delim("TCGA.tumor-specific.IF_FS.exitron.txt",header = TRUE,stringsAsFactors = FALSE)

cohort<-c('ACC','BLCA','BRCA', 'CESC', 'CHOL','COAD','DLBC','ESCA', 'GBM', 'HNSC', 'KICH', 'KIRC', 'KIRP', 'LAML', 'LGG', 'LIHC',  'LUAD', 'LUSC',
          'MESO', 'OV', 'PAAD','PCPG', 'PRAD','READ','SARC','SKCM','STAD', 'TGCT', 'THCA', 'THYM','UCEC', 'UCS', 'UVM')
color<-c('#7b6a4a', '#8049d8',  '#5dd74f', '#cf45cc','#dde73f','#646acd', '#a0d64b','#d3458d', '#67d88e','#dd4529', '#7fdfcd', '#d34058', '#568734',
         '#cd83d9','#bebe49', '#964891', '#dfb23f', '#729bd9','#db842d','#51acc0','#b25c35','#a7cbda', '#8f7b30', '#516590','#c5d994', '#9d5356',
         '#569973', '#dd8999', '#4b7267','#daa476','#d0aed7', '#ccbeaa','#886b82')
#color_bc<-c(rep(c(1,2),16),1)
cohortcolor<-data_frame(cohort=cohort,color=color)

exitronplot<-left_join(exitron.iffs,cohortcolor,by="cohort")
exitronplot.1<-exitronplot[,c(1,2,3,5,6)]
exitronplot.2<-exitronplot[,c(1,2,4,5,6)]

##create rank of number and generate median of it within each cohort##
exitronplot.1m<-exitronplot.1 %>%
  group_by(cohort) %>%
  mutate(my_ranks = order(order(number, decreasing=FALSE)))%>%
  mutate(my_median = median(number))%>%
  mutate(my_mean_ranks=mean(my_ranks))%>%
  mutate(my_median_ranks=median(my_ranks))%>%
  mutate(ranks = order(order(inframe, decreasing=FALSE)))%>%
  mutate(condition = "inframe")%>%
  rename(numbers_bygroup=inframe)
  

exitronplot.2m<-exitronplot.2 %>%
  group_by(cohort) %>%
  mutate(my_ranks = order(order(number, decreasing=FALSE)))%>%
  mutate(my_median = median(number))%>%
  mutate(my_mean_ranks=mean(my_ranks))%>%
  mutate(my_median_ranks=median(my_ranks))%>%
  mutate(ranks = order(order(frameshift, decreasing=FALSE)))%>%
  mutate(condition = "frameshift")%>%
  rename(numbers_bygroup=frameshift)

exitronplot<-rbind(exitronplot.1m,exitronplot.2m)

exitronplot<-exitronplot%>%
  group_by(cohort,condition)%>%
  mutate(my_medianbygroup=median(numbers_bygroup))%>%
  mutate(my_median_ranksbygroup=median(ranks))
##order by median##
exitronplot<-as.data.frame(exitronplot)
exitronplot<-exitronplot[order(exitronplot$my_median),]
##create factor variable of cohort, rename color##
exitronplot$cohort_f<-factor(exitronplot$cohort,levels=unique(exitronplot$cohort))
exitronplot<-rename(exitronplot,my_color=color)
exitronplot$condition<-as.factor(exitronplot$condition)


##plot##
p<-ggplot(exitronplot,aes(x=ranks, y=numbers_bygroup,color=cohort_f))+
  geom_point(size=ifelse(exitronplot$condition=="inframe",1.5,0.1),shape=ifelse(exitronplot$condition=="inframe",16,1))+
  scale_colour_manual(values=unique(exitronplot$my_color))+
  labs(y="Number of Exitron")

#pp<-ggplot(exitronplot,aes(x=my_ranks, y=frameshift,color=cohort_f))+
#  geom_point(size=0.5,alpha=0.5)+
#  scale_colour_manual(values=unique(exitronplot$my_color))+
#  labs(y=NULL)

#ppp<-p+pp

p1<-p+facet_grid(.~cohort_f+condition,switch="x",scales = "free_x")+
  #geom_hline(data=data2,aes(yintercept=my_median),linetype="dotted", color = "red")+
  geom_segment(aes(x=my_median_ranksbygroup-my_median_ranksbygroup*0.3,xend=my_median_ranksbygroup+my_median_ranksbygroup*0.3,y=my_medianbygroup,yend=my_medianbygroup),colour="red",size=0.5)+
  theme_bw()+
  theme(axis.title.x=element_blank(),axis.ticks.x=element_blank(),axis.text.x=element_blank(),legend.position = "none",
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.spacing = unit(0, "mm"),panel.border = element_blank(),                      
        strip.background = element_blank(),axis.line.x = element_line(),axis.line.y = element_line(),strip.text.x = element_text(angle = 90),strip.placement = "outside")

p1

ggsave("tcga.inframe.frameshift.scatterplot.pdf",width = 8,height = 8,dpi = 300)


###########color by frameshift or inframe and facet by cohort##

p<-ggplot(exitronplot,aes(x=ranks, y=numbers_bygroup,color=ifelse(condition=="inframe","red","green")))+
 # scale_colour_manual(values=c("green","red"))+
  geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, fill=(exitronplot$gridcolor)))+
  scale_fill_manual( values = c("grey"="grey","white"="white"))+
  geom_point(size=0.5)+
  labs(y="Number of Exitron")



p1<-p+facet_grid(.~cohort_f+condition,switch="x",scales = "free_x")+
  
  #geom_hline(data=data2,aes(yintercept=my_median),linetype="dotted", color = "red")+
  geom_segment(aes(x=my_median_ranksbygroup-my_median_ranksbygroup*0.3,xend=my_median_ranksbygroup+my_median_ranksbygroup*0.3,y=my_medianbygroup,yend=my_medianbygroup),colour="black",size=0.5)+
  theme_bw()+
  theme(axis.title.x=element_blank(),axis.ticks.x=element_blank(),axis.text.x=element_blank(),legend.position = "none",
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.spacing = unit(0, "mm"),panel.border = element_blank(),                      
        strip.background = element_blank(),axis.line.x = element_line(),axis.line.y = element_line(),strip.text.x = element_text(angle = 90),strip.placement = "outside")

p1

ggsave("tcga.inframe.frameshift.scatterplot.twocolor.pdf",width = 8,height = 8,dpi = 300)


##add grid line##
odd_cohortf<-unique(exitronplot$cohort_f)[seq(1,length(unique(exitronplot$cohort_f)),2)]
even_cohortf<-unique(exitronplot$cohort_f)[seq(2,length(unique(exitronplot$cohort_f)),2)]
for (i in 1:nrow(exitronplot)){
exitronplot$gridcolor[i]<-ifelse(exitronplot$cohort[i]%in%odd_cohortf,"grey","white")
}

p1<-p+
  facet_grid(.~cohort_f+condition,switch="x",scales = "free_x")+
  #geom_hline(data=data2,aes(yintercept=my_median),linetype="dotted", color = "red")+
  geom_segment(aes(x=my_median_ranksbygroup-my_median_ranksbygroup*0.3,xend=my_median_ranksbygroup+my_median_ranksbygroup*0.3,y=my_medianbygroup,yend=my_medianbygroup),colour="black",size=0.5)+
  theme_bw(base_size = 10)+
  theme(axis.title.x=element_blank(),axis.ticks.x=element_blank(),axis.text.x=element_blank(),legend.position = "none",
        panel.spacing = unit(0, "mm"),panel.grid.minor = element_blank(), panel.grid.major = element_blank(), panel.border = element_blank(),              
        strip.background = element_blank(),axis.line.x = element_line(),axis.line.y = element_line(),strip.text.x = element_text(angle = 90),strip.placement = "outside") 
  

p1

ggsave("tcga.inframe.frameshift.scatterplot.gridline.pdf",width = 8,height = 8,dpi = 300)

#################Finally, decide to use this one: test if this one can remove facet grid line################
p<-ggplot(exitronplot,aes(x=ranks, y=numbers_bygroup))+
  facet_grid(.~cohort_f+condition,switch="x",scales = "free_x")+
  # scale_colour_manual(values=c("green","red"))+
  geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, fill=(exitronplot$gridcolor)))+
  scale_fill_manual( values = c("grey"="grey","white"="white",alpha=0.001))+
  geom_point(size=0.5,color=ifelse(exitronplot$condition=="inframe","red","green"))+
  labs(y="Number of Exitron")

p1<-p+geom_segment(aes(x=my_median_ranksbygroup-my_median_ranksbygroup*0.3,xend=my_median_ranksbygroup+my_median_ranksbygroup*0.3,y=my_medianbygroup,yend=my_medianbygroup),colour="black",size=0.5)+
  theme_bw(base_size = 12)+
  theme(axis.title.x=element_blank(),axis.ticks.x=element_blank(),axis.text.x=element_blank(),legend.position = "none",
        panel.spacing = unit(0, "mm"), panel.border = element_rect(size = 0.1),                    
        strip.background = element_blank(),axis.line.x = element_line(),axis.line.y = element_line(),strip.text.x = element_text(angle = 90),strip.placement = "outside")

p1

ggsave("tcga.inframe.frameshift.scatterplot.greywhite.pdf",width = 13,height = 8,dpi = 300)





