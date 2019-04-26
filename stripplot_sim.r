####################2019-04-25#######################################################
##create https://www.cell.com/cancer-cell/fulltext/S1535-6108(18)30306-4 Figure 4C###
##choose 10 cohort from 33 cohort with 33 defined color in len_boxplot.py############
#####################################################################################

library(dplyr)
library(tidyr)
library(ggplot2)
library(gridExtra)

exitron<-read.delim("TCGA.tumor-specific.exitron.txt",header=TRUE,stringsAsFactors = FALSE)
cohort<-c('ACC','BLCA','BRCA', 'CESC', 'CHOL','COAD','DLBC','ESCA', 'GBM', 'HNSC', 'KICH', 'KIRC', 'KIRP', 'LAML', 'LGG', 'LIHC',  'LUAD', 'LUSC',
          'MESO', 'OV', 'PAAD','PCPG', 'PRAD','READ','SARC','SKCM','STAD', 'TGCT', 'THCA', 'THYM','UCEC', 'UCS', 'UVM')
color<-c('#7b6a4a', '#8049d8',  '#5dd74f', '#cf45cc','#dde73f','#646acd', '#a0d64b','#d3458d', '#67d88e','#dd4529', '#7fdfcd', '#d34058', '#568734',
         '#cd83d9','#bebe49', '#964891', '#dfb23f', '#729bd9','#db842d','#51acc0','#b25c35','#a7cbda', '#8f7b30', '#516590','#c5d994', '#9d5356',
         '#569973', '#dd8999', '#4b7267','#daa476','#d0aed7', '#ccbeaa','#886b82')
cohortcolor<-data_frame(cohort=cohort,color=color)
exitronplot<-left_join(exitron,cohortcolor,by="cohort")
exitronplot<-na.omit(exitronplot)



##create rank of number and generate median of it within each cohort##
exitronplot<-exitronplot %>%
  group_by(cohort) %>%
  mutate(my_ranks = order(order(number, decreasing=FALSE)))%>%
  mutate(my_median = median(number))%>%
  mutate(my_mean_ranks=mean(my_ranks))
##order by median##
exitronplot<-as.data.frame(exitronplot)
exitronplot<-exitronplot[order(exitronplot$my_median),]
##create factor variable of cohort, rename color##
exitronplot$cohort_f<-factor(exitronplot$cohort,levels=unique(exitronplot$cohort))
exitronplot<-rename(exitronplot,my_color=color)
##remove two ends ranks##
exitronplot<-exitronplot%>%group_by(cohort_f)%>%
  filter(my_ranks>5&my_ranks<max(my_ranks)-20)
##select 10 chort for sim plot##
cohortselect<-c("UCEC","BRCA","LUAD","LIHC","BLCA","PRAD","THCA","LGG","SKCM","LUSC")
cohortcolor<-filter(cohortcolor,cohort%in%cohortselect)
##choose LIHC data as example to replicate it 10 times of rank, but add 2,4,...18to next 9 replicated samples##
exitronplot.LIHC<-filter(exitronplot,cohort=="LIHC")


##function to make sim data frame##
simcohort<-function(add,study,color){
  exitronplot.temp<-exitronplot.LIHC%>%mutate(number=number+add,cohort=study,my_color=color)
  exitronplot.temp$cohort_f<-NULL
  exitronplot.temp$cohort_f<-study
  return(exitronplot.temp)
}

exitronplot.UCEC<-simcohort(2,"UCEC","#d0aed7")
exitronplot.BLCA<-simcohort(4,"BLCA","#8049d8")
exitronplot.BRCA<-simcohort(6,"BRCA","#5dd74f")
exitronplot.LGG<-simcohort(8,"LGG","#bebe49")
exitronplot.LUAD<-simcohort(10,"LUAD","#dfb23f")
exitronplot.PRAD<-simcohort(12,"PRAD","#8f7b30")
exitronplot.THCA<-simcohort(14,"THCA","#4b7267")
exitronplot.SKCM<-simcohort(16,"SKCM","#9d5356")
exitronplot.LUSC<-simcohort(18,"LUSC","#729bd9")

##rbind together##
exitronplot.sim<-rbind(exitronplot.LIHC,exitronplot.UCEC,exitronplot.BLCA,exitronplot.BRCA,exitronplot.LGG,exitronplot.LUAD,exitronplot.PRAD,exitronplot.THCA,exitronplot.SKCM,exitronplot.LUSC)
exitronplot.sim$cohort_f<-factor(exitronplot.sim$cohort_f,levels=c("LIHC","UCEC","BLCA","BRCA","LGG","LUAD","PRAD","THCA","SKCM","LUSC"))


##plot##
p<-ggplot(exitronplot.sim,aes(x=my_ranks, y=number,color=cohort_f))+
  geom_point(size=0.5)+
  scale_colour_manual(values=unique(exitronplot.sim$my_color))+
  labs(y="",x="")


p1<-p+facet_grid(.~cohort_f,switch="x",scales = "free_x")+
  #geom_hline(data=data2,aes(yintercept=my_median),linetype="dotted", color = "red")+
  #geom_segment(aes(x=median(my_ranks)-100,xend=median(my_ranks)+50,y=my_median,yend=my_median),colour="red",size=1)+
  theme_bw()+
  theme(axis.title.x=element_blank(),axis.ticks.x=element_blank(),axis.ticks.y=element_blank(),axis.text.y=element_blank(),axis.text.x=element_blank(),legend.position = "none",
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.spacing = unit(0, "mm"),panel.border = element_blank(),                      
        strip.background = element_blank(),axis.line.x = element_line(),axis.line.y = element_line(),strip.text.x = element_blank(),strip.placement = "outside")

p1
#p2<-p1+scale_x_continuous(breaks=exitronplot$my_mean_ranks, labels=exitronplot$cohort)+  
#  theme(axis.text.x=element_text(angle=90, vjust=.5))
#p2



#######################test###################
x<-list(add=seq(2,18,2),study=c("UCEC","BLCA","BRCA","LGG","LUAD","PRAD","THCA","SKCM","LUSC"),color=c("#d0aed7","#8049d8","#5dd74f","#bebe49","#dfb23f","#8f7b30","#4b7267","#9d5356","#729bd9"))

df<-data.frame(add=seq(2,18,2),study=c("UCEC","BLCA","BRCA","LGG","LUAD","PRAD","THCA","SKCM","LUSC"),color=c("#d0aed7","#8049d8","#5dd74f","#bebe49","#dfb23f","#8f7b30","#4b7267","#9d5356","#729bd9"))
simtest<-function(df){
  exitronplot.temp<-exitronplot.LIHC%>%mutate(number=number+df[,1],cohort=df[,2],my_color=df[,3])
  exitronplot.temp$cohort_f<-NULL
  exitronplot.temp$cohort_f<-df[,2]
  return(exitronplot.temp)
}

