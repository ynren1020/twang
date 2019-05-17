##################2019-05-17#######################
##boxplot for exitron paper########################
###################################################

library(tidyverse)
library(ggpubr)

input1<-"NEFH.exp.txt"
input2<-"NEFH.T.N.txt"

rnaseq<-read.delim(input1,header=TRUE,stringsAsFactors = FALSE) #Tumor (500) #Normal (172)
exitron<-read.delim(input2,header=TRUE,stringsAsFactors = FALSE) #Tumor(422) Normal (114)

rnaseq$logrpkm<-log10(rnaseq$RPKM)
##plot logrpkm##
p <- ggboxplot(rnaseq, x = "category", y = "logrpkm",
               color = "category", palette =c("#00AFBB", "#E7B800"),
               add = "jitter",
               xlab=FALSE,
               ylab="log(RPKM)")

# Add p-values comparing groups
# Specify the comparisons you want
#my_comparisons <- list( c("AR_V7 negtive", "AR_V7 positive (high)"), c("AR_V7 positive (high)", "AR_V7 positive (low)"), c("AR_V7 negtive", "AR_V7 positive (low)") )
p1<-p + stat_compare_means(label.x = 2.0,
  method.args = list(alternative = "less"),
                           ref.group = "Normal"
                           ) # Add pairwise comparisons p-value

##change tick mark labels##
p2<-p1+scale_x_discrete(labels=c("Tumor" = "Tumor\n(N=500)", "Normal" = "Normal\n(N=172)"))
                                
##exitron##
p <- ggboxplot(exitron, x = "cohort", y = "number",
               color = "cohort", palette =c("#00AFBB", "#E7B800"),
               add = "jitter",
               xlab=FALSE,
               ylab="Number of Exitron")
p
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



