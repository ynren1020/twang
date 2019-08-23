####################2019-08-23#######################
##calcualte igf1r-as1 tissue specificity#############
##tau from "A benchmark of gene expression tissue-###
##specificity metrics"                  #############
#####################################################
BiocManager::install("preprocessCore")
library("preprocessCore") #did not use quantile normalization#
library(tidyverse)
###+++###
#Function requires data frame to be normalized
#1. All 0 are set to NA, to exclude them from quatile normalization
#2. Data are quantile normalized
#3. 0 values (the one set to NA) are set back to 0
fQN <- function(x) #
{
  x[x==0] <- NA
  x_m <- as.matrix(x)
  x <- normalize.quantiles(x_m)
  x[is.na(x)] <- 0
  return(data.frame(x))
}	
###***###***###

#Function require a vector with expression of one gene in different tissues.
#If expression for one tissue is not known, gene specificity for this gene is NA
#Minimum 2 tissues
fTau <- function(x)
{
  if(all(!is.na(x)))
  {
    if(min(x, na.rm=TRUE) >= 0)
    {
      if(max(x)!=0)
      {
        x <- (1-(x/max(x)))
        res <- sum(x, na.rm=TRUE)
        res <- res/(length(x)-1)
      } else {
        res <- 0
      }
    } else {
      res <- NA
      #print("Expression values have to be positive!")
    } 
  } else {
    res <- NA
    #print("No data for this gene avalable.")
  } 
  return(res)
}
###***###***###

###+++###
#Z-score
#Function require a vector with expression of one gene in different tissues.
#If expression for one tissue is not known, gene specificity for this gene is NA
fZ <- function(x)
{	
  res<-(x-mean(x))/sd(x)
  return(res)
}
###***###***###	


##tcga cohort tumor patient's igf1r-as1##
tumor<-read.delim("tcga.tumor.rpkm.id.txt",header = FALSE,stringsAsFactors = FALSE)
tumor<-na.omit(tumor)

df<-tumor%>%group_by(V3)%>%
  summarize(fpkmSum=sum(V2))

##no igf1ras1 expression cohort as 0##
for (i in 1:length(unique(meta$PROJECT))){
cohorts[i]<-strsplit(unique(meta$PROJECT)[i],"-")[[1]][2]
}

df2<-data.frame(V3=cohorts,V2=NA)

df.join<-full_join(df,df2,by="V3")
df.join$V2<-NULL
df.join[is.na(df.join)]<-0

##mean obtained by sum/allsamplesize##
##nofilter, but use all sample size to average##
meta.tumor<-filter(meta,TUMOR_TYPE!="Solid Tissue Normal")
meta.tumor.sum<-meta.tumor%>%group_by(PROJECT)%>%
  summarise(count=n())

for (i in 1:nrow(meta.tumor.sum)){
  meta.tumor.sum$V3[i]<-strsplit(meta.tumor.sum$PROJECT[i],"-")[[1]][2]
}

df.join.meta<-full_join(df.join,meta.tumor.sum,by="V3")
##mean by dividing all samplesize##
df.join.meta$fpkmMean<-df.join.meta$fpkmSum/df.join.meta$count
df.join.meta$tau<-fTau(df.join.meta$fpkmMean) #0.8823
df.join.meta$zscore<-fZ(df.join.meta$fpkmMean) #BRCA 3.56841391 SKCM 2.58196487 PRAD 2.44558988

write.table(df.join,"IGF1RAS1-FPKM-allsample-specificity.txt",quote = FALSE,row.names = FALSE,col.names = TRUE,sep = "\t")


##logtransform##
tumor$logV2<-log2(tumor$V2)
tumor<-tumor[tumor$logV2!="-Inf",]
df<-tumor%>%group_by(V3)%>%
  summarize(fpkmSum=sum(logV2),fpkmMean=mean(logV2))
df$tau<-fTau(-df$fpkmMean)

########logtransform fpkm>1 all samplesize#######
tumor.f<-filter(tumor,V2>=1)
tumor.f$logV2<-log2(tumor.f$V2)
df<-tumor.f%>%group_by(V3)%>%
  summarize(fpkmSum=sum(logV2),fpkmMeansub=mean(logV2))

df.join<-full_join(df,meta.tumor.sum,by="V3")
df.join<-na.omit(df.join)
df.join$fpkmMean<-df.join$fpkmSum/df.join$count
##all sample size counted for Mean, PRAD specific##
df.join$tau<-fTau(df.join$fpkmMean)
df.join$zscore<-fZ(df.join$fpkmMean)
##only fpkm>1 sample size counted for meansub## BREA specific##
df.join$tausub<-fTau(df.join$fpkmMeansub)
df.join$zscoresub<-fZ(df.join$fpkmMeansub)


##tau specificity for raw fpkm##
tau.igf1ras1<-fTau(df.join$fpkmMean)  #0.881
##zscore##
zscore.igf1ras1<-fZ(df.join$fpkmMean)

##output##
df.join$tau<-tau.igf1ras1
df.join$zscore<-zscore.igf1ras1
df.join<-df.join%>%arrange(desc(zscore))

write.table(df.join,"IGF1RAS1-specificity.txt",quote = FALSE,row.names = FALSE,col.names = TRUE,sep = "\t")


##filter out low expression patients by fpkm <1##
tumor.f<-filter(tumor,V2>=1)

df<-tumor.f%>%group_by(V3)%>%
  summarize(fpkmMean=mean(V2))

tau.igf1ras1<-fTau(df.join$fpkmMean) #0.9795021
zscore.igf1ras1<-fZ(df.join$fpkmMean)

df.join$tau<-tau.igf1ras1
df.join$zscore<-zscore.igf1ras1
df.join<-df.join%>%arrange(desc(zscore))

write.table(df.join,"IGF1RAS1-specificity.filter.txt",quote = FALSE,row.names = FALSE,col.names = TRUE,sep = "\t")


##junction matched##
input2<-"tcga.junction.txt"
junction<-read.delim(input2,header = TRUE,stringsAsFactors = FALSE,sep = "\t")
##join##
junction.meta<-full_join(junction,meta.tumor,by=c("sample"="FILE_ID"))
junction.meta<-na.omit(junction.meta)
junction.meta.tumor<-full_join(junction.meta,tumor,by=c("ALIQUOT_BARCODE"="V1"))
junction.meta.tumor<-na.omit(junction.meta.tumor)

##specificity by V2, FPKM##
df<-junction.meta.tumor%>%group_by(V3)%>%
  summarize(fpkmSum=sum(V2))

df2<-data.frame(V3=cohorts,V2=NA)

df.join<-full_join(df,df2,by="V3")
df.join$V2<-NULL
df.join[is.na(df.join)]<-0

df.join.meta<-full_join(df.join,meta.tumor.sum,by="V3")
##mean by dividing all samplesize##
df.join.meta$fpkmMean<-df.join.meta$fpkmSum/df.join.meta$count
df.join.meta$tau<-fTau(df.join.meta$fpkmMean) #0.8823
df.join.meta$zscore<-fZ(df.join.meta$fpkmMean) #SKCM 3.90106452 PRAD 3.64559388
df.join.meta<-df.join.meta%>%arrange(desc(zscore))
##final results (use junction matched FPKM for tissue specific calculation)

write.table(df.join.meta,"IGF1RAS1-junctionmatched-specificity.txt",quote = FALSE,row.names = FALSE,col.names = TRUE,sep = "\t")



