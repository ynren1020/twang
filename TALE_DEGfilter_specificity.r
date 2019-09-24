#############################2019-09-23###############################################
##After TALE DEG analysis, calculate tissue specificity using TCGA pan-cancer data####
######################################################################################

##function for specificity##
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

##data preparation##
args <- commandArgs(TRUE)
#dat<-read.delim("tcga_target_no_normal_rsem_gene_tpm_taleDEG_withsample_first2.txt",header = TRUE,stringsAsFactors = FALSE,check.names = FALSE)
dat<-read.delim("tcga_target_no_normal_rsem_gene_tpm_taleDEG_withsample.txt",header = TRUE,stringsAsFactors = FALSE,check.names = FALSE)
meta<-read.delim("TCGA.metainfo.txt",header = TRUE,stringsAsFactors = FALSE,check.names = FALSE)

##rowname##
rownames(dat)<-dat$sample
dat$sample<-NULL

##transform back to tpm##
fexp<-function(x){
  y<-2^x-0.001
  return(y)
}

dat<-apply(dat,2,fexp)
##transpose##
datT<-as.data.frame(t(dat))
datT$sample<-rownames(datT)

##meta##
meta$sample<-str_sub(meta$ALIQUOT_BARCODE,1,15)
meta.sub<-meta[,4:6]

##join meta and datT##
datT.join<-full_join(datT,meta.sub,by="sample")
datT.join<-na.omit(datT.join)
datT.join$gene<-names(datT.join)[1]

##summary##
# Refer to column names stored as strings with the `.data` pronoun:
#var <- "ENSG00000146083.11"
var <- args[1]
#summarise(starwars, avg = mean(.data[[var]], na.rm = TRUE))
datT.join.sum<-datT.join%>%group_by(PROJECT)%>%
  summarise(count=n(),
            fpkmSum=sum(.data[[var]]))

##mean by dividing all samplesize;specificity##
datT.join.sum$fpkmMean<-datT.join.sum$fpkmSum/datT.join.sum$count
datT.join.sum$tau<-fTau(datT.join.sum$fpkmMean) 
datT.join.sum$zscore<-fZ(datT.join.sum$fpkmMean) 
datT.join.sum$gene<-var

write.table(datT.join.sum,paste0("TCGA_TALE_DEGgenes_specificity_",var),col.names = FALSE,row.names = FALSE,sep = "\t")



