###########################2019-08-08#########################
##TCGA featureCounts output into each cohort folder###########
##############################################################
library("filesstrings")
library(tidyverse)

#dir.create("My_directory")
#file.create("My_file.txt")
#file.move("My_file.txt", "My_directory")

input<-"TCGA.metainfo.txt"
meta<-read.delim(input,header = TRUE,stringsAsFactors = FALSE,sep = "\t")
meta$fileid<-paste0(meta$FILE_ID,".bam.txt")
##Solid Tissue Normal--normal,others--tumor##
for (i in 1:nrow(meta)){
  meta$group[i]<-ifelse(meta$TUMOR_TYPE[i]=="Solid Tissue Normal","Normal","Tumor")
}

##cohort name##
for (i in 1:nrow(meta)){
  meta$cohort[i]<-strsplit(meta$PROJECT[i],"-")[[1]][2]
}

input2<-"tmp1_bam.txt"
fcfile<-read.delim(input2,header = FALSE,stringsAsFactors = FALSE,sep = "\t")

##join fcfile and meta
fcfile.join<-left_join(fcfile,meta,by=c("V1"="fileid"))
fcfile.join$newfileid<-paste0(fcfile.join$FILE_ID,".",fcfile.join$group,".",fcfile.join$cohort,".txt")
fcfile.join$newpath<-paste0("/data/ryang/yren/pacbio/tcga/",fcfile.join$cohort,"/",fcfile.join$group)
fcfile.join$oldpath<-paste0("/data/ryang/yren/pacbio/tcga/tmp_1/",fcfile.join$V1)
##rename file##
file.rename(from, to)

fcfile.join.cmd<-data.frame(cmd="filesstrings::file.rename(",from=fcfile.join$V1,to=fcfile.join$newfileid,end=")")

##file move##
#dir.create("My_directory")
file.move("My_file.txt", "My_directory")
fcfile.join.movecmd<-data.frame(cmd="filesstrings::file.move(",from=paste0('"',fcfile.join$oldpath,'",'),to=paste0('"',fcfile.join$newpath,'"'),end=")")

write.table(fcfile.join.movecmd,"fcfile.join.movecmd.tmp1.txt",col.names = FALSE,row.names = FALSE,quote = FALSE,sep = " ")
