######################2019-05-15#################
##apply overlaptest function to *.exitron and ###
##*.txt to get final output for plot#############
#################################################

source("exitronload.splicefunc.r")
cohort<-read.delim("cohort.exitron",header = FALSE,stringsAsFactors = FALSE)
count<-read.delim("count.txt",header = FALSE,stringsAsFactors = FALSE)

cohort$V2<-str_replace_all(cohort$V1,".exitron","")
count$V2<-str_replace_all(count$V1,".txt","")

cohort.count<-full_join(cohort,count,by="V2")
cohort.count<-na.omit(cohort.count)

test<-mapply(overlaptest,cohort.count$V1.x,cohort.count$V1.y)
testT<-t(test)
sink("testT.txt")
testT<-as.data.frame(testT)
print(testT)
sink()
#write.table(testT,"exitronload.splice.pval4plot.txt",row.names = TRUE,col.names = TRUE,quote = FALSE,sep="\t")