######################2019-09-06#######################
##prad log2(x+1) data from xena match with meta.info###
#######################################################

#args<-commandArgs(TRUE)

meta<-read.delim("TCGA.metainfo.txt",header=TRUE,stringsAsFactors = FALSE,check.names = FALSE)
meta.prad<-meta[meta$PROJECT=="TCGA-PRAD",c(1,4)]
meta.prad$sample<-substr(meta.prad$ALIQUOT_BARCODE,1,15)

prad<-read.delim("HiSeqV2_prad.sub.txt",header=TRUE,stringsAsFactors = FALSE,check.names = FALSE)

#View(prad)
rownames(prad)<-prad$sample
prad$sample<-NULL

##transpose##

pradT<-as.data.frame(t(prad))
pradT$sample<-rownames(pradT)

##join##
pradT.join<-dplyr::full_join(pradT,meta.prad,by="sample")
pradT.join<-pradT.join[!duplicated(pradT.join$sample),]

pradT.join$type<-ifelse(pradT.join$TUMOR_TYPE=="Solid Tissue Normal",0,1)
pradT.join<-pradT.join%>%arrange(type)
rownames(pradT.join)<-pradT.join$sample

ex<-c("sample","ALIQUOT_BARCODE","TUMOR_TYPE","type")
pradT.join.sub<-dplyr::select(pradT.join,-ex)

##transpose back##
pradT.join.subT<-as.data.frame(t(pradT.join.sub))

pradT.join.subT$sample<-rownames(pradT.join.subT)
pradT.join.subT<-pradT.join.subT[,c(ncol(pradT.join.subT),1:ncol(pradT.join.subT)-1)]

write.table(pradT.join.subT,"TCGA_prad_RSEM_normcount.txt",quote=FALSE,row.names = FALSE,col.names = TRUE,sep = "\t")





