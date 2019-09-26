#####################################2019-09-25#########################################
########################mich:MICH.gene.wilcox.pvalue.sig.txt############################
########################tale:TALE.gene.wilcox.pvalue.sig.txt############################
##prad:prad.tumor.vs.normal.degs.normalfiltered.loose_updatesu2c101_groupCorrect.txt####
########################################################################################

library(dplyr)
library(purrr)
library(tidyverse)

mich<-read.delim("MICH.gene.wilcox.pvalue.sig.txt",header = TRUE,stringsAsFactors = FALSE,check.names = FALSE)
tale<-read.delim("TALE.gene.wilcox.pvalue.sig.txt",header = TRUE,stringsAsFactors = FALSE,check.names = FALSE)
prad<-read.delim("prad.tumor.vs.normal.degs.normalfiltered.loose_updatesu2c101_groupCorrect.txt",header = TRUE,stringsAsFactors = FALSE,check.names = FALSE)

##filter by significance##
#mich.sub<-mich[mich$primary.vs.met.p.adj<=0.05&mich$normal.vs.met.p.adj<=0.05&mich$normal.vs.primary.p.adj<=0.05,] #3380
#tale.sub<-tale[tale$primary.vs.met.p.adj<=0.05&tale$normal.vs.met.p.adj<=0.05&tale$normal.vs.primary.p.adj<=0.05,] #397
prad.sub<-prad[prad$adj.P.Val<=0.05,] #22565

##join##
#mich.sub<-mich%>%rename("gene"="gene.mich")
#tale.sub<-tale%>%rename("gene"="gene.tale")
#prad.sub<-prad.sub%>%rename("geneid"="gene.prad")

mich$gene.mich<-mich$gene
tale$gene.tale<-tale$gene
prad.sub$gene.prad<-prad.sub$geneid

mich.sub<-mich[,c("gene","gene.mich")]
tale.sub<-tale[,c("gene","gene.tale")]
prad.sub<-prad.sub[,c("geneid","gene.prad")]

#list(mich.sub.sub, tale.sub.sub, prad.sub.sub) %>% reduce(full_join, by = "gene")
two.join<-full_join(mich.sub,tale.sub,by="gene")
three.join<-full_join(two.join,prad.sub,by=c("gene"="geneid"))

##at least two studies share common genes##
keep<-rowSums(is.na(three.join))<=1
three.join.sub<-three.join[keep,]

##all three share common genes##
keep<-rowSums(is.na(three.join))<1
three.join.sub.all<-three.join[keep,]

write.table(three.join.sub,"mich.tale.prad.twomatch.wilcox.txt",quote = FALSE,col.names = TRUE,row.names = FALSE,sep="\t")
write.table(three.join.sub.all,"mich.tale.prad.allmatch.wilcox.txt",quote = FALSE,col.names = TRUE,row.names = FALSE,sep="\t")

##further filter by known gene##
known<-read.delim("prostate_knowngene.txt",header = TRUE,stringsAsFactors = FALSE)
known$gene<-known$gene_known

##may need to convert to ensemble id ##
##genename to geneid##
library("org.Hs.eg.db") # remember to install it if you don't have it already
symbols<-known$gene_known
known$geneid<-mapIds(org.Hs.eg.db, symbols, 'ENSEMBL', 'SYMBOL') #symbol to geneid##

  gene_known      gene geneid
47      MRE11A    MRE11A   ENSG00000020922
  137   C10orf54  C10orf54   ENSG00000107738
  193  C14orf105 C14orf105   ENSG00000100557
  240     GPR123    GPR123   ENSG00000197177
  244       LHFP      LHFP   ENSG00000183722
  260    FAM194A   FAM194A   ENSG00000163645
  278      INADL     INADL   ENSG00000132849
  281   KIAA1984  KIAA1984   ENSG00000213213
  
known$geneid[known$gene_known=="MRE11A"]<-"ENSG00000020922"
known$geneid[known$gene_known=="C10orf54"]<-"ENSG00000107738"
known$geneid[known$gene_known=="C14orf105"]<-"ENSG00000100557"
known$geneid[known$gene_known=="GPR123"]<-"ENSG00000197177"
known$geneid[known$gene_known=="LHFP"]<-"ENSG00000183722"
known$geneid[known$gene_known=="FAM194A"]<-"ENSG00000163645"
known$geneid[known$gene_known=="INADL"]<-"ENSG00000132849"
known$geneid[known$gene_known=="KIAA1984"]<-"ENSG00000213213"
  


symbols<-three.join.sub.all$gene
three.join.sub.all$geneid<-mapIds(org.Hs.eg.db, symbols, 'ENSEMBL', 'SYMBOL') #symbol to geneid##

gene gene.mich gene.tale gene.prad geneid
1823 MEIS3P1   MEIS3P1   MEIS3P1   MEIS3P1   ENSG00000179277

three.join.sub.all$geneid[three.join.sub.all$gene=="MEIS3P1"]<-"ENSG00000179277"

##join three and known, keep unknown genes##
three.join.sub.all.known<-full_join(three.join.sub.all,known,by="geneid")

three.join.sub.all.unknown<-three.join.sub.all.known[is.na(three.join.sub.all.known$gene_known),]

write.table(three.join.sub.all.unknown,"mich.tale.prad.allmatch.unknown.wilcox.txt",quote = FALSE,col.names = TRUE,row.names = FALSE,sep="\t")

##tale specificity genes##
tale.spec<-read.delim("TCGA_TALE_DEGgenes_specificity_all.txt",header = FALSE,stringsAsFactors = FALSE)
tale.spec$geneid<-str_sub(tale.spec$V7,1,15)

##join three unknown with tale.spec##
tale.spec.join<-full_join(tale.spec,three.join.sub.all.unknown,by="geneid")

tale.spec.join.sub<-tale.spec.join[!is.na(tale.spec.join$gene.x),]

#zscore >=2#
tale.spec.join.sub.sub<-tale.spec.join.sub[tale.spec.join.sub$V6>=2,]
tale.spec.join.sub.prad<-tale.spec.join.sub.sub[tale.spec.join.sub.sub$V1=="TCGA-PRAD",]
tale.spec.join.sub.prad<-tale.spec.join.sub.prad[1:12,]

write.table(tale.spec.join.sub.prad,"prad.mich.tale_spec.filtered.gene.txt",quote = FALSE,col.names = TRUE,row.names = FALSE,sep="\t")

#12464 TCGA-PRAD 498 13894 27.9 0.827 2.01 ENSG00000103175.10 ENSG00000103175  WFDC1     WFDC1     WFDC1     WFDC1
#41372 TCGA-PRAD 498 22440 45.1 0.886 2.81 ENSG00000154330.12 ENSG00000154330   PGM5      PGM5      PGM5      PGM5


