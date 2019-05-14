###########################2019-05-14############################
##how many rows for each unique id is called exitron load########
##separate patient into two groups based on exitron load#########
##top 25% and bottom 25%#########################################
#################################################################

##load package##
library(tidyverse)
library(data.table)
library("edgeR")
library("limma")

##read data##
input<-"OV.exitron"
load<-read.delim(input,header=FALSE,stringsAsFactors = FALSE)

##calculate exitron load##
loadstat<-load%>%
  group_by(V16)%>%
  summarize(load=n())

##select top and bottom 25%##
load.q75<-quantile(loadstat$load,probs = 0.75,na.rm = TRUE) #388
load.q25<-quantile(loadstat$load,probs = 0.25,na.rm = TRUE) #192

loadstat.sub<-loadstat[loadstat$load>=load.q75|loadstat$load<=load.q25,]

##top 25->high, bottom 25->low##
for (i in 1:nrow(loadstat.sub)){
  loadstat.sub$group[i]<-ifelse(loadstat.sub$load[i]>=388,"high","low")
}

highload<-loadstat.sub$V16[loadstat.sub$group=="high"]
lowload<-loadstat.sub$V16[loadstat.sub$group=="low"]

##read gene count file## 
input2<-"OV.txt"
count<-read.delim(input2,header = TRUE,row.names="GeneID",stringsAsFactors = FALSE,check.names = FALSE)
countsub<-count[,colnames(count)%in%loadstat.sub$V16|colnames(count)=="Length"]
##DEG analysis##
data <- countsub[,seq(2,ncol(countsub))]
##low first (65) high second (65)
setcolorder(data, c(lowload,highload))
group <- factor(c(rep(1,65),rep(2,65)))

##create DEG object##
y <- DGEList(counts=data, group=group)
y$genes <- data.frame(Length=countsub$Length)
rownames(y$genes) <- rownames(y$counts)
  
# filtering out low expressed genes
#  the minimum number of samples in each group is 65, over here.
keep <- rowSums(cpm(y)>1) >= 65
table(keep)
y <- y[keep, keep.lib.sizes=FALSE]

design <- model.matrix(~0+group)
rownames(design) <- colnames(y)
design

# normalization by the library sizes
y <- calcNormFactors(y)
y$samples

write.table(rpkm(y), file=paste0('OV.rpkm','.exitronload.txt'), sep='\t', row.names = TRUE, quote = FALSE)

# normalise the read counts with 'voom' function
v <- voom(y,design,plot = TRUE)
# extract the normalised read counts
counts.voom <- v$E
# save normalised expression data into output dir
write.table(counts.voom,file="counts.voom.txt",row.names=T,quote=F,sep="\t")

# fit linear model for each gene given a series of libraries
fit <- lmFit(v, design)
# construct the contrast matrix corresponding to specified contrasts of a set of parameters
matrix.2vs1 <- makeContrasts(group2-group1,levels=design)
# compute estimated coefficients and standard errors for a given set of contrasts
fit.2vs1 <- contrasts.fit(fit, matrix.2vs1)

# compute moderated t-statistics of differential expression by empirical Bayes moderation of the standard errors
fit.2vs1 <- eBayes(fit.2vs1)
options(digits=3)

#colnames(fit$coefficients)

summary(decideTests(fit.2vs1, p.value=0.05,lfc=1))
#group2 - group1  (high vs low)
#Down               295
#NotSig           13323
#Up                 259 (upregulated)

num = length(fit.2vs1$genes$Length)
degs.2vs1 <- topTable(fit.2vs1, coef="group2 - group1", confint=TRUE, number = num)
write.table(degs.2vs1, file=paste0('OV.highvslow','exitronload.degs.txt'), sep='\t',row.names = TRUE, quote = FALSE)

##read in splicing factor list##
input3<-"splicing_factor_genes.txt"
splicing<-read.delim(input3,header = TRUE,stringsAsFactors = FALSE,check.names = FALSE) #"Gene Symbol"

degs.2vs1$gene<-rownames(degs.2vs1)

degs.2vs1.splice<-full_join(degs.2vs1,splicing,by=c("gene"="Gene Symbol"))
degs.2vs1.splice.sub<-na.omit(degs.2vs1.splice)

degs.2vs1.splice.sub.pos<-degs.2vs1.splice.sub[degs.2vs1.splice.sub$logFC>0&degs.2vs1.splice.sub$adj.P.Val<=0.05,] #64
degs.2vs1.splice.sub.neg<-degs.2vs1.splice.sub[degs.2vs1.splice.sub$logFC<0&degs.2vs1.splice.sub$adj.P.Val<=0.05,] #120


##From that we can select two list of genes to test for overlap significance using my homemade approach (fisher exact test, fisher.test):

geneset1 <- degs.2vs1$gene[degs.2vs1$logFC>0&degs.2vs1$adj.P.Val<=0.05] #3224
geneset2 <- splicing$`Gene Symbol` #404

universe <- length(
  unique(c(degs.2vs1$gene,geneset2))      #13901
)

common <- length(
  intersect(
    unique(geneset1),               #64
    unique(geneset2)
  )
)


mat <- matrix(
  c(
    universe - length(union(geneset1, geneset2)),
    length(setdiff(geneset1, geneset2)),
    length(setdiff(geneset2, geneset1)),
    length(intersect(geneset1, geneset2))
  ),
  nrow=2
)

fr <- fisher.test(mat, alternative="greater")
fr
 
                

##geneset1 <- degs.2vs1$gene[degs.2vs1$logFC<0&degs.2vs1$adj.P.Val<=0.05] 
##	Fisher's Exact Test for Count Data
#data:  mat
#p-value = 2e-07
#alternative hypothesis: true odds ratio is greater than 1
#95 percent confidence interval:
#  1.5 Inf
#sample estimates:
#  odds ratio 
#1.81 
