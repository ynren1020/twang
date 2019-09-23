#################################2019-09-23##################################
##try to do differential expression analysis on TALE data to find genes######
#############################################################################

tale.filter<-read.delim("MSKCC_PCa_mRNA_data_sig.txt",header=TRUE,stringsAsFactors = FALSE)
genes<-read.delim("TALE.gene.ttest.pvalue.sig_genename.txt",header = FALSE,stringsAsFactors = FALSE)

tale.filter.sub<-tale.filter[tale.filter$GeneSymbol%in%genes$V1,-c(182:187)]
tale.filter.sub$GeneID<-NULL
rownames(tale.filter.sub)<-tale.filter.sub$GeneSymbol
tale.filter.sub$GeneSymbol<-NULL
tale.filter.sub<-tale.filter.sub[,names(tale.filter.sub)%in%tale.surv.sub.sub$`Sample ID`]

##order columns by column names##
tale.filter.sub.order<-tale.filter.sub[,order(names(tale.filter.sub))] #>=PCA0182 are met sample, 131 primary,19mets

##DEG analysis##
library("edgeR")
#library("sva")
library("R.utils")
library("limma")

group <- factor(c(rep(1,131),rep(2,19)))


y <- DGEList(counts=tale.filter.sub.order, group=group);
#y$genes <- data.frame(Length=2890);
#rownames(y$genes) <- rownames(y$counts);

design <- model.matrix(~0+group);
rownames(design) <- colnames(y)
design

v <- voom(y,design,plot = TRUE)
# extract the normalised read counts
counts.voom <- v$E

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
# group2 - group1
# Down               170
# NotSig           11685
# Up                 646

num = length(fit.2vs1$F.p.value)
degs.2vs1 <- topTable(fit.2vs1, coef="group2 - group1", confint=TRUE, number = num)
degs.2vs1$genes<-rownames(degs.2vs1)
write.table(degs.2vs1,"TALE_DEGgenes.primaryvsMET.txt", sep='\t',row.names = FALSE, col.names = TRUE, quote = FALSE);

##negative logFC## tumor suppressor
degs.2vs1.neg<-degs.2vs1[degs.2vs1$logFC<0&degs.2vs1$adj.P.Val<0.05,] #1240 10
write.table(degs.2vs1.neg,"TALE_DEGgenes.primaryvsMET.neg.potentialSuppressor.txt", sep='\t',row.names = FALSE, col.names = TRUE, quote = FALSE);

##positive logFC## onco gene
degs.2vs1.pos<-degs.2vs1[degs.2vs1$logFC>0&degs.2vs1$adj.P.Val<0.05,] #1528 10
write.table(degs.2vs1.pos,"TALE_DEGgenes.primaryvsMET.neg.potentialOnco.txt", sep='\t',row.names = FALSE, col.names = TRUE, quote = FALSE);

##gene name to ensemble id##
BiocManager::install("org.Hs.eg.db")
library("org.Hs.eg.db") # remember to install it if you don't have it already
symbols<-degs.2vs1$genes
geneid<-mapIds(org.Hs.eg.db, symbols, 'ENSEMBL', 'SYMBOL') #symbol to geneid##
##to entranz##
entrezid<-mapIds(org.Hs.eg.db, symbols, 'ENTREZID', 'SYMBOL')
##gene id to symbols##
#symbols <- mapIds(org.Hs.eg, keys = ensemblsIDS, keytype = "ENSEMBL", column="SYMBOL")

##use biomart to convert gene name to ensemble id##
library(biomaRt)
genes <- degs.2vs1$genes
mart <- useMart("ensembl", dataset="hsapiens_gene_ensembl")
attributes=c('ensembl_gene_id','hgnc_symbol')
G_list<- getBM(attributes=c('ensembl_gene_id','hgnc_symbol'), filters="hgnc_symbol",values=genes,
               mart=mart, uniqueRows=T)

#df.id.name.match<-unique(G_list[,c(1,3)])

df.id.name.match.DEG<-full_join(degs.2vs1,G_list,by=c("genes"="hgnc_symbol")) #one symbol two id ??##
########
geneid.df<-data.frame(name=names(geneid),id=geneid)
df.id.name.match.DEG<-full_join(degs.2vs1,geneid.df,by=c("genes"="name"))

write.table(df.id.name.match.DEG,"TALE_DEGgenes.primaryvsMET.idmatchname.txt", sep='\t',row.names = FALSE, col.names = TRUE, quote = FALSE);

df.id.name.match.DEG.wona<-df.id.name.match.DEG[!is.na(df.id.name.match.DEG$id),]

geneid4spec<-df.id.name.match.DEG.wona$id
write.table(geneid4spec,"TALE_DEGgenes.primaryvsMET.idmatchname.idonly.txt",sep='\t',row.names = FALSE, col.names = TRUE, quote = FALSE)
