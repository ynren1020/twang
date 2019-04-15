library(ggplot2)
library(ggrepel)
library(sets)


###############################################################



luad<- read.csv('./LUAD.scatter.txt', header=TRUE, sep="\t", row.names="gene")

##cancer genes##
genes <- as.vector(read.csv('./LUAD_genes_in_cosmic.txt', header = FALSE)$V1); #29

cancer.genes <- as.vector(read.csv('./cancer_genes.hc_list.txt', header = TRUE,sep="\t")$GeneName);#831

gene.list <- c(genes, cancer.genes);

##cancer genes in area I##
big <- as.set(rownames(subset(luad, exitron_freq >2 & mutation_freq >2 )))
gene.list.I <- unlist(set_intersection(big, as.set(gene.list))); #10
#{"BCL9L", "COL1A1", "DHX9", "ERBB3", "FLNA", "KMT2D", "MUC4", "MUC6", "MYH9", "NUMA1"}
##cancer genes in area II or IV##
bigII<-as.set(rownames(subset(luad, exitron_freq <2 & mutation_freq >10 )))
gene.list.II<-c(unlist(set_intersection(bigII, as.set(genes))),"NOTCH1","RB1") # "CSMD3" "EGFR"  "KRAS"  "PTPRD" "TP53" 
##add another two genes##"NOTCH1""RB1" 
##areaIV##
bigIV<-as.set(rownames(subset(luad, exitron_freq >2 & exitron_freq <50&mutation_freq <2 )))
gene.list.IV<-unlist(set_intersection(bigIV, as.set(gene.list))) 
gene.list.IV<-c("CTNND1","FUS","EWSR1","MUC1","EEF2")





luad$label <- ifelse(rownames(luad) %in% gene.list.II | rownames(luad) %in% gene.list.IV, rownames(luad), NA);


ggplot(luad, aes(exitron_freq, mutation_freq)) +
    #scale_y_continuous(trans='log10') +
    #scale_x_continuous(trans='log10') +
    geom_vline(xintercept = 2, linetype='dashed', color='grey') +
    geom_hline(yintercept = 2, linetype='dashed', color='grey') +
    geom_point(color = ifelse(is.na(luad$label), "#3987cc", "red"),size = ifelse(is.na(luad$label), 1.0, 2.0)) +
    
    geom_label_repel(
        aes(exitron_freq, mutation_freq, label=luad$label),
        arrow = arrow(length = unit(0.0001, "npc"), type = "closed", ends = "first"),
        force = 8,
        size = 5,
        max.iter=1000
    )+
    labs(title="",
         x="Exitron splicing frequency (%)", y = "Mutation frequency (%)")+
    xlim(0, 50)+
    theme_classic()+
    theme(axis.text.x = element_text(size=16),
          axis.text.y = element_text(size=16),
          axis.title=element_text(size=18))
################################################################
