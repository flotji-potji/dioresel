### ENVIRONMENT
if (!require("topGO")) {
  install.packages("BiocManager")
  BiocManager::install("topGO")
}
library(tidyverse)
library(topGO)

### VARIABLES
shared_genes_file <- snakemake@input$shared
fr_genes_file <- snakemake@input$fr_genes
sr_genes_file <- snakemake@input$sr_genes

go_df_file <- snakemake@output[[1]]
go_table_file <- snakemake@output[[2]]

### EXECUTION
shared_genes <- read.delim(file = shared_genes_file, header = FALSE)
fr_genes <- read.delim(file = fr_genes_file, header = FALSE)
sr_genes <- read.delim(file = sr_genes_file, header = FALSE)

gene_annotation <- data.frame(
  category = c(rep("shared_genes", dim(shared_genes)[1]),
               rep("fr_genes", dim(fr_genes)[1]),
               rep("sr_genes", dim(sr_genes)[1])),
  gene_name = c(shared_genes$V9, fr_genes$V9, sr_genes$V1),
  gene_desc = c(shared_genes$V10, fr_genes$V10, sr_genes$V10),
  GO = c(shared_genes$V11, fr_genes$V11, sr_genes$V11)
)

Unfold <- gene_annotation %>%
  mutate(GO = strsplit(as.character(GO), ",")) %>%
  unnest(GO)
geneID2GO <- Unfold %>% split(x = .$GO, f = .$gene_name)
geneNames <- names(geneID2GO)

GOdata = NULL
myInterestingGenes = NULL
geneList = NULL
for(i in 1:length(unique(all.markers$cluster)))
{
    #identify the GOI from the DE list
    myInterestingGenes[[i]] <- all.markers$gene[all.markers$cluster==unique(all.markers$cluster)[i]] #list of genes you want to perform GO enrichment for
    geneList[[i]] <- factor(as.integer(geneNames %in% myInterestingGenes[[i]]))
    names(geneList[[i]]) <- geneNames
    GOdata[[i]] <- new("topGOdata",ontology = "BP", allGenes = geneList[[i]],annot = annFUN.gene2GO, gene2GO = geneID2GO)
}

# filter and generate figures.
p=NULL
resultFis=NULL
for(i in 1:length(unique(all.markers$cluster)))
{
    #run the test...
    resultFis[[i]] <- runTest(GOdata[[i]], algorithm = "classic", statistic = "fisher") 
    pvalFis <- score(resultFis[[i]])
    #filter for only >0.05
    pvalFis = pvalFis[pvalFis>=0.05]
    allRes_intgenes<- GenTable(GOdata[[i]], pvalues = resultFis[[i]], orderBy = "pvalues", topNodes=30)
    allRes_intgenes$pvalues<-as.numeric(allRes_intgenes$pvalues)
    #convert NA values to zero; only if p value is so small it cannot be displayed in r
    allRes_intgenes[is.na(allRes_intgenes)]<-0.00000001
    
    #plot GOenrichment 
    
    colCP=c("darkred")
    p[[i]]=
    ggplot2::ggplot(allRes_intgenes, ggplot2::aes(x=(reorder(Term,(-log10(pvalues)))), y=(-log10(pvalues)))) +
    ggplot2::stat_summary(geom = "bar", fun = mean, position = "dodge",
                            col=colCP,fill=colCP) +
    ggplot2::coord_flip()+
    ggplot2::xlab("Biological Process") +
    ggplot2::ylab("Enrichment -log10 p-value") +
    ggplot2::labs(title="GO enrichment",subtitle=levels(all.markers$cluster)[i])
}
