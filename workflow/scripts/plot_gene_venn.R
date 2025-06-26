### ENVIRONMENT
library(tidyverse)
library(ggVennDiagram)

### VARIABLES
input_shared_genes_file <- snakemake@input[[1]]
input_fr_genes_file <- snakemake@input[[2]]
input_sr_genes_file <- snakemake@input[[3]]
output_venn_file <- snakemake@output[[1]]

### EXECUTION
venn_df <- list(
  shared_genes = read.delim(file = input_shared_genes_file)[[9]],
  fr_genes = read.delim(file = input_fr_genes_file)[[9]],
  sr_genes = read.delim(file = input_sr_genes_file)[[1]]
)

attributes(venn_df) <- list(names = names(venn_df),
  row.names = 1:max(length(venn_df$shared_genes), length(venn_df$fr_genes), length(venn_df$fr_genes)),
  class = 'data.frame'
)

jpeg(
  output_venn_file,
  width = 10,
  height = 10,
  units = "in",
  res = 100
)
ggVennDiagram(venn_df)
dev.off()
