### ENVIRONMENT
library(tidyverse)

### FUNCTIONS

### VARIABLES
load(file = snakemake@input[[1]])
MKT <- gene.summary

output.mkt_df <- snakemake@output[[1]]
output.mkt_tsv <- snakemake@output[[2]]
output_prefix <- snakemake@params$output_prefix

### EXECUTION

# calculate MK ratios
MKT$pN.pS <- MKT$pN/MKT$pS
MKT$dN.dS <- MKT$dN/MKT$dS

# create new column for p-values
MKT$fisher.test.P <- 99
for (i in 1:nrow(MKT)) {
  # calculate fisher exact test and copy p-value for every contig
  MKT$fisher.test.P[i] <- fisher.test(matrix(as.numeric(MKT[i, c(2, 4, 3, 5)]),
                                             ncol = 2))$p.value
  # this lines assigns an NA to all p-values that are meaningless,
  # because the contingency table was incomplete
  if ((MKT$pN[i] == 0 && MKT$dN[i] == 0) ||
        (MKT$pS[i] == 0 && MKT$dS[i] == 0) ||
        (MKT$pS[i] == 0 && MKT$pN[i] == 0) ||
        (MKT$dS[i] == 0 && MKT$dN[i] == 0)) {
    MKT$fisher.test.P[i] <- NA
  }
  # only use cases where total number of SNPs is higher than or equal to 3
  if (sum(as.numeric(MKT[i, c(2, 4, 3, 5)])) < 3) {
    MKT$fisher.test.P[i] <- NA
  }
}

# multiple hypothesis testing
# correct p-values using Benjamini & Hochberg (1985) FDR
MKT$p.adj <- p.adjust(MKT$fisher.test.P, method = "BH")

# store MKT-dataframe as R object and tsv
ifelse(!dir.exists(file.path(output_prefix)),
        dir.create(file.path(output_prefix)),
        "Directory Exists")

save(MKT, file = output.mkt_df)

write.table(MKT, output.mkt_tsv,
            sep = "\t", row.names = FALSE, quote = FALSE)
