### ENVIRONMENT
library(tidyverse)

### VARIABLES
input_picmin_file <- snakemake@input[[1]]
output_bed_file <- snakemake@output[[1]]

## EXECUTION
load(file = input_picmin_file)

output_bed <- data.frame(picMin_results$redundan,
                         picMin_results$scaffold,
                         picMin_results$scaffold + 9999,
                         picMin_results$numLin,
                         picMin_results$n_est,
                         picMin_results$p,
                         picMin_results$q,
                         picMin_results$pooled_q)

write_tsv(output_bed, file = output_bed_file, col_names = FALSE)
