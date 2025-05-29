### ENVIRONMENT
library(tidyverse)
if (!require(PicMin)) {
  install.packages("devtools", repos = "https://cloud.r-project.org")
  library(devtools)
  install_github("TBooker/PicMin", force = TRUE)
  install.packages("poolr", repos = "https://cloud.r-project.org")
  library(PicMin)
  library(poolr)
} else {
  library(PicMin)
  library(poolr)
}
options(error = quote({
  dump.frames(to.file=T, dumpto='last.dump')
  load('last.dump.rda')
  print(last.dump)
  q()
}))

### VARIABLES
input_bed <- snakemake@input[[1]]
output_picmin_df_file <- snakemake@output[[1]]

pval_small <- snakemake@params$pval_small

### EXECUTION
lineage_df <- read.delim(input_bed, header = FALSE)
lineage_df$emp_p <- PicMin:::EmpiricalPs(lineage_df$V4,
                                         large_i_small_p = pval_small)
colnames(lineage_df)[1] <- "scaff"
colnames(lineage_df)[2] <- "start"
lineage_df$name <- PicMin:::get_names(lineage_df)
lineage_out <- PicMin:::min_lin(lineage_df,
                                paste0("sp_", snakemake@wildcards$sp2))
save(lineage_out, file = output_picmin_df_file)