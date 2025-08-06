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
input_files <- snakemake@input
output_picmin_results_file <- snakemake@output[[1]]

numReps <- snakemake@params$num_reps
num_genes <- 27000

### EXECUTION
df_list <- list()
for (i in seq(length(input_files))) {
  load(file = input_files[[i]])
  df_list[[i]] <- lineage_out
}

all_lins <- df_list %>% reduce(full_join, by = "window")
all_lins_p <- all_lins[ , !(names(all_lins) %in% c("window"))]
rownames(all_lins_p) <- all_lins$window

count <- 0
results <- list()
nLins <- length(df_list)
missingDataLevels <- seq(3, nLins)

for (n in missingDataLevels){
  count <- count + 1
  emp_p_null_dat <- t(replicate(40000,
                                PicMin:::GenerateNullData(1.0, n, 0.5,
                                                          3, num_genes)))
  emp_p_null_dat_unscaled <- t(apply(emp_p_null_dat, 1,
                                     PicMin:::orderStatsPValues))
  # Use those p-values to construct the correlation matrix
  null_pMax_cor_unscaled <- cor(emp_p_null_dat_unscaled)


  # Screen out gene with no evidence for adaptation
  lins_p_n <-  as.matrix(all_lins_p[rowSums(is.na(all_lins_p)) == nLins - n, ])

  if (dim(lins_p_n)[1] == 0) {
    next
  }
  res_p <- rep(-1, nrow(lins_p_n))
  res_n <- rep(-1, nrow(lins_p_n))

  for (i in seq(nrow(lins_p_n))) {
    test_result <- PicMin:::PicMin(na.omit(lins_p_n[i, ]),
                                   null_pMax_cor_unscaled, numReps = numReps)
    res_p[i] <- test_result$p
    res_n[i] <- test_result$config_est
  }
  results[[count]] <- data.frame(numLin = n,
                                 p = res_p,
                                 q = p.adjust(res_p, method = "fdr"),
                                 n_est = res_n,
                                 locus = row.names(lins_p_n))

}

picMin_results <- do.call(rbind, results)

picMin_results$pooled_q <- p.adjust(picMin_results$p, method = "fdr")

picMin_results <- cbind(picMin_results,
  read.csv(text = picMin_results$locus, header = FALSE,
           sep = "_",
           col.names = c('redundan', 'scaffold', 'start'))
)

save(picMin_results, file = output_picmin_results_file)