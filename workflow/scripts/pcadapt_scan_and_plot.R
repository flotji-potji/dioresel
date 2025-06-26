### ENVIRONMENT
library(tidyverse)
if (!require(pcadapt)) {
  install.packages("pcadapt", repos = "https://cran.wu.ac.at/")
  library(pcadapt)
} else {
  library(pcadapt)
}
options(error = quote({
  dump.frames(to.file=T, dumpto='last.dump')
  load('last.dump.rda')
  print(last.dump)
  q()
}))

### FUNCTIONS
get_sample_names <- function(file) {
  samples.array <- read_delim(file, col_select=1, 
                              col_names = c("fem", "mal"))$fem
  split_list <- strsplit(samples.array, split = "(?<=[a-z])\\s*(?=[0-9A-Z])",
                        perl = TRUE)
  res <- c()
  for (arr in split_list) {
    res <- append(res, arr[1])
  }
  return(res)
}

get_bim_vars <- function(file) {
  bim.df <- read_delim(file, col_names = c("chr", "varID", "genDist", 
                                            "position", "ref", "alt"))
  bim.df <- bim.df[, !names(bim.df) %in% c("varID", "genDist")]
  return(bim.df)
}

manhattan_data_frame <- function(df) {
  colnames(df) <- c("CHR", "BP", "REF", "ALT", "pval", "padj")
  don <- df %>% 
    # Compute chromosome size
    group_by(CHR) %>% 
    summarise(chr_len=max(BP)) %>% 
    # Calculate cumulative position of each chromosome
    mutate(tot=cumsum(chr_len)-chr_len) %>%
    select(-chr_len) %>%
    # Add this info to the initial dataset
    left_join(df, ., by=c("CHR"="CHR")) %>%
    # Add a cumulative position of each SNP
    arrange(CHR, BP) %>%
    mutate(BPcum=BP+tot)

  return(don)
}

plot_manhattan <- function(df, pvalue, title="") {
  axisdf = df %>%
    group_by(CHR) %>%
    summarize(center=(max(BPcum) + min(BPcum)) / 2)

  ggplot(df, aes(x=BPcum, y=-log10(df[[pvalue]]))) +
    # Show all points
    geom_point(aes(color=as.factor(CHR)), alpha=0.8, size=1.3) +
    scale_color_manual(values = rep(c("black", "grey"), 31)) +
    # custom X axis:
    scale_x_continuous(label = axisdf$CHR, breaks= axisdf$center ) +
    scale_y_continuous(expand = c(0, 0) ) +
    # custom labels
    labs(y = "-log10(p)",
         x = "Scaffold",
         title = title) +
    # Custom the theme:
    theme_bw() +
    theme( 
      legend.position="none",
      panel.border = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)
    )
}

### VARIABLES
input_bed <- snakemake@input[[1]]
input_prefix <- snakemake@params$input_prefix
output_file <- snakemake@output[[1]]
raw_output_prefix <- snakemake@params$raw_output_prefix
pcadapt.k_param <- snakemake@params$k_param
plot.title <- snakemake@params$plot_title
genome_size <- snakemake@params$vie_genome_size

# create output dirs if not already present
ifelse(!dir.exists(file.path(raw_output_prefix)),
        dir.create(file.path(raw_output_prefix)),
        "Directory Exists")

# pcadapt variables
pcadapt.file <- read.pcadapt(input = input_bed, type = "bed")
window_size <- 10000
ld_thin.threshold <- 0.1
output.pcadapt_object <- paste0(raw_output_prefix, "pcadapt_object.Robject")

# plotting output-variables
output.scree <- paste0(raw_output_prefix, "plot_scree.jpg")
output.pca <- paste0(raw_output_prefix, "plot_pca.jpg")
output.pval <- paste0(raw_output_prefix, "plot_manhattan_pval.jpg")
output.padj <- paste0(raw_output_prefix, "plot_manhattan_padj.jpg")

### EXECUTION
# What do I want to include in the analysis?
# - scree and pca plot
# - 2 manhattan plots of physical SNP location with p and padj
# - R object of pcadapt object and list of SNPs with p and padj values

# create dataframe of SNP positions 
pcadapt.df <- get_bim_vars(paste0(input_prefix, ".bim"))
num_total_snps <- length(pcadapt.df)

# Produce a pcadapt object with 20 or so k-mers to plot pca and scree
if (snakemake@params$ld_thin == "yes") {
  window_num_snps <- (num_total_snps / genome_size) * window_size
  pcadapt.diffK <- pcadapt(input = pcadapt.file, K = 20,
                           LD.clumping = list(size = window_num_snps,
                                              thr = ld_thin.threshold))
} else {
  print("in pcadapt")
  pcadapt.diffK <- pcadapt(input = pcadapt.file, K = 20)
  print("out pcadapt")
}


# load sample name file
population.list <- get_sample_names(paste0(input_prefix, ".nosex"))

# plot distribution of PC share
plot.scree <- plot(pcadapt.diffK, option = "screeplot")
# ... and the first to PCs
plot.pca <- plot(pcadapt.diffK,
                 option = "scores", i = 1, j = 2,
                 pop = population.list)

# save plots contained in one jpg
jpeg(
  output.scree,
  width = 25,
  height = 10,
  units = "in",
  res = 100
)
plot.scree
dev.off()

jpeg(
  output.pca,
  width = 25,
  height = 10,
  units = "in",
  res = 100
)
plot.pca
dev.off()

# create pcadapt object for analysis with k_params
if (snakemake@params$ld_thin == "yes") {
  window_num_snps <- (num_total_snps / genome_size) * window_size
  pcadapt.scan <- pcadapt(input = pcadapt.file, K = pcadapt.k_param,
                           LD.clumping = list(size = window_num_snps,
                                              thr = ld_thin.threshold))
} else {
  pcadapt.scan <- pcadapt(input = pcadapt.file, K = pcadapt.k_param)
}

# add p and padj values to SNP-position dataframe
pcadapt.df$pvalue <- pcadapt.scan$pvalues
pcadapt.df$padj <- p.adjust(pcadapt.scan$pvalues, method = "fdr")

# plot pcadapt results as manhattan plots
don <- manhattan_data_frame(pcadapt.df)

# manhattan plot of pcadapt p-values
jpeg(
  output.pval,
  width = 25,
  height = 10,
  units = "in",
  res = 100
)
plot_manhattan(don, "pval", paste("pcadapt - ", plot.title, ": p-value")) +
  geom_hline(yintercept = -log10(0.05), linetype = "dotted") +
  geom_hline(yintercept = -log10(0.01)) +
  theme(text = element_text(size = 15),
        axis.text = element_text(color="black", size=14),
        axis.title = element_text(color = "black", size=15, face = "bold"))
dev.off()

# manhattan plot of p-adjusted
jpeg(
  output.padj,
  width = 25,
  height = 10,
  units = "in",
  res = 100
)
plot_manhattan(don, "padj", paste("pcadapt - ", plot.title, ": p-adj")) +
  geom_hline(yintercept = -log10(0.05), linetype = "dotted") +
  geom_hline(yintercept = -log10(0.01)) +
  theme(text = element_text(size = 25),
        axis.text = element_text(color="black", size=25),
        axis.title = element_text(color = "black", size=25, face = "bold"))
dev.off()

# store pcadapt results as tsv file
write.table(pcadapt.df, output_file,
            sep = "\t", row.names = FALSE, quote = FALSE)

save(pcadapt.scan, file = output.pcadapt_object)
