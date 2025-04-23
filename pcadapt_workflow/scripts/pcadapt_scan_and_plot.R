### ENVIRONMENT
library(pcadapt)
library(tidyverse)

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
results_output_prefix <- snakemake@params$results_output_prefix
pcadapt.k_param <- snakemake@params$k_params
plot.title <- snakemake@params$plot_title

# create output dirs if not already present
ifelse(!dir.exists(file.path(raw_output_prefix)),
        dir.create(file.path(raw_output_prefix)),
        "Directory Exists")

# pcadapt variables
pcadapt.file <- read.pcadapt(input = input_bed, type = "bed")

# plotting output-variables
output.scree_pca <- paste0(raw_output_prefix, "scree_pca.jpg")
output.pval <- paste0(raw_output_prefix, "manhattan_pval.jpg")
output.padj <- paste0(raw_output_prefix, "manhattan_padj.jpg")

### EXECUTION
# What do I want to include in the analysis?
# - scree and pca plot
# - 2 manhattan plots of physical SNP location with p and padj
# - R object of pcadapt object and list of SNPs with p and padj values

# Produce a pcadapt object with 20 or so k-mers to plot pca and scree
pcadapt.diffK <- pcadapt(input = pcadapt.file, K = 20)

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
  output.scree_pca,
  width = 25,
  height = 10,
  units = "in",
  res = 100
)
ggarrange(plot.scree, plot.pca, nrow = 2, ncol = 2)
dev.off()

# create pcadapt object for analysis with k_params
pcadapt.scan <- pcadapt(pcadapt.file, K = pcadapt.k_param)

# create dataframe of SNP positions, p and padj values
pcadapt.df <- get_bim_vars(paste0(input_prefix, ".bim"))
pcadapt.df$pvalue <- pcadapt.scan$pvalues
pcadapt.df$padj <- p.adjust(pcadapt.scan$pvalues, method = "bonferroni")

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
