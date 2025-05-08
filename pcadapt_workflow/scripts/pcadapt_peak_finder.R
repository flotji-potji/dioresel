### ENVIRONMENT
library(pcadapt)
library(tidyverse)

### FUNCTIONS
manhattan_data_frame <- function(df, snps_of_interest) {
  colnames(df) <- c("CHR", "BP", "REF", "ALT", "pval", "padj", "snpID")

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
    mutate(BPcum=BP+tot) %>% 
    mutate( is_highlight=ifelse(snpID %in% snps_of_interest$snpID, "yes", "no"))

  return(don)
}

plot_manhattan <- function(df, pvalue, title="") {
  axisdf = don %>%
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
input_peaks <- snakemake@input[[1]]
output_file <- snakemake@output[[1]]
pcadapt.df.file <- snakemake@params$pcadapt_res
raw_output_prefix <- snakemake@params$raw_output_prefix
plot_title <- snakemake@params$plot_title
pval_threshold <- snakemake@params$pval_threshold
pval_type <- snakemake@params$pval_type

output.peak_find <- paste0(raw_output_prefix, "_manhattan_plot.jpg")

### EXECUTION
pcadapt.df <- read.delim(pcadapt.df.file)
colnames(pcadapt.df) <- c("chrom", "position", "ref", "alt", "pvalue", "padj")
peaks <- read.delim(input_peaks)

pcadapt.df <- pcadapt.df %>%
  mutate(snpID = row_number())

peaks <- peaks %>%
  filter(!is.na(max)) %>%
  separate(range, into = c("rangeLower", "rangeUpper"),
           sep = "-", convert = TRUE)

peaks_df_intersect <- peaks %>%
  inner_join(pcadapt.df, by = "chrom") %>%
  filter(position >= rangeLower & position <= rangeUpper)

snps_of_interest <- peaks_df_intersect %>%
  filter(pvalue < pval_threshold) %>%
  filter(range.1 < quantile(peaks_df_intersect$range.1, 0.75)) %>%
  filter(spacing < quantile(peaks_df_intersect$spacing, 0.75))

# plot manhattan plot of all scaffolds

# Bonferroni corrected manhattan plot
don <- manhattan_data_frame(pcadapt.df, snps_of_interest)

jpeg(
  output.peak_find,
  width = 25,
  height = 10,
  units = "in",
  res = 100
)
plot_manhattan(don, pval_type, 
               paste("pcadapt -", plot_title, ": peak highlights")) +
  geom_hline(yintercept = -log10(0.05), linetype = "dotted") +
  geom_hline(yintercept = -log10(0.01)) +
  geom_point(data=subset(don, is_highlight=="yes"), color="red", size=2) +
  theme(text = element_text(size = 25),
        axis.text = element_text(color="black", size=25),
        axis.title = element_text(color = "black", size=25, face = "bold"))
dev.off()

# store results as tsv and R object file
write.table(snps_of_interest, output_file,
            sep = "\t", row.names = FALSE, quote = FALSE)