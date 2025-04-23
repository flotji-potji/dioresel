#!/usr/bin/env R

library(tidyverse)

# helper functions
plot_manhattan_fst <- function(df, title="") {
  ggplot(df, aes(x=BPcum, y=pval)) +
    
    # Show all points
    geom_point(aes(color=as.factor(CHR)), alpha=0.8, size=1.3) +
    scale_color_manual(values = rep(c("black", "grey"), 31)) +
    
    # custom X axis:
    scale_x_continuous(label = axisdf$CHR, breaks= axisdf$center ) +
    scale_y_continuous(expand = c(0, 0) ) +     # remove space between plot area and x axis
    
    # custom labels
    labs(y = "Fst",
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

fst_val <- function(file) {
  fst_pos <- read_delim(file, col_names = c("chr", "pos", "fst"),
                        skip = 1)
  return(fst_pos)
}

# main part
# find fst files and sort
dir <- "/lisc/scratch/botany/frschmidt/fst/"
setwd(dir)

file <- "raw_data/vcftools/fst_step_window/vieref_fla_spn.10000.step.windowed.weir.fst"
fst_vals <- read_delim(file, col_names = c("chr", "bin_start", "bin_end", 
                                           "n_variants", "wFst", "mFst"), 
                       skip = 1)

# combine results into data frame
df <- data.frame(fst_vals$chr, 
                 fst_vals$bin_end - ((fst_vals$bin_end-fst_vals$bin_start)/2), 
                 fst_vals$mFst)
names(df) <- c("CHR", "BP", "pval")

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

axisdf = don %>%
  group_by(CHR) %>%
  summarize(center=(max(BPcum) + min(BPcum)) / 2)

# All p-value manhattan plot
jpeg(
  "all_scaffold_fst.jpg",
  width = 25,
  height = 10,
  units = "in",
  res = 100
)
plot_manhattan_fst(don, "Fst - D. calciphila vs. D. sp PicN'ga") + 
  ylim(0, 1)  + 
  theme(text = element_text(size = 15),
        axis.text = element_text(color="black", size=14),
        axis.title = element_text(color = "black", size=15, face = "bold"))
dev.off()
















