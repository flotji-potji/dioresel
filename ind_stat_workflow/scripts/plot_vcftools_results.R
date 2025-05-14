### ENVIRONMENT
library(tidyverse)

### FUNCTIONS
df_res_table <- function(file, col) {
  file_df <- read.delim(file)
  file_df <- file_df[c("CHROM", "BIN_START", col)]
  return(file_df)
}

manhattan_data_frame <- function(df) {
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

  return(don)
}

plot_manhattan_fst <- function(df, col_type, title="") {
  axisdf = df %>%
    group_by(CHR) %>%
    summarize(center=(max(BPcum) + min(BPcum)) / 2)

  ggplot(df, aes(x=BPcum, y=pval)) +
    # Show all points
    geom_point(aes(color=as.factor(CHR)), alpha=0.8, size=1.3) +
    scale_color_manual(values = rep(c("black", "grey"), 31)) +
    
    # custom X axis:
    scale_x_continuous(label = axisdf$CHR, breaks= axisdf$center ) +
    scale_y_continuous(expand = c(0, 0) ) +     # remove space between plot area and x axis
    
    # custom labels
    labs(y = col_type,
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
input_res_table <- snakemake@input[[1]]
output_jpg <- snakemake@output[[1]]
col_type <- snakemake@params$col_type

### EXECUTION

manhattan_df <- manhattan_data_frame(df_res_table(input_res_table, col_type))

jpeg(
  output_jpg,
  width = 25,
  height = 10,
  units = "in",
  res = 100
)
plot_manhattan_fst(manhattan_df, col_type,
                   paste(snakemake@params$plot_title,
                         snakemake@params$window_info)) +
  theme(text = element_text(size = 15),
        axis.text = element_text(color = "black", size = 14),
        axis.title = element_text(color = "black", size = 15, face = "bold"))
dev.off()
