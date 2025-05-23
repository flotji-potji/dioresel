### ENVIRONMENT
library(tidyverse)

### VARIABLES
pi_intersect_file <- snakemake@input[[1]]
output_file <- snakemake@output[[1]]

plot_title <- snakemake@params$plot_title

### EXECUTION
pi_intersect_df <- read.delim(pi_intersect_file, header = FALSE,
                              col.names = c("test", "count", "pair"))

pi_intersect_df$pair <- factor(pi_intersect_df$pair,
                               levels = unique(pi_intersect_df$pair))

jpeg(
  output_file,
  width = 20,
  height = 15,
  units = "in",
  res = 100
)
ggplot(pi_intersect_df, aes(x=pair, y=count, color = test, group = test)) +
    # Show all points
    geom_point(aes(size = 15)) +
    geom_line(aes(linewidth = 15)) +
    # custom labels
    labs(y = "Num. overlapping windows",
         x = "Divergence",
         title = plot_title) +
    # Custom the theme:
    guides(fill=guide_legend(override.aes = list(alpha=1))) +
    theme_bw() +
    theme( 
      panel.border = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1, margin = margin(20,0,0,0)),
      text = element_text(size = 25),
      axis.text = element_text(color="black", size=25),
      axis.title = element_text(color = "black", size=25, face = "bold"),
    )
dev.off()
