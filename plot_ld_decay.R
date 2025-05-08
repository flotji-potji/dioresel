#!/usr/bin/env R

# load libraries
library(tidyverse)

# extract script arguments
args <- commandArgs(trailingOnly = TRUE)

if(length(args) == 0) {
  stop("Requires to supply one text file")
}

out_file <- split(args[1], "[.]")[[1]][1]
my_bins <- args[1]

ld_bins<- read_tsv(my_bins)

pdf(
  paste0(out_file, ".ld_decay.pdf")
)
ggplot(ld_bins, aes(distance, avg_R2)) + 
  geom_line() + 
  geom_hline(yintercept = 0.1, colour = "black") + 
  xlab("Distance (bp)") + 
  ylab(expression(italic(r)^2))  + 
  labs(title = args[1]) +
  theme_light()
dev.off()
