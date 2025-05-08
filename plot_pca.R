#!/usr/bin/env R

# load libraries
library(tidyverse)
library(ggpubr)
library(pcadapt)

# functions
extract_file_class <- function(arr, ind, directory) {
  for (file_name in list.files(directory)) {
    file_path <- file.path(directory, file_name)
    conn <- file(file_path, open = "r")
    lines <- readLines(conn)
    arr[grep(paste(lines, collapse = "|"), ind)] <- file_name
    close(conn)
  }
  return(arr)
}

plot_pca <- function(pc_pair, pca, pve) {
  ggplot(pca, aes_string(x = pc_pair[1], 
                         y = pc_pair[2], 
                         col = "spp", 
                         shape = "other")) +
    geom_point(size = 3) + 
    scale_colour_manual(values = c("red", "blue")) + 
    xlab(paste0(pc_pair[1], 
                " (", 
                signif(pve$pve[grep(pc_pair[1], pve$pc_names)], 3), 
                "%)")) +
    ylab(paste0(pc_pair[2], 
                " (", 
                signif(pve$pve[grep(pc_pair[2], pve$pc_names)], 3), 
                "%)")) +
    theme_light() 
}

# extract script arguments
args <- commandArgs(trailingOnly = TRUE)

if(length(args) == 0) {
  stop("Requires to supply one text file")
}

out_file <- split(args[1], "[.]")[[1]][1]

pca <- read_table(args[1], col_names = FALSE)
eigenval <- scan(args[2])

# sort out the pca data
# remove nuisance column
pca <- pca[,-1]
# set names
names(pca)[1] <- "ind"
names(pca)[2:ncol(pca)] <- paste0("PC", 1:(ncol(pca)-1))

# sort out the individual species and pops
# spp
spp <- rep(NA, length(pca$ind))
spp_dir <- args[3]
spp <- extract_file_class(spp, pca$ind, spp_dir)
pca <- as.tibble(data.frame(pca, spp))

if (length(args) == 4) {
  other_dir <- args[4]
  other <- rep(NA, length(pca$ind))
  other <- extract_file_class(other, pca$ind, other_dir)  
  pca <- as.tibble(data.frame(pca, other))
}
print(other)
# first convert to percentage variance explained
pve <- data.frame(PC = 1:length(eigenval), pve = eigenval/sum(eigenval)*100)

# plot pve distribution
pdf(
  paste0(out_file, ".pve_distribution.pdf")
)
ggplot(pve, aes(PC, pve)) +
  geom_bar(stat = "identity") +
  ylab("Percentage variance explained") +
  theme_light()
dev.off()

# plot PCA scatter plots
pc_names <- names(pca)[grep("PC", names(pca))]
pve$pc_names <- pc_names
num_maj_pcs <- which((cumsum(pve$pve) > 50) == T)[1]

pc_pairs <- combn(pc_names[1:num_maj_pcs], 2, simplify = FALSE)
pc_plot_list <- lapply(pc_pairs, plot_pca, pca = pca, pve = pve)

pdf(
  paste0(out_file, ".pca.pdf"),
  width = 8.3,
  height = 11.7
)
ggarrange(plotlist = pc_plot_list, 
          nrow = ceiling(choose(num_maj_pcs, 2)/2), 
          ncol = 2, common.legend = TRUE, legend = "bottom")
dev.off()



