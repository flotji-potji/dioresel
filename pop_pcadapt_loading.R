#!/usr/bin/env R

# load modules
library(pcadapt)
library(ggplot2)
library(ggpubr)
library(gtools)

# helper functions
plot_loading <- function(file) {
  file_path <- paste0(file)
  pc_file <- read.pcadapt(file_path, type = "bed")
  # run pcadapt for K=1
  x <- pcadapt(input = pc_file, K = 1)
  
  df <- data.frame(x$loadings[, 1], seq(1, length(x$loadings)))
  colnames(df) <- c("Loadings", "Index")
  
  ggplot(df, aes(y=Loadings, x = Index)) + 
    geom_point() + 
    labs(title = unlist(strsplit(file, "/"))[1],
         y = paste0("Loadings PC", 1)) + 
    theme_light()
}

plot_loading_fixed <- function(file, ylim=NULL) {
  plot_loading(file) +
    ylim(ylim[1], ylim[2])
}

# main part
dir <- "/lisc/scratch/botany/frschmidt/pcadapt_pca_ldDecay_cs/"
setwd(dir)
pattern <- "*.vcf.gz_all_cal_spn_ind.bed"
files <- list.files(dir, pattern, recursive = T)
files <- mixedsort(files)

jpeg(
  "pcadapt_loadings.jpeg",
  width = 20,
  height = 23,
  units = "in",
  res = 100
)
plots <- lapply(files, plot_loading)
ggarrange(plotlist = plots)
dev.off()

jpeg(
  "pcadapt_loadings_fixed.jpeg",
  width = 20,
  height = 23,
  units = "in",
  res = 100
)
plots <- lapply(files, plot_loading_fixed, ylim=c(-0.02, 0.02))
ggarrange(plotlist = plots)
dev.off()
