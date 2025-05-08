#!/usr/bin/env R

# load libraries
library(tidyverse)
library(pcadapt)
library(ggpubr)
library(qvalue)

# helper functions
plot_simple_manhattan <- function(pc_file, K) {
  x <- pcadapt(input = pc_file, K = K)
  
  padj <- p.adjust(x$pvalues, method = "bonferroni")
  alpha <- 0.05
  outliers <- which(padj < alpha)
  
  plot(x, option = "manhattan") + 
    geom_hline(
      yintercept = abs(log10(max(x$pvalues[outliers]))), colour = "red"
    ) + 
    labs(title = paste0("Manhattan Plot ", "K=", K))
}

# extract script arguments
args <- commandArgs(trailingOnly = TRUE)

if(length(args) == 0) {
  stop("Requires to supply one text file")
}

out_file <- split(args[1], "[.]")[[1]][1]

# load population text file
file_path <- args[2]
conn <- file(file_path, open = "r")
lines <- readLines(conn)
split_list <- strsplit(lines, split = "(?<=[a-z])\\s*(?=[0-9A-Z])", perl = TRUE)
poplist.names <- c()
for (arr in split_list) {
  poplist.names <- append(poplist.names, arr[1])
}
close(conn)

# file path to bed file
file_path <- paste0(args[1], ".bed")
pc_file <- read.pcadapt(file_path, type = "bed")

# do preliminary pcadapt plots
x <- pcadapt(input = pc_file, K = 20)

# Score plot
h <- plot(x, option = "scores", i = 1, j = 2, pop = poplist.names)
# Scree plot
g <- plot(x, option = "screeplot")
d1 <- diff(x$singular.values^2)
cutoff <- which.max((d1 / d1[1]) < 0.1) - 1
g <- g + geom_vline(xintercept = cutoff)

pdf(
  paste0(out_file, "_scree_pca.pdf")
)
ggarrange(g, h, nrow = 1, ncol = 2)
dev.off()

# plot manhattan plots of different K pcadapts
manhattan_list <- lapply(1:cutoff, plot_simple_manhattan, pc_file = pc_file)

pdf(
  paste0(out_file, "_diffK_manhattan.pdf"),
  width = 8.3,
  height = 11.7
)
ggarrange(plotlist = manhattan_list)
dev.off()

# do ld thinning
res <- pcadapt(pc_file, K = 2, LD.clumping = list(size = 200, thr = 0.1))

padj <- p.adjust(res$pvalues, method = "bonferroni")
alpha <- 0.05
outliers <- which(padj < alpha)

pdf(
  paste0(out_file, "_ldthin_manhattan.pdf")
)
plot(res) + 
  geom_hline(
    yintercept = abs(log10(max(res$pvalues[outliers]))), colour = "red"
  ) + 
  labs(title = args[1])
dev.off()

# plot ld thinned pcadapt SNPs on physical location
snp_pos <- read_delim(paste0(args[1], ".bim"),
                      col_names = c("chr", "na", "num", "position", "ref", "alt"))
res$pos <- snp_pos$position
df <- as.tibble(data.frame(res$pos, res$pvalues))
df$pvaluelog10 <- abs(log10(df$res.pvalues))

pdf(
  paste0(out_file, "_ldthin_physloc.pdf")
)
ggplot(df, aes(x = res.pos, y = pvaluelog10)) + 
  geom_point() + 
  labs(x = "Position (bp)", y = "-log10(p-values)") + 
  labs(title = args[1])
dev.off()

