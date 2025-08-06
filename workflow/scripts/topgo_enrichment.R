### ENVIRONMENT
if (!require("topGO")) {
  install.packages("BiocManager")
  BiocManager::install("topGO")
}
library(tidyverse)
library(topGO)
library(patchwork)
load(file = "data/rdata/stand_var.Robject")
load(file = "data/rdata/gos_of_interest.Robject")

### FUNCTIONS
get_go_enrichment <- function(genes_of_interest) {
  gene_list <- factor(as.integer(gene_names %in% genes_of_interest))
  names(gene_list) <- gene_names
  go_data <- new("topGOdata", ontology = "BP", allGenes = gene_list,
                 annot = annFUN.gene2GO, gene2GO = gene_2_go)
  results_fis <- runTest(go_data, algorithm = "elim", statistic = "fisher")
  all_go <- usedGO(go_data)
  gen_table <- GenTable(go_data, pvalues = results_fis,
                        orderBy = "pvalues", topNodes = length(all_go))
  gen_table <- gen_table %>% filter(Annotated > 1) %>% filter(Annotated < 1000)
  gen_table$pvalues <- as.numeric(gen_table$pvalues)
  gen_table$padj <- p.adjust(gen_table$pvalues, method = "fdr")
  return(gen_table)
}

wrap.it <- function(x, len) { 
  sapply(x, function(y) paste(strwrap(y, len), 
                              collapse = "\n"), 
         USE.NAMES = FALSE)
}

wrap.labels <- function(x, len) {
  if (is.list(x))
  {
    lapply(x, wrap.it, len)
  } else {
    wrap.it(x, len)
  }
}

### VARIABLES
fg_genes_file <- snakemake@input$fg_genes
bg_genes_file <- snakemake@input$bg_genes
output_files <- snakemake@output

### EXECUTION
fg_genes <- read.delim(fg_genes_file, header = FALSE)
bg_genes <- read.delim(bg_genes_file, header = FALSE,
                       col.names = c("gene_id", "gene_desc", "GO", "accession",
                                     "misc", "db_origin"))

fg_genes <- fg_genes %>% filter(V6 != "-")
fg_genes_list <- NULL
fg_genes_col_names <- c("chr", "start", "end", "gene_id",
                        "gene_desc", "GO", "accession", "misc",
                        "db_origin", "spp_pair")
if (sum(unique(fg_genes$V10) %in% names(pair_list)) > 1) {
  fg_genes$V10 <- factor(fg_genes$V10, levels = names(pair_list))
  levels(fg_genes$V10) <- as.character(pair_list)
} else {
  fg_genes$V10 <- factor(fg_genes$V10, levels = unique(fg_genes$V10))
}
if (dim(fg_genes)[2] == 11) {
  names(fg_genes) <- c(fg_genes_col_names, "X")
  matched_entries <- match(unique(fg_genes$X), names(test_list))
  fg_genes$X <- factor(fg_genes$X,
                       levels = names(test_list[matched_entries]))
  levels(fg_genes$X) <- as.character(test_list[matched_entries])
  fg_genes_list <- split(fg_genes, f = fg_genes$X)
} else {
  names(fg_genes) <- fg_genes_col_names
  fg_genes_list <- list(unfiltered = fg_genes)
}

bg_genes <- bg_genes %>%
  mutate(GO = strsplit(as.character(GO), ",")) %>%
  unnest(GO) %>%
  mutate(GO = trimws(GO))
gene_2_go <- bg_genes %>% split(x = .$GO, f = .$gene_id)
gene_names <- names(gene_2_go)

if (!file.exists(output_files$table)) {
  gen_table_list <- lapply(names(fg_genes_list), function(x) {
    spp_pairs <- c(levels(fg_genes_list[[x]]$spp_pair), "all")
    out <- lapply(spp_pairs, function(y) {
      genes_of_interest <- if (y != "all") filter(fg_genes_list[[x]], spp_pair == y) else fg_genes_list[[x]]
      if (dim(genes_of_interest)[1] == 0)
        return(data.frame())
      genes_of_interest <- unique(genes_of_interest$gene_id)
      get_go_enrichment(genes_of_interest)
    })
    names(out) <- spp_pairs
    out
  })
  names(gen_table_list) <- names(fg_genes_list)

  ifelse(!dir.exists(dirname(output_files$table)),
        dir.create(dirname(output_files$table)),
        "Directory Exists")
  save(gen_table_list, file = output_files$table)
} else {
  load(file = output_files$table)
}

pval_threshold <- 0.05

gen_plots <- lapply(names(gen_table_list), function(x) {
  spp_names <- names(gen_table_list[[x]])
  out <- lapply(spp_names, function(y) {
    if (dim(gen_table_list[[x]][[y]])[1] == 0)
      return(ggplot())
    gen_subset <- gen_table_list[[x]][[y]] %>%
      mutate(interest = if_else(GO.ID %in% gos_of_interest, "GOs of interest", "General GOs")) %>%
      mutate(interest = factor(interest, levels = unique(interest))) %>%
      filter(pvalues < pval_threshold) %>%
      arrange(pvalues) %>%
      mutate(row_id = ave(seq_len(nrow(.)), interest, FUN = seq_along))
    gen_subset <- gen_subset[order(gen_subset$row_id,
                                   as.numeric(as.character(gen_subset$interest))), ] %>%
      head(16)

    #gen_subset <- gen_subset %>%
    #  mutate(Term = wrap.labels(Term, 15))

    ggplot(gen_subset, aes(x = reorder(Term, Significant),
                           y = (-log10(pvalues)),
                           fill = Significant)) +
      geom_col(position = position_dodge2(preserve = "single")) +
      facet_wrap(~ interest, scale = "free_y", ncol = 1) +
      scale_fill_gradient2(
        name = "Sign. genes",
        low = "blue",
        high = "darkblue"
      ) +
      coord_flip() +
      labs(x = NULL,
           y = "Enrichment -log10 p-value",
           title = paste("subset of", y, "pair(s)")) +
      theme_bw() +
      theme(
        panel.border = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        text = element_text(size = 30),
        axis.text = element_text(color="black", size=30),
        axis.title = element_text(color = "black", size=30, face = "bold"),
      )
  })
  names(out) <- spp_names
  out
})
names(gen_plots) <- names(gen_table_list)

for (name in names(gen_plots)) {
  out_file <- output_files[[1]]
  if (name != "unfiltered") {
    out_file <- output_files[[names(test_list[test_list == name])]]
  }
  ifelse(!dir.exists(dirname(out_file)),
         dir.create(dirname(out_file)),
         "Directory Exists")
  jpeg(
    out_file,
    width = 40,
    height = 40,
    unit = "in",
    res = 100,
    family = "ArialMT"
  )
  g <- gen_plots[[name]] %>%
    wrap_plots(ncol = 1) +
    plot_annotation(tag_levels = "A") &
    theme(plot.tag = element_text(size = 35, face = "bold"))
  print(g)
  dev.off()
}

