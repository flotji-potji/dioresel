### ENVIRONMENT
if (!require(SuperExactTest))
  install.packages("SuperExactTest")
library(tidyverse)
library(SuperExactTest)
library(ComplexHeatmap)
load(file = "data/rdata/stand_var.Robject")
load(file = "data/rdata/comb2group.Robject")

### FUNCTIONS
get_pairwise_subsets <- function(binary_str) {
  bits <- strsplit(binary_str, "")[[1]]
  ones <- which(bits == "1")
  if (length(ones) < 2) return(NULL)
  combs <- combn(ones, 2)
  apply(combs, 2, function(pair) {
    temp <- rep("0", length(bits))
    temp[pair] <- "1"
    paste(temp, collapse = "")
  })
}

generate_combinations <- function(n_pairs = 5, min_bits = 3) {
  all <- lapply(min_bits:n_pairs, function(k) {
    combn(n_pairs, k, function(idxs) {
      temp <- rep("0", n_pairs)
      temp[idxs] <- "1"
      paste(temp, collapse = "")
    })
  })
  unlist(all)
}

gen_comb2group <- function() {
  multi_pair_combinations <- generate_combinations()
  dxy_table <- read.delim(file = "results/pixy/dxy_means.tsv", header = FALSE,
                          col.names = c("pair", "dxy"))
  known_values <- dxy_table$dxy
  names(known_values) <- c("01100", "01010", "01001", "00110", "00101",
                           "11000", "10100", "10010", "10001", "00011")
  multi_pair_values <- sapply(multi_pair_combinations, function(bin_str) {
    subsets <- get_pairwise_subsets(bin_str)
    valid_subsets <- subsets[subsets %in% names(known_values)]
    mean(known_values[valid_subsets])
  })
  comb2val <- append(
    list("10000" = NA, "01000" = NA, "00100" = NA, "00010" = NA, "00001" = NA),
    lapply(split(known_values, names(known_values)), unname)
  )
  comb2val <- append(comb2val, multi_pair_values)
  comb2group <- cut(as.numeric(comb2val), breaks = 4,
                    labels = c("1", "2", "3", "4"))
  names(comb2group) <- names(comb2val)
  return(comb2group)
}

get_supertest_table <- function(input_mat, n = 27000) {
  Result <- supertest(input_mat, n)
  tab <- summary(Result)$Table
  tab$comb <- rownames(tab)
  tab$group <- comb2group[rownames(tab)]
  return(tab)
}

get_signif <- function(p) {
  if (is.na(p)) return("")
  if (p < 0.05) return("*")
  else return("")
}

### VARIABLES
fg_genes_file <- snakemake@input$fg_genes
gen_table_file <- snakemake@input$gen_table
output_files <- snakemake@output
gen_table_col <- snakemake@params$col

type2label <- list(gene = "Gene", func = "Function")
group_cols <- c(`1` = "steelblue1", `2` = "steelblue2",
                `3` = "steelblue3", `4` = "steelblue4")

### EXECUTION
fg_genes <- read.delim(fg_genes_file, header = FALSE)
load(file = gen_table_file)

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

fg_split <- split(fg_genes_list[[gen_table_col]],
                  f = fg_genes_list[[gen_table_col]]$spp_pair)

input_mats <- list(
  gene = sapply(as.character(pair_list), function(x) {
    fg_split[[x]] %>%
      pull(gene_id) %>%
      unique(.)
  }),
  func = sapply(as.character(pair_list), function(x) {
    if (x != "all") {
      gen_table_list[[gen_table_col]][[x]] %>%
        filter(pvalues < 0.05) %>%
        pull(GO.ID)
    }
  })
)

supertests <- lapply(names(input_mats), function(x) {
  n <- 27000
  if (x == "func") n <- dim(gen_table_list[[gen_table_col]][[1]])[1]
  out <- get_supertest_table(input_mats[[x]], n)
  out
})
names(supertests) <- names(input_mats)

for (name in names(supertests)) {

  mat <- input_mats[[name]]
  tab <- supertests[[name]]
  type_label <- type2label[[name]]

  comb_mat <- make_comb_mat(mat, mode = "intersect")
  comb_codes <- comb_name(comb_mat)
  comb_groups <- comb2group[comb_codes]
  comb_cols   <- group_cols[comb_groups]
  filtered_degree <- comb_degree(comb_mat) > 1
  filtered_mat <- comb_mat[filtered_degree]
  intersection_sizes <- comb_size(filtered_mat)

  pvals <- tab$P.value
  names(pvals) <- tab$comb
  pvals <- pvals[names(sort(-intersection_sizes))]
  signif_labels <- vapply(pvals, get_signif, FUN.VALUE = "")

  ht <- UpSet(
    filtered_mat,
    comb_order = order(-intersection_sizes),
    set_order = names(mat), column_names_gp = gpar(fontsize = 10),
    pt_size = unit(3, "mm"),
    lwd = 2,
    column_title = paste(type_label, "intersections"),
    row_title = "",
    top_annotation = HeatmapAnnotation(
      "Intersection\nsize" = anno_barplot(
        intersection_sizes,
        border = FALSE,
        gp = gpar(fill = comb_cols[filtered_degree]),
        height = unit(3.5, "cm"),
        axis_param = list(side = "left")
      ),
      annotation_name_rot = 90,
      annotation_name_side = "left"
    ),
    heatmap_width = unit(10, "cm"),
    heatmap_height = unit(6, "cm")
  )

  pdf(
    output_files[[name]],
    family = "ArialMT",
  )
  draw(ht)
  decorate_annotation("Intersection\nsize", {
    for (i in seq_along(signif_labels)) {
      size <- sort(intersection_sizes, decreasing = TRUE)[i]
      grid.text(
        signif_labels[i],
        x = unit(i, "native"),
        y = unit(size, "native") + unit(1, "mm")
      )
    }
  })
  dev.off()

}

