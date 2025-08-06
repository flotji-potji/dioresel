### ENVIRONMENT
if (!require("ComplexHeatmap")) {
  library(devtools)
  install_github("jokergoo/ComplexHeatmap")
}
if (!require("latex2exp"))
  install.packages("latex2exp")
library(ComplexHeatmap)
library(circlize)
library(tidyverse)
library(latex2exp)
library(svglite)
load(file = "data/rdata/stand_var.Robject")

### FUNCTIONS

### VARIABLES
input_files <- snakemake@input
output_heatmap_file <- snakemake@output[[1]]

test_names <- c("ngenes", "pcadapt", "mktest", "fst", "pi", "pi_bottom",
                "pair_pi", "pair_pi_bottom", "tajimad", "pair_tajimad")

col_fun <- colorRamp2(c(0, 25, 75, 100), c("white", "blue", "blue3",
                                           "darkblue"))

### EXECUTION
pair_matrices <- sapply(unlist(input_files), function(file) {
  intersect_df <- read.delim(file, header = FALSE)
  if (dim(intersect_df)[1] == 9) {
    intersect_df <- intersect_df %>%
      add_row(!!! setNames(as.list(rep(NA, 11)), names(.)), .after = 2)
  }
  intersect_df <- intersect_df[, -1]
  rownames(intersect_df) <- test_names
  colnames(intersect_df) <- test_names
  data.matrix(intersect_df)
}, simplify = FALSE)

# add additional matrix for mean window intersections
pair_matrices$all <- Reduce("+", lapply(pair_matrices, function(x) {
  x[is.na(x)] <- 0
  x
}))
pair_matrices$all <- pair_matrices$all / 5

for (heatmap_file in names(pair_matrices)) {
  pair_mat <- pair_matrices[[heatmap_file]]
  main_mat <- pair_mat[2:10, 2:10]
  gene_vec <- pair_mat[1, 2:10]
  
  test_max_vals <- apply(main_mat, 2, max, na.rm = TRUE)
  main_mat_prop <- (main_mat / test_max_vals) * 100
  gene_vec_prop <- (gene_vec / test_max_vals) * 100
  
  column_ha <- HeatmapAnnotation(`Proportion of\nwindows\nin genes` =
                                   anno_barplot(
                                    gene_vec_prop,
                                    labels_gp = gpar(fontsize = 8),
                                  ),
                                annotation_name_gp = gpar(fontsize = 8),
                                 height = unit(10, "mm"))
  ht <- Heatmap(main_mat_prop, bottom_annotation = column_ha,
                cluster_rows = FALSE, cluster_columns = FALSE,
                row_names_gp = gpar(fontsize = 8), row_names_side = "left",
                column_names_side = "top", row_title = rep("", 9),
                column_order = sort(colnames(main_mat)),
                row_labels = test_list_latex,
                column_names_gp = gpar(fontsize = 8), column_names_rot = 45,
                row_split = factor(rownames(main_mat),
                                   levels = sort(rownames(main_mat))),
                row_gap = unit(1, "mm"), column_labels = test_list_latex,
                layer_fun = function(j, i, x, y, width, height, fill) {
                  v <- pindex(main_mat, i, j)
                  grid.text(sprintf("%.1f", v), x, y, gp = gpar(fontsize = 6))
                  l <- fill <= "#0000FF"
                  if (sum(l) != 0)
                    grid.text(sprintf("%.1f", v[l]), x[l], y[l],
                              gp = gpar(fontsize = 6, col = "white"))
                  grid.rect(gp = gpar(lwd = 1, fill = "transparent"))
                }, col = col_fun,
                heatmap_legend_param = list(
                  title = "% overlap", legend_width = unit(5, "mm"),
                  legend_height = unit(15, "mm"),
                  labels_gp = gpar(fontsize = 6),
                  title_gp = gpar(fontsize = 7)
                ))

  output_file_name <- NULL
  if (heatmap_file != "all") {
    output_file_name <- file.path(dirname(output_heatmap_file),
                                  "pairs",
                                  paste0(str_split(heatmap_file, "/")[[1]][3],
                                         ".svg"))
    ifelse(!dir.exists(dirname(output_file_name)),
            dir.create(dirname(output_file_name)),
            "Directory Exists")
  } else {
    output_file_name <- output_heatmap_file
  }

  svglite(
    output_file_name,
    height = 4,
    width = 5,
    system_font = list(sans = "Nimbus Sans")
  )
  draw(ht, newpage = FALSE)
  dev.off()
}
