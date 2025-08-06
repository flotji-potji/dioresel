### ENVIRONMENT
if (!require("ComplexHeatmap")) {
  library(devtools)
  install_github("jokergoo/ComplexHeatmap")
}
library(ComplexHeatmap)
library(circlize)
library(tidyverse)
library(svglite)
load(file = "data/rdata/stand_var.Robject")

### FUNCTIONS

### VARIABLES
input_count_file <- snakemake@input[[1]]
output_heatmap_file <- snakemake@output[[1]]

number_of_genes <- 39

### EXECUTION
# data table modifications
goi_count <- read.delim(input_count_file, header = FALSE,
                        col.names = c("test", "count", "pair"))
goi_count <- goi_count %>%
  pivot_wider(names_from = test, values_from = count, values_fill = 0)
goi_count <- goi_count %>% arrange(match(pair, names(pair_list)))
goi_count$mktest <- 0
goi_count <- goi_count %>% column_to_rownames(var = "pair")
goi_count$mktest[1:2] <- NA

# matrix modifications
mat <- data.matrix(goi_count)
mat <- (mat / number_of_genes) * 100
mat <- t(mat)
colnames(mat) <- as.character(pair_list)
mat <- mat[match(names(test_list), rownames(mat)), ]
rownames(mat) <- unlist(test_list)

# heatmap color configuration
col_fun <- colorRamp2(c(0, 25, 100), c("white", "blue", "darkblue"))

# plotting and heatmap modifications
svglite(
  output_heatmap_file,
  width = 2.7,
  height = 3.5,
  system_font = list(sans = "Nimbus Sans")
)
ht <- Heatmap(mat, cluster_rows = FALSE, cluster_columns = FALSE,
              column_gap = unit(1, "mm"), column_names_rot = 45,
              row_names_gp = gpar(fontsize = 8, font = "Helvetica"),
              row_names_side = "left", column_names_side = "top",
              row_labels = test_list_latex, 
              column_names_gp = gpar(fontsize = 8, font = "Helvetica"),
              column_split = factor(colnames(mat), levels = colnames(mat)),
              row_order = sort(rownames(mat)), column_title = rep("", 5),
              layer_fun = function(j, i, x, y, width, height, fill) {
                v <- pindex(mat, i, j)
                grid.text(sprintf("%.1f", v), x, y, gp = gpar(fontsize = 8))
                grid.rect(gp = gpar(lwd = 1, fill = "transparent",
                                    font = "Helvetica"))
              }, col = col_fun,
              heatmap_legend_param = list(
                title = "% overlap", legend_width = unit(1, "inch"),
                legend_height = unit(1.5, "inch"), title_position = "lefttop",
                direction = "horizontal", position = "bottom",
                labels_gp = gpar(fontsize = 6, font = "Helvetica"),
                title_gp = gpar(fontsize = 7, font = "Helvetica")
              ))
draw(ht, heatmap_legend_side = "bottom")
dev.off()
