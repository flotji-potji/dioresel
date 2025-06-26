### ENVIRONMENT
library(tidyverse)
if (!require(ggVennDiagram)) {
  devtools::install_github("gaospecial/ggVennDiagram")
}
library(ggVennDiagram)

### VARIABLES
input_picmin_file <- snakemake@input[[1]]
output_pval_file <- snakemake@output[[1]]
output_padj_file <- snakemake@output[[2]]
output_lin_venn <- snakemake@output[[3]]
output_est_venn <- snakemake@output[[4]]
output_lin_est <- snakemake@output[[5]]


### EXECUTION
load(file = input_picmin_file)

picMin_results <- picMin_results %>%
  mutate(start = as.numeric(scaffold),
         scaffold = factor(redundan, levels = unique(redundan)))

# numLin venn diagram
jpeg(
  output_lin_venn,
  width = 10,
  height = 10,
  units = "in",
  res = 100
)
picmin_venn <- picMin_results %>%
  pivot_wider(values_from = locus, names_from = numLin) %>%
  select(tail(names(.), 3))
ggVennDiagram(picmin_venn)
dev.off()

# n_est venn diagram
jpeg(
  output_est_venn,
  width = 10,
  height = 10,
  units = "in",
  res = 100
)
picmin_venn <- picMin_results %>%
  pivot_wider(values_from = locus, names_from = n_est) %>%
  select(tail(names(.), 4))
ggVennDiagram(picmin_venn)
dev.off()

# numLin x n_est heatmap
jpeg(
  output_lin_est,
  width = 15,
  height = 15,
  units = "in",
  res = 100
)
freq_table <- picMin_results %>%
  count(numLin, n_est)
ggplot(freq_table, aes(x = factor(numLin), y = factor(n_est), fill = n)) +
  geom_tile(color = "white") +
  geom_text(aes(label = n), color = "white", size = 5) +
  scale_fill_gradient(low = "lightblue", high = "darkblue") +
  labs(title = "Frequencies of (Col1, Col2) Combinations",
       x = "numLin", y = "n_est", fill = "Count") +
  theme_minimal()
dev.off()

col_pal <- c("white", "#8ec641", "#897696", "#e93826", "#13a4f5", "#f89b56")

# pval manhattan plot
jpeg(
  output_pval_file,
  width = 20,
  height = 25,
  units = "in",
  res = 100
)
ggplot(data = picMin_results,
       aes(x = start/1e6,
           y = -log10(p),
           fill = factor(n_est))) +
  facet_wrap(~ scaffold, scales = "free_x", ncol = 4) +
  geom_point(shape = 21,
             size = 4) +
  geom_hline(aes(yintercept = -log10(0.05)),
             lty=2)+
  scale_fill_manual("Number of\nLineages\n",values = col_pal)+
  scale_y_continuous(expression(-log[10]*"(q-value)"))+
  scale_x_continuous("Position in Scaffold (Mbp)")+
  theme_classic() + 
  theme(strip.background = element_blank(),
        strip.text = element_text(face = "bold"))
dev.off()

# padj manhattan plot
jpeg(
  output_padj_file,
  width = 20,
  height = 25,
  units = "in",
  res = 100
)
ggplot(data = picMin_results,
       aes(x = start/1e6,
           y = -log10(pooled_q),
           fill = factor(n_est))) +
  facet_wrap(~ scaffold, scales = "free_x") +
  geom_point(shape = 21,
             size = 4) +
  geom_hline(aes(yintercept = -log10(0.05)),
             lty=2)+
  scale_fill_manual("Number of\nLineages\n",values = col_pal)+
  scale_y_continuous(expression(-log[10]*"(q-value)"))+
  scale_x_continuous("Position in Scaffold (Mbp)")+
  theme_classic() + 
  theme(strip.background = element_blank(),
        strip.text = element_text(face = "bold"))
dev.off()