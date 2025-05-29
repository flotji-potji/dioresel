### ENVIRONMENT
library(tidyverse)

### VARIABLES
input_picmin_file <- snakemake@input[[1]]
output_pval_file <- snakemake@output[[1]]
output_padj_file <- snakemake@output[[2]]

### EXECUTION
load(file = input_picmin_file)

picMin_results <- picMin_results %>%
  mutate(start = as.numeric(scaffold),
         scaffold = factor(redundan, levels = unique(redundan)))

col_pal <- c("white", "#8ec641", "#897696", "#e93826", "#13a4f5", "#f89b56")

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