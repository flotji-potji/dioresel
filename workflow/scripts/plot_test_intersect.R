### ENVIRONMENT
library(tidyverse)
if (!require("patchwork")) {
  install.packages("patchwork")
}
library(patchwork)
### FUNCTIONS

### VARIABLES
test_outputs_file <- snakemake@input$test_output
test_summary_file <- snakemake@input$test_summary
output_files <- snakemake@output

# ordered according to divergence
ordered_species_names <- c("pair_cal_spn", "pair_cal_ruf", "pair_cal_eru",
                           "pair_cal_umb", "pair_cal_vie")
ordered_test_name <- c("pcadapt", "mktest")
summary_statistics <- c("fst", "pair_pi", "pi", "tajimad")

### EXECUTION
test_outputs <- read.delim(test_outputs_file, header = TRUE, na.strings = "")
test_summary <- read.delim(test_summary_file, header = FALSE,
                           col.names = c("test", "count", "pair", "row", "num"))

for (colname in names(test_outputs)) {
  new_cols <- unlist(strsplit(colname, "\\."))
  test_outputs <- test_outputs %>%
    separate(colname, new_cols, "\\|")
  test_outputs[, new_cols[1]] <- as.numeric(test_outputs[, new_cols[1]])
}
test_outputs <- test_outputs %>%
  filter(pi_sp != "pi_sp")

test_summary <- test_summary %>%
  complete(test, pair, row, fill = list(count = 0, num = 0))
test_summary$pair <- factor(test_summary$pair,
                            levels = ordered_species_names)
test_summary$row <- factor(test_summary$row,
                           levels = ordered_test_name)
for (i in 1:length(output_files)) {
  jpeg(
    output_files[[i]],
    width = 20,
    height = 30,
    units = "in",
    res = 100
  )
  stat_name <- summary_statistics[i]
  line_stat_name <- if (stat_name == "pi") "pi_bottom" else stat_name
  if (grepl("pair", stat_name)) {
    distribution_stat_name <- paste0(
      unlist(strsplit(stat_name, "_"))[2], "_pair_sp"
    )
    distribution_pair_name <- paste0(
      unlist(strsplit(stat_name, "_"))[2], "_pair"
    )
  } else {
    distribution_stat_name <- paste0(stat_name, "_sp")
    distribution_pair_name <- stat_name
  }
  g <- test_summary %>%
    filter(test == line_stat_name) %>%
    arrange(match(pair, ordered_species_names)) %>%
    mutate(num_windows = paste0(pair, "\n", "n=", num)) %>%
    mutate(num_windows = factor(num_windows, levels = num_windows)) %>%
    ggplot(aes(x = num_windows, y = count, group = test)) +
    facet_wrap(~ row, scales = "free", ncol = 1) +
    geom_point(aes(size = 15)) +
    geom_line(aes(linewidth = 15)) +
    # custom labels
    labs(y = paste0("Num. ,", line_stat_name, " windows"),
        x = "") +
    theme_bw() +
    theme( 
      panel.border = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1, margin = margin(20,0,0,0)),
      text = element_text(size = 25),
      axis.text = element_text(color="black", size=25),
      axis.title = element_text(color = "black", size=25, face = "bold"),
    )
  set.seed(42)
  test_outputs_subset <- test_outputs[sample(nrow(test_outputs), 10000), ]
  h <- test_outputs_subset %>%
    filter(!is.na(distribution_stat_name)) %>%
    ggplot(aes(x = .data[[distribution_stat_name]],
              y = .data[[distribution_pair_name]],
              fill = .data[[distribution_stat_name]])) +
    geom_violin() +
    labs(y = stat_name,
        x = "Divergence") +
    theme_bw() +
    theme( 
      panel.border = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1, margin = margin(20,0,0,0)),
      text = element_text(size = 25),
      axis.text = element_text(color="black", size=25),
      axis.title = element_text(color = "black", size=25, face = "bold"),
    )
  print(g / h)
  dev.off()
}
