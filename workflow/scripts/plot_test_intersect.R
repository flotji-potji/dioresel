### ENVIRONMENT
library(tidyverse)
if (!require("patchwork")) {
  install.packages("patchwork")
}
library(patchwork)
library(svglite)
load(file = "data/rdata/stand_var.Robject")

### FUNCTIONS

### VARIABLES
test_outputs_file <- snakemake@input$test_output
test_summary_file <- snakemake@input$test_summary
output_files <- snakemake@output

# ordered according to divergence
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
                            levels = names(pair_list))
levels(test_summary$pair) <- as.character(pair_list)
test_summary$row <- factor(test_summary$row,
                           levels = ordered_test_name)
levels(test_summary$row) <- as.character(test_list[ordered_test_name])
test_summary$test_prop <- (test_summary$count / test_summary$num) * 100

for (i in 1:length(output_files)) {
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
    arrange(match(pair, names(pair_list))) %>%
    mutate(label = case_when(is.na(test_prop) ~ "NA", test_prop > 0 ~ "",
                             TRUE ~ "0")) %>%
    mutate(count = as.numeric(count)) %>%
    ggplot(aes(x = pair, y = test_prop, group = test, fill = count)) +
    facet_wrap(~ row, scales = "free_y", ncol = 1) +
    geom_bar(stat = "identity") +
    # custom labels
    labs(y = paste("% Prop.", test_list[[line_stat_name]], "windows"),
         x = "") +
    geom_text(
      aes(y = 0, label = label, group = test),
      position = position_dodge(.9),
      vjust = 0, size = 2
    ) +
    scale_fill_gradient2(
      name = "Window\noverlap",
      midpoint = 500,
      low = "darkgreen",
      mid = "yellow",
      high = "orange",
    ) +
    theme_minimal() +
    theme(
      panel.border = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      axis.title.x = element_blank(),
      axis.text.x = element_blank(),
      text = element_text(size = 8),
      axis.text = element_text(color="black", size=6),
      axis.title = element_text(color = "black", size=8, face = "bold"),
      strip.text = element_text(size = 7, face = "bold"),
      legend.key.height = unit(3, "mm"),
      legend.key.width  = unit(1.5, "mm"),
      legend.title = element_text(size = 6),
      legend.text = element_text(size = 5),
    )

  to_subset <- test_outputs[, c(distribution_stat_name,
                                distribution_pair_name)]
  to_subset <- to_subset %>%
    filter(.data[[distribution_pair_name]] != "NA")
  to_subset$stat_order <- factor(to_subset[[distribution_stat_name]],
                                 levels = names(pair_list))
  levels(to_subset$stat_order) <- as.character(pair_list)
  n_group <- to_subset %>%
    group_by(.data[[distribution_stat_name]]) %>%
    summarize(n_test = n())
  to_subset2 <- merge(to_subset, n_group)

  #set.seed(42)
  #test_outputs_subset <- test_outputs[sample(nrow(to_subset2), 10000), ]
  h <- to_subset2 %>%
    filter(!is.na(distribution_stat_name)) %>%
    filter(.data[[distribution_pair_name]] != "NA") %>%
    mutate(n_test = as.numeric(n_test)) %>%
    ggplot(aes(x = stat_order,
               y = .data[[distribution_pair_name]])) +
    geom_boxplot(
      aes(fill = n_test),
      outlier.size = 0.04,
      lwd = 0.4
    ) +
    scale_fill_gradient2(
      name = "Distribution\nsize",
    ) +
    labs(y = test_list_latex[[stat_name]]) +
    theme_minimal() +
    theme(
      panel.border = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      axis.title.x = element_blank(),
      text = element_text(size = 8),
      axis.text = element_text(color="black", size=6),
      axis.title = element_text(color = "black", size=8, face = "bold"),
      legend.key.height = unit(3, "mm"),
      legend.key.width  = unit(1.5, "mm"),
      legend.title = element_text(size = 6),
      legend.text = element_text(size = 5)
    )
  svglite(
    output_files[[i]],
    width = 2.7,
    height = 3.5,
    system_font = list(sans = "Nimbus Sans")
  )
  print(g / h)
  dev.off()
}
