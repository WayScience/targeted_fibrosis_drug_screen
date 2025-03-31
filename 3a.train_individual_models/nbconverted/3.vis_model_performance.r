suppressPackageStartupMessages({
    library(ggplot2)
    library(tidyr)
    library(dplyr)
    library(arrow)
    library(ggridges)
    library(RColorBrewer)
})

output_dir <- "figures"
if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
}

performance_metrics <- read_parquet("./performance_metrics/performance_metrics.parquet")

head(performance_metrics)

width = 12
height = 8
options(repr.plot.width = 12, repr.plot.height = 8)  # Adjust width and height as desired

# Filter the data for Train and Test only
filtered_data <- performance_metrics[performance_metrics$dataset %in% c("train", "test"), ]

# Create the ggplot for PR curves
test_train_all_models <- ggplot(filtered_data, aes(x = recall, y = precision, color = dataset, linetype = model_type)) +
    geom_line() +
    facet_wrap(plate_trained ~ .) +
    labs(
        x = "Recall",
        y = "Precision",
        color = "Data Type",
        linetype = "Model Type"
    ) +
    theme_bw() +
    theme(
        text = element_text(size = 16),  # Increase font size for all text
        axis.title = element_text(size = 18),  # Increase font size for axis titles
        axis.text = element_text(size = 14),  # Increase font size for axis text
        legend.title = element_text(size = 16),  # Increase font size for legend title
        legend.text = element_text(size = 14),  # Increase font size for legend text
        strip.text = element_text(size = 16)  # Increase font size for facet labels
    )

# Save the plot to the output directory
ggsave(file.path(output_dir, "test_train_all_models_pr_curves.png"), test_train_all_models, dpi = 500, height = height, width = width)

test_train_all_models

width = 12
height = 8
options(repr.plot.width = 12, repr.plot.height = 8)  # Adjust width and height as desired

# Filter out everything that is not Train and Test, and exclude those with 'test_' as a prefix
filtered_data <- performance_metrics[!(performance_metrics$dataset %in% c("train", "test") | grepl("^test_", performance_metrics$dataset)), ]

# Create the ggplot for PR curves with larger font sizes
holdout_pr_curves <- ggplot(filtered_data, aes(x = recall, y = precision, color = dataset, linetype = model_type)) +
    geom_line() +
    facet_wrap(plate_trained ~ .) +
    labs(
        x = "Recall",
        y = "Precision",
        color = "Data split",
        linetype = "Model type"
    ) +
    theme_bw() +
    theme(
        text = element_text(size = 16),  # Increase font size for all text
        axis.title = element_text(size = 18),  # Increase font size for axis titles
        axis.text = element_text(size = 14),  # Increase font size for axis text
        legend.title = element_text(size = 16),  # Increase font size for legend title
        legend.text = element_text(size = 14),  # Increase font size for legend text
        strip.text = element_text(size = 16)  # Increase font size for facet labels
    )

# Save the plot to the output directory
ggsave(file.path(output_dir, "holdout_plate_pr_curves.png"), holdout_pr_curves, dpi = 500, height = height, width = width)

holdout_pr_curves

width = 12
height = 8
options(repr.plot.width = 12, repr.plot.height = 8)  # Adjust width and height as desired

# Filter to only keep rows with 'test_' as a prefix in the dataset column
filtered_data <- performance_metrics[grepl("^test_", performance_metrics$dataset), ]

# Create the ggplot for PR curves with larger font sizes
holdout_combined_model_curves <- ggplot(filtered_data, aes(x = recall, y = precision, color = dataset, linetype = model_type)) +
    geom_line() +
    labs(
        x = "Recall",
        y = "Precision",
        color = "Data split",
        linetype = "Model type"
    ) +
    theme_bw() +
    theme(
        text = element_text(size = 16),  # Increase font size for all text
        axis.title = element_text(size = 18),  # Increase font size for axis titles
        axis.text = element_text(size = 14),  # Increase font size for axis text
        legend.title = element_text(size = 16),  # Increase font size for legend title
        legend.text = element_text(size = 14),  # Increase font size for legend text
        strip.text = element_text(size = 16)  # Increase font size for facet labels
    )

# Save the plot to the output directory
ggsave(file.path(output_dir, "holdout_combined_model_pr_curves.png"), holdout_combined_model_curves, dpi = 500, height = height, width = width)

holdout_combined_model_curves

# Filter to only include DMSO treatment
DMSO_performance_df <- performance_metrics %>% filter(Metadata_treatment == "DMSO")

# Filter the data for Train and Test only
DMSO_performance_df <- DMSO_performance_df[DMSO_performance_df$dataset %in% c("test"), ]

height <- 8
width <- 16
options(repr.plot.width = width, repr.plot.height = height)

ridge_plot_controls <- ggplot(DMSO_performance_df, aes(x = predicted_probability, y = actual_label, fill = Metadata_treatment)) +
  geom_density_ridges(aes(fill = Metadata_treatment), alpha = 0.7, scale = 2, rel_min_height = 0.01, bandwidth = 0.1) +
  scale_fill_manual(values = c("DMSO" = brewer.pal(8, "Dark2")[7])) +  # Only include DMSO color
  geom_vline(xintercept = 1, linetype = "dashed", color = "black") +
  scale_x_continuous(breaks = seq(0, 1, 0.5)) +
  facet_grid(model_type ~ plate_trained, scales = "free_y") + 
  labs(x = "Probability of healthy prediction", y = "Heart Type") +  # Update x-axis label
  theme_bw() +
  theme(
    legend.position = "none",
    axis.text = element_text(size = 20),
    axis.text.x = element_text(size = 20),
    axis.title = element_text(size = 24),
    strip.text = element_text(size = 18),
    strip.background = element_rect(
      colour = "black",
      fill = "#fdfff4"
    )
  )

# Save the plot to the output directory
ggsave(file.path(output_dir, "DMSO_controls_prob_ridge_plot_testing_data.png"), ridge_plot_controls, dpi = 500, height = height, width = width)

ridge_plot_controls
