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

pr_results <- read_parquet("./performance_metrics/batch1_pr_curve_results.parquet")

head(pr_results)

# Create a mapping for localhost names
localhost_names <- setdiff(unique(pr_results$plate_trained), "combined_batch1")
plate_mapping <- setNames(paste0("plate", seq_along(localhost_names)), localhost_names)
plate_mapping["combined_batch1"] <- "combined_batch1" # keep the original name for combined_batch1

# Apply the mapping to the dataframe
pr_results$plate_trained <- plate_mapping[pr_results$plate_trained]

# Convert named vector to a data frame
plate_mapping_df <- data.frame(
  original = names(plate_mapping),
  mapped = unname(plate_mapping)
)

# Save as CSV
write.csv(plate_mapping_df, "batch1_plate_name_mapping.csv", row.names = FALSE)

# Display the mapping
print(plate_mapping)
head(pr_results)

# Copy original dataset values
dataset_original <- pr_results$dataset

# Extract base name from 'test_localhostX' or 'holdout_localhostY'
base_names <- sub("^(test_|holdout_)", "", dataset_original)

# Map base names if they exist in plate_mapping, otherwise keep original
mapped_base <- ifelse(base_names %in% names(plate_mapping),
                      plate_mapping[base_names],
                      base_names)

# Reattach prefixes (test_/holdout_/none)
pr_results$dataset <- ifelse(
  grepl("^test_", dataset_original),
  paste0("test_", mapped_base),
  ifelse(
    grepl("^holdout_", dataset_original),
    paste0("holdout_", mapped_base),
    mapped_base
  )
)

# Display the updated dataframe and check that the mapping worked
unique(pr_results$dataset)
head(pr_results)

width = 12
height = 8
options(repr.plot.width = 12, repr.plot.height = 8)  # Adjust width and height as desired

# Filter the data for Train and Test only
filtered_data <- pr_results[pr_results$dataset %in% c("train", "test"), ]

# Create the ggplot for PR curves
test_train_all_models <- ggplot(filtered_data, aes(x = recall, y = precision, color = dataset, linetype = model_type)) +
    geom_line(linewidth = 1.15) +
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

width = 10
height = 8
options(repr.plot.width = 12, repr.plot.height = 8)  # Adjust width and height as desired

# Filter out everything that is not Train and Test, and exclude those with 'test_' as a prefix
filtered_data <- pr_results[!(pr_results$dataset %in% c("train", "test") | grepl("^test_", pr_results$dataset)), ]

# Create the ggplot for PR curves with larger font sizes
holdout_pr_curves <- ggplot(filtered_data, aes(x = recall, y = precision, color = dataset, linetype = model_type)) +
    geom_line(linewidth = 1.15) +
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

width = 8
height = 6
options(repr.plot.width = 12, repr.plot.height = 8)  # Adjust width and height as desired

# Filter to only keep rows with 'test_' as a prefix in the dataset column
filtered_data <- pr_results[grepl("^test_", pr_results$dataset), ]

# Create the ggplot for PR curves with larger font sizes
holdout_combined_model_curves <- ggplot(filtered_data, aes(x = recall, y = precision, color = dataset, linetype = model_type)) +
    geom_line(linewidth = 1.15) +
    labs(
        x = "Recall",
        y = "Precision",
        color = "Data split",
        linetype = "Model type"
    ) +
    theme_bw() +
    ylim(0,1) +
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

probability_results <- read_parquet("./performance_metrics/batch1_probabilities_DMSO_results.parquet")

# Apply the mapping to the dataframe
probability_results$plate_trained <- plate_mapping[probability_results$plate_trained]

head(probability_results)

# # Filter to only include DMSO treatment
# DMSO_performance_df <- performance_metrics %>% filter(Metadata_treatment == "DMSO")

# Filter the data for test data only
DMSO_test_performance_df <- probability_results[probability_results$dataset %in% c("test"), ]

height <- 8
width <- 16
options(repr.plot.width = width, repr.plot.height = height)

ridge_plot_controls <- ggplot(DMSO_test_performance_df, aes(x = predicted_probability, y = actual_label, fill = Metadata_treatment)) +
  geom_density_ridges(aes(fill = Metadata_treatment), alpha = 0.7, scale = 2, rel_min_height = 0.01, bandwidth = 0.1) +
  scale_fill_manual(values = c("DMSO" = brewer.pal(8, "Dark2")[7])) +  # Only include DMSO color
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
