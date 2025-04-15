suppressPackageStartupMessages(suppressWarnings(library(ggplot2))) # plotting
suppressPackageStartupMessages(suppressWarnings(library(dplyr))) # data manipulation
suppressPackageStartupMessages(suppressWarnings(library(ggridges))) # ridge line plots
suppressPackageStartupMessages(suppressWarnings(library(RColorBrewer))) # color palettes
suppressPackageStartupMessages(suppressWarnings(library(arrow))) # parquet files
suppressPackageStartupMessages(suppressWarnings(library(ggrepel))) # add names to plots

# Specify the path for the figures directory
figures_dir <- "./figures"

# Create the figures directory if it doesn't exist
if (!dir.exists(figures_dir)) {
  dir.create(figures_dir, recursive = TRUE)
}

# Load in the probabilities
combined_probabilities_path <- file.path(
    "./performance_metrics/batch1_compound_predictions.parquet"
)

# Read in the data from the parquet file
combined_probabilities_df <- arrow::read_parquet(
    combined_probabilities_path
)

# Load the plate name mapping file
plate_name_mapping_path <- file.path(
    "./batch1_plate_name_mapping.csv"
)

plate_name_mapping <- read.csv(plate_name_mapping_path)

# Map the model_name to the new name
combined_probabilities_df <- combined_probabilities_df %>%
    left_join(plate_name_mapping, by = c("model_name" = "original")) %>%
    mutate(model_name = ifelse(!is.na(mapped), mapped, model_name)) %>%
    select(-mapped)

# Display dimensions and a preview of the data
dim(combined_probabilities_df)
head(combined_probabilities_df, 2)

# Custom color palette inspired by Dark2
custom_palette <- c(
  "Angiogenesis" = "#1b9e77",
  "Apoptosis" = "#d95f02",
  "DNA Damage" = "#7570b3",
  "Endocrinology & Hormones" = "#e7298a",
  "Epigenetics" = "#66a61e",
  "MAPK" = "#e6ab02",
  "Metabolism" = "#a6761d",
  "Neuronal Signaling" = "#667665",
  "Others" = "#b3b3b3",
  "PI3K/Akt/mTOR" = "#8dd3c7",
  "Stem Cells &  Wnt" = "#fb8072"
)

# Filter out shuffle models
filtered_df <- combined_probabilities_df %>% filter(model_type != "shuffled")

# Step 1: Calculate the median Healthy_probas for each treatment (you can use mean instead)
treatment_order <- filtered_df %>%
  group_by(Metadata_treatment) %>%
  summarise(median_healthy_proba = median(predicted_probas)) %>%
  arrange(desc(median_healthy_proba))  # Sort in descending order of median

# Step 2: Reorder the 'Metadata_treatment' factor levels in reverse order for top-to-bottom plotting
filtered_df <- filtered_df %>%
  mutate(Metadata_treatment = factor(Metadata_treatment, levels = rev(treatment_order$Metadata_treatment)))

# Step 3: Create the ridge plot with reordered treatments
height <- 15
width <- 24
options(repr.plot.width = width, repr.plot.height = height)

# Create a ridge plot excluding DMSO and using all other treatments on the y-axis
ridge_plot_non_DMSO <- ggplot(filtered_df, aes(x = predicted_probas, y = Metadata_treatment, fill = Metadata_Pathway)) +
  geom_density_ridges(alpha = 0.7, scale = 2.25, rel_min_height = 0.01, bandwidth = 0.1) +  # Ridges colored by Pathway
  geom_vline(xintercept = 1, linetype = "dashed", color = "black") +
  scale_x_continuous(breaks = seq(0, 1, 0.5)) +
  labs(x = "Probability of healthy prediction", y = "Treatment", color = "Pathway") +  # Update axis labels
  scale_fill_manual(values = custom_palette) +  # Use the custom color palette
  theme_bw() +
  theme(
    legend.position = "right",  # Position the legend on the right
    legend.title = element_text("Pathway", size = 22),
    axis.text = element_text(size = 22),
    axis.text.x = element_text(size = 22),
    axis.title = element_text(size = 24),
    strip.text = element_text(size = 18),
    strip.background = element_rect(
      colour = "black",
      fill = "#fdfff4"
    ),
    legend.text = element_text(size = 18),  # Adjust legend text size if needed
    legend.key.size = unit(1.5, "lines")  # Adjust the size of the legend keys
  ) +
  facet_wrap(~ model_name, ncol = 5)  # Facet by model applied

# Save figure
ggsave(paste0(figures_dir, "/treatments_only_ridge_plot.png"), ridge_plot_non_DMSO, height = height, width = width, dpi = 500)

ridge_plot_non_DMSO

# Get the median probability and total cell counts per compound for only the combined model
summary_data <- filtered_df %>%
    filter(model_name == "combined_batch1") %>%
    group_by(Metadata_treatment) %>%
    summarise(
        median_predicted_proba = median(predicted_probas),
        cell_count = n(),
        .groups = "drop"
    ) %>%
    left_join(
        filtered_df %>% distinct(Metadata_treatment, Metadata_Pathway),
        by = "Metadata_treatment"
    ) %>%
    mutate(
        high_count_and_proba = cell_count > 300 & median_predicted_proba > 0.65
    )

# Display the summary data
dim(summary_data)
head(summary_data, 2)

# Custom color palette inspired by Dark2 (includes one more color for a pathway missing in above plot)
custom_palette <- c(
  "Angiogenesis" = "#1b9e77",
  "Apoptosis" = "#d95f02",
  "DNA Damage" = "#7570b3",
  "Endocrinology & Hormones" = "#e7298a",
  "Epigenetics" = "#66a61e",
  "MAPK" = "#e6ab02",
  "Metabolism" = "#a6761d",
  "Neuronal Signaling" = "#667665",
  "Others" = "#b3b3b3",
  "PI3K/Akt/mTOR" = "#8dd3c7",
  "Stem Cells &  Wnt" = "#fb8072",
  "GPCR & G Protein" = "#984ea3"
)

height <- 10
width <- 14
options(repr.plot.width = width, repr.plot.height = height)

# Generate scatterplot for probability and cell count
count_probas_plot <- ggplot(summary_data, aes(x = cell_count, y = median_predicted_proba, color = Metadata_Pathway)) +
    geom_point(size = 4, alpha = 0.7) +
    geom_text_repel(data = subset(summary_data, high_count_and_proba == TRUE), aes(label = Metadata_treatment), 
                    size = 6, box.padding = 0.35, point.padding = 0.5, max.overlaps = 20, show.legend = FALSE) +  # Add labels
    scale_color_manual(values = custom_palette) +  # Use the custom color palette
    labs(
        x = "Cell count\n(across plates)",
        y = "Median predicted probability\n(combined model)",
        color = "Treatment",
    ) +
    theme_bw() +
    theme(
        legend.position = "right",
        legend.title = element_text(size = 18),
        legend.text = element_text(size = 16),
        axis.text = element_text(size = 16),
        axis.title = element_text(size = 18)
    )

# Save the scatterplot
ggsave(paste0(figures_dir, "/batch1_median_proba_vs_cell_count_compounds.png"), count_probas_plot, height = height, width = width, dpi = 500)

count_probas_plot

# Read in the mAP scores CSV file
mAP_scores_path <- "./_all_mAP_scores.csv"
mAP_scores <- read.csv(mAP_scores_path)
mAP_scores <- dplyr::filter(mAP_scores, shuffled == "False") # Filter only for non-shuffled results

# Set treatment as a character to match for merging
summary_data$Metadata_treatment <- as.character(summary_data$Metadata_treatment)

# Merge the mAP scores with the summary_data using the Metadata_treatment column
summary_data_w_mAP <- summary_data %>%
    left_join(mAP_scores %>% select(Metadata_treatment, negative_mean_average_precision, positive_mean_average_precision), 
        by = "Metadata_treatment")

# Drop rows with NaNs in the mAP scores
summary_data_w_mAP <- summary_data_w_mAP %>%
  filter(!is.na(negative_mean_average_precision) & !is.na(positive_mean_average_precision))


# Display the updated summary_data
dim(summary_data_w_mAP)
head(summary_data_w_mAP)

setdiff(summary_data$Metadata_treatment, mAP_scores$Metadata_treatment)

height <- 10
width <- 14
options(repr.plot.width = width, repr.plot.height = height)

# Generate scatterplot for probability and cell count
count_probas_plot <- ggplot(summary_data_w_mAP, aes(x = cell_count, y = median_predicted_proba, color = Metadata_Pathway)) +
    geom_point(aes(size = negative_mean_average_precision), alpha = 0.7) +
    geom_text_repel(data = subset(summary_data, high_count_and_proba == TRUE), aes(label = Metadata_treatment), 
                    size = 6, box.padding = 0.35, point.padding = 0.5, max.overlaps = 20, show.legend = FALSE) +  # Add labels
    scale_color_manual(values = custom_palette) +  # Use the custom color palette
    scale_size_continuous(name = "Negative mAP score") +
    labs(
        x = "Cell count\n(across plates)",
        y = "Median predicted probability\n(combined model)",
        color = "Treatment",
    ) +
    theme_bw() +
    theme(
        legend.position = "right",
        legend.title = element_text(size = 18),
        legend.text = element_text(size = 16),
        axis.text = element_text(size = 16),
        axis.title = element_text(size = 18)
    )

# Save the scatterplot
ggsave(paste0(figures_dir, "/batch1_neg_mAP_median_proba_vs_cell_count_compounds.png"), count_probas_plot, height = height, width = width, dpi = 500)

count_probas_plot

height <- 10
width <- 14
options(repr.plot.width = width, repr.plot.height = height)

# Generate scatterplot for probability and cell count
count_probas_plot <- ggplot(summary_data_w_mAP, aes(x = cell_count, y = median_predicted_proba, color = Metadata_Pathway)) +
    geom_point(aes(size = positive_mean_average_precision), alpha = 0.7) +
    geom_text_repel(data = subset(summary_data, high_count_and_proba == TRUE), aes(label = Metadata_treatment), 
                    size = 6, box.padding = 0.35, point.padding = 0.5, max.overlaps = 20, show.legend = FALSE) +  # Add labels
    scale_color_manual(values = custom_palette) +  # Use the custom color palette
    scale_size_continuous(name = "Positive mAP score") +
    labs(
        x = "Cell count\n(across plates)",
        y = "Median predicted probability\n(combined model)",
        color = "Treatment",
    ) +
    theme_bw() +
    theme(
        legend.position = "right",
        legend.title = element_text(size = 18),
        legend.text = element_text(size = 16),
        axis.text = element_text(size = 16),
        axis.title = element_text(size = 18)
    )

# Save the scatterplot
ggsave(paste0(figures_dir, "/batch1_pos_mAP_median_proba_vs_cell_count_compounds.png"), count_probas_plot, height = height, width = width, dpi = 500)

count_probas_plot
