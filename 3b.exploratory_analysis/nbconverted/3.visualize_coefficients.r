suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(dplyr))


# Set the output directory for figures
output_dir <- "figures/coefficient_plots"

# Create the directory if it doesn't exist
if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
}

cat("Output directory set to:", output_dir, "\n")


# instantiate a list to hold the coefficients
coefficients_dict <- list()

# load all CSV files paths from the "coefficients" folder
coefficients_folder <- "coefficients"
csv_files <- list.files(coefficients_folder, pattern = "\\.csv$", full.names = TRUE)

# read in coefficients per model
for (file in csv_files) {
    file_name <- basename(file)

    if (startsWith(file_name, "merged")) {
        next
    }

    df <- read.csv(file)

    if (!"Model" %in% names(df)) {
        warning(paste("Skipping", file, "- no 'Model' column found."))
        next
    }

    model_name <- unique(df$Model)

    if (length(model_name) != 1) {
        warning(paste("Skipping", file, "- multiple or no unique model names:", paste(model_name, collapse = ", ")))
        next
    }

    coefficients_dict[[model_name]] <- df
}

# print the names of the models
cat(names(coefficients_dict), sep = "\n")


# Process each model's coefficients
processed_dict <- lapply(coefficients_dict, function(df) {
    df %>%
        arrange(desc(abs(Coefficient))) %>%
        tidyr::separate(
            Feature,
            into = c(
                "compartment", "feature_group", "measurement",
                "channel", "parameter1", "parameter2", "parameter3"
            ),
            sep = "_",
            remove = FALSE,
            fill = "right"
        ) %>%
        mutate(
            Coefficient = abs(Coefficient),
            channel_cleaned = case_when(
                channel == "Hoechst" ~ "Nucleus",
                channel == "ER" ~ "ER",
                channel == "Actin" ~ "Actin",
                channel == "Mitochondria" ~ "Mito",
                channel == "PM" ~ "PM",
                TRUE ~ "other"
            )
        ) %>%
        filter(channel_cleaned %in% c("Mito", "Nucleus", "PM", "ER", "Actin", "other")) %>%
        group_by(feature_group, channel_cleaned, compartment) %>%
        slice_max(order_by = Coefficient, n = 1) %>%
        ungroup()
})

# Print one example cleaned output — the first model in the dict
example_name <- names(processed_dict)[1]
cat("Model:", example_name, "\n")
head(processed_dict[[example_name]], 5)


# Set channel order
channel_order <- c("Nucleus", "ER", "PM", "Mito", "Actin", "other")

# Set height and width of plot for visualizing
width <- 12
height <- 12
options(repr.plot.width = width, repr.plot.height = height)

# Loop through each model in the dictionary and make a plot
plots <- lapply(names(processed_dict), function(model_name) {
    df <- processed_dict[[model_name]] %>%
        mutate(
            rounded_coeff = round(Coefficient, 2),
            channel_cleaned = factor(channel_cleaned, levels = channel_order)
        )

    gg <- ggplot(df, aes(x = channel_cleaned, y = feature_group)) +
        geom_point(aes(fill = Coefficient), pch = 22, size = 14) +
        geom_text(aes(label = rounded_coeff), size = 4) +
        facet_wrap("~compartment", ncol = 3) +
        theme_bw() +
        scale_fill_distiller(
            name = "Top absolute\nvalue weight\nfrom model",
            palette = "YlGn",
            direction = 1
        ) +
        xlab("Channel") +
        ylab("Feature") +
        theme(
            axis.text = element_text(size = 12),
            axis.text.x = element_text(angle = 45, size = 13, vjust = 0.7, hjust = 0.5),
            axis.title = element_text(size = 14),
            legend.text = element_text(size = 12),
            legend.title = element_text(size = 14),
            strip.text = element_text(size = 14),
            strip.background = element_rect(colour = "black", fill = "#fdfff4"),
            legend.position = "right"
        )

    # Save coefficient plots
    ggsave(filename = paste0(output_dir, "/", model_name, "_coef_plot.png"), plot = gg, height = 5, width = 11.5, dpi = 500)

    return(gg)
})

# Print first plot as an example
plots[[1]]


# Load the merged coefficients file
merged_coefficients_file <- "coefficients/merged_coefficients_circ_model_combined_batch1.csv"

if (file.exists(merged_coefficients_file)) {
    merged_coefficients <- read.csv(merged_coefficients_file)
    cat("Merged coefficients file loaded successfully.\n")
} else {
    stop("File not found:", merged_coefficients_file)
}

# Split the feature column into parts to add back as columns
feature_parts <- strsplit(merged_coefficients$Feature, "_")

# Add feature parts as columns
merged_coefficients$compartment <- sapply(feature_parts, `[`, 1)
merged_coefficients$feature_group <- sapply(feature_parts, `[`, 2)
merged_coefficients$measurement <- sapply(feature_parts, `[`, 3)
merged_coefficients$organelle <- sapply(feature_parts, `[`, 4)
merged_coefficients$parameter1 <- sapply(feature_parts, `[`, 5)
merged_coefficients$parameter2 <- sapply(feature_parts, `[`, 6)
merged_coefficients$parameter3 <- sapply(feature_parts, `[`, 7)

# Replace invalid organelle values with "other"
merged_coefficients$organelle <- ifelse(
    is.na(merged_coefficients$organelle) |
        grepl("^[0-9]+$", merged_coefficients$organelle) |
        merged_coefficients$organelle %in% c("Adjacent", "X", "Y"),
    "other",
    merged_coefficients$organelle
)

# Display the first few rows of the loaded data
head(merged_coefficients)


# Set height and width of plot for visualizing
width <- 10
height <- 8
options(repr.plot.width = width, repr.plot.height = height)

# Generate scatterplot of coefficients
scatterplot_models <- ggplot(merged_coefficients, aes(
    x = Coefficient_orig_circ_model, y = Coefficient_combined_batch1_model,
    color = feature_group, shape = organelle
)) +
    geom_vline(xintercept = 0, linetype = "dotted", color = "black") +
    geom_hline(yintercept = 0, linetype = "dotted", color = "black") +
    geom_point(size = 4, alpha = 0.5) +
    theme_bw() +
    theme(
        text = element_text(size = 18),
        axis.title = element_text(size = 20),
        axis.text = element_text(size = 16),
        legend.title = element_text(size = 18),
        legend.text = element_text(size = 16)
    ) +
    labs(
        x = "Coefficient: Circulation model",
        y = "Coefficient: Combined batch 1 model",
        color = "Feature group",
        shape = "Organelle"
    )

# Save coefficient plots
ggsave(
    filename = paste0(output_dir, "/circ_model_combined_batch1_scatterplot.png"), plot = scatterplot_models,
    height = height, width = width, dpi = 500
)

scatterplot_models
