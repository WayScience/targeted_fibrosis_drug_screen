suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(platetools))
suppressPackageStartupMessages(library(patchwork))

# Set the name of batch being processed
batch_name <- "platemap_11"

# Output directory for the plots
output_directory <- file.path("./qc_plots", batch_name)
# Create the output directory if it doesn't exist
if (!dir.exists(output_directory)) {
  dir.create(output_directory, recursive = TRUE)
}

# Define the directory containing the folders for each plate
qc_results <- paste0("./qc_results/", batch_name)

# List all CSV files in the directory and its subdirectories
csv_files <- list.files(qc_results, pattern = "\\Image.csv$", full.names = TRUE, recursive = TRUE)

# Read and concatenate all CSV files, updating Metadata_Plate to match folder name
illum_data <- bind_rows(lapply(csv_files, function(file) {
  df <- read.csv(file)
  
  # Extract plate name from folder
  plate_name <- basename(dirname(file))
  
  # Warn if Metadata_Plate is inconsistent
  if (!all(df$Metadata_Plate == plate_name)) {
    warning(paste0("Metadata_Plate mismatch in file: ", file, 
                   ". Updating to folder name '", plate_name, "'"))
  }
  
  # Update Metadata_Plate to match folder
  df$Metadata_Plate <- plate_name
  
  return(df)
}))

# Display the first few rows of the concatenated data
head(illum_data)


well_qc_data <- illum_data %>%
  mutate(Failing_FOV = (Metadata_Blur_Flag == 1 | Metadata_Saturation_Flag == 1)) %>%
  group_by(Metadata_Plate, Metadata_Well) %>%
  summarise(
    Total_FOVs = n(),
    Count_Failing = sum(Failing_FOV),
    Percent_Failing = 100 * Count_Failing / Total_FOVs,
    .groups = "drop"
  )

# Add Metadata_treatment information based on Metadata_Well
control_wells <- c("B02", "B05", "B08", "B11", "E02", "E05", "E08", "E11")
well_qc_data <- well_qc_data %>%
  mutate(Metadata_treatment = ifelse(Metadata_Well %in% control_wells, "control", "compound"))

head(well_qc_data)


library(dplyr)
library(stringr)

# Use the actual plate names in well_qc_data
unique_plates <- unique(well_qc_data$Metadata_Plate)

# Extract numeric part for sorting (handles leading zeros)
plate_numbers <- str_extract(unique_plates, "(?<=CARD-CelIns-CX7)\\d+") %>%
    as.numeric()

# Order by numeric value
sorted_indices <- order(plate_numbers)
old_names <- unique_plates[sorted_indices]

# New names
new_names <- paste0("Plate_", seq_along(old_names))

# Create a mapping and update well_qc_data
plate_mapping <- setNames(new_names, old_names)

well_qc_data <- well_qc_data %>%
    mutate(Metadata_Plate = recode(Metadata_Plate, !!!plate_mapping))

# Check
unique(well_qc_data$Metadata_Plate)
head(well_qc_data)


# Generate platemaps for each Plate
unique_plates <- unique(well_qc_data$Metadata_Plate)

for (plate in unique_plates) {
    plate_data <- well_qc_data %>%
        filter(Metadata_Plate == plate)

    fov_platemap <- platetools::raw_map(
        data = plate_data$Percent_Failing,
        well = plate_data$Metadata_Well,
        plate = 96,
        size = 8
    ) +
        ggtitle(paste(plate, "(25 FOVs per well)")) +
        theme(plot.title = element_text(size = 10, face = "bold")) +
        scale_fill_gradientn(
            name = "Percent failing\nFOVs",
            colors = c("#1a9850", "#fee08b", "#d73027"), # green → yellow → red
            values = scales::rescale(c(0, 40, 100)),
            limits = c(0, 100)
        ) +
        geom_point(aes(shape = plate_data$Metadata_treatment), size = 2.5) +
        scale_shape_discrete(name = "Treatment")

    print(fov_platemap)
}


# Generate platemaps for each Plate
unique_plates <- unique(well_qc_data$Metadata_Plate)

# Create a list to store the plots
plate_plots <- lapply(unique_plates, function(plate) {
    plate_data <- well_qc_data %>%
        filter(Metadata_Plate == plate)

    platemap <- platetools::raw_map(
        data = plate_data$Percent_Failing,
        well = plate_data$Metadata_Well,
        plate = 96,
        size = 10
    ) +
        ggtitle(paste(plate, "(25 FOVs per well)")) +
        theme(plot.title = element_text(size = 10, face = "bold")) +
        scale_fill_gradientn(
            name = "Percent failing\nFOVs",
            colors = c("#1a9850", "#fee08b", "#d73027"),
            values = scales::rescale(c(0, 40, 100)),
            limits = c(0, 100)
        )+
        geom_point(aes(shape = plate_data$Metadata_treatment), size = 2.5) +
        scale_shape_discrete(name = "Treatment" ) +
        theme(
            plot.title = element_text(size = 16),
            legend.title = element_text(size = 14),
            legend.text = element_text(size = 12),
            axis.text = element_text(size = 12),
            axis.title = element_text(size = 14),
        )
    
    return(platemap)
})

# Combine all plots in a 2x2 grid
combined_plot <- wrap_plots(plate_plots, ncol = 2)
combined_plot

# Save the combined plot to a file
ggsave(filename = file.path(output_directory, paste0(batch_name, "_per_plate_fov_platemaps.png")),
       plot = combined_plot,
       width = 15, height = 9, dpi = 600)


merged_well_qc_data <- illum_data %>%
  mutate(Failing_FOV = (Metadata_Blur_Flag == 1 | Metadata_Saturation_Flag == 1)) %>%
  group_by(Metadata_Well) %>%
  summarise(
    Total_FOVs = n(),
    Count_Failing = sum(Failing_FOV),
    Percent_Failing = 100 * Count_Failing / Total_FOVs,
    .groups = "drop"
  )

head(merged_well_qc_data)


# Extract the plate number from your batch_name
plate_number <- sub("platemap_", "", batch_name)

# Construct the correct platemap filename
platemap_file <- paste0("../../metadata/original_platemaps/Target_Selective_Library_Screen_Plate_", plate_number, ".csv")

# Read in the platemap
platemap_data <- read.csv(platemap_file)

merged_well_qc_data <- illum_data %>%
  mutate(Failing_FOV = (Metadata_Blur_Flag == 1 | Metadata_Saturation_Flag == 1)) %>%
  group_by(Metadata_Well) %>%
  summarise(
    Total_FOVs = n(),
    Count_Failing = sum(Failing_FOV),
    Percent_Failing = 100 * Count_Failing / Total_FOVs,
    .groups = "drop"
  )

# Merge with the well QC data via well_position and only include treatment column
merged_well_qc_data <- merged_well_qc_data %>%
  inner_join(platemap_data %>% select(well_position, treatment), 
             by = c("Metadata_Well" = "well_position"))

# Update treatment column to be either control if DMSO or if UCD- prefix then compound
merged_well_qc_data <- merged_well_qc_data %>%
  mutate(treatment = ifelse(treatment == "DMSO", "control",
                            ifelse(startsWith(treatment, "UCD-"), "compound", treatment)))

head(merged_well_qc_data)


merged_fov_platemap <- platetools::raw_map(
    data = merged_well_qc_data$Percent_Failing,
    well = merged_well_qc_data$Metadata_Well,
    plate = 96,
    size = 12
) +
    ggtitle(paste("All plates in batch (25 FOVs per well)")) +
    theme(plot.title = element_text(size = 10, face = "bold")) +
    scale_fill_gradientn(
        name = "Percent failing\nFOVs",
        colors = c("#1a9850", "#fee08b", "#d73027"), # green → yellow → red
        values = scales::rescale(c(0, 40, 100)),
        limits = c(0, 100)
    ) +
    geom_point(aes(shape = merged_well_qc_data$treatment), size = 2.5) +
    scale_shape_discrete(name = "Treatment", guide = guide_legend(override.aes = list(size = 2.5)))

# Save the merged plot to the output directory, including the batch name in the filename
output_file_merged <- file.path(output_directory, paste0("merged_", batch_name, "_fov_platemap.png"))
ggsave(output_file_merged, plot = merged_fov_platemap, width = 8, height = 6)

print(merged_fov_platemap)

