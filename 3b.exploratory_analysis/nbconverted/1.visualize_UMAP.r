suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(arrow))
suppressPackageStartupMessages(library(ggExtra))

# Set up output directory for UMAP figures
dir.create("./figures", showWarnings = FALSE)

# Set directory and file structure
umap_dir <- "results"
umap_files <- list.files(umap_dir, full.names = TRUE)

output_fig_dir <- "figures"
plate_suffix <- ".parquet"

# Define output figure paths
output_umap_files <- setNames(
  file.path(
    output_fig_dir, 
    stringr::str_remove(basename(umap_files), plate_suffix) # Remove only .parquet
  ),
  basename(umap_files) # Use full original filenames as names
)

# Print the mapping in a cleaner format
cat("Mapping of input files to output paths:\n")
formatted_output <- data.frame(
  Original_File = basename(umap_files),
  Output_Path = file.path(output_fig_dir, stringr::str_remove(basename(umap_files), plate_suffix))
)
print(formatted_output, row.names = FALSE)

# Load data
umap_cp_df <- list()

for (plate in names(output_umap_files)) {
    # Find the umap file associated with the plate
    umap_file <- umap_files[stringr::str_detect(umap_files, plate)]
    
    if (length(umap_file) > 0) {
        # Load the umap data directly from Parquet file
        df <- arrow::read_parquet(umap_file)
         
        # Group by Metadata_Well and count cells
        cell_count_df <- df %>%
            dplyr::group_by(Metadata_Well) %>%
            dplyr::count() %>%
            dplyr::rename(Metadata_Cell_Count = n)
        
        # Merge the cell count data with the original dataframe
        umap_cp_df[[plate]] <- df %>%
            dplyr::left_join(cell_count_df, by = "Metadata_Well")
        
        # Update 'Endocrinology & Hormones' in Metadata_Pathway
        umap_cp_df[[plate]] <- umap_cp_df[[plate]] %>%
            dplyr::mutate(Metadata_Pathway = dplyr::recode(Metadata_Pathway,
                                                           "Endocrinology & Hormones" = "Endocrinology &\nHormones"))
            
    } else {
        message(paste("No file found for plate:", plate))
    }
}

# Inspect the first processed plate's data and print its dimensions
if (length(umap_cp_df) > 0) {
    plate_to_inspect <- names(umap_cp_df)[1]
    df_to_inspect <- umap_cp_df[[plate_to_inspect]]
    print(paste("Inspecting plate:", plate_to_inspect))
    print(paste("Dimensions:", dim(df_to_inspect)[1], "rows x", dim(df_to_inspect)[2], "columns"))
    head(df_to_inspect)
}


for (plate in names(umap_cp_df)) {
    # cell type UMAP
    output_file <- output_umap_files[[plate]]
    output_file <- paste0(output_file, "_cell_type.png")
    
    umap_dose_gg <- (
        ggplot(umap_cp_df[[plate]], aes(x = UMAP0, y = UMAP1))
        + geom_point(
            aes(color = Metadata_treatment_type), size = 0.4, alpha = 0.7
        )
        + facet_grid(Metadata_treatment_type ~ .)
        + theme_bw()
        + scale_color_brewer(palette = "Dark2", name = "Cell type")
        + theme(legend.position = "none")

    )
    
    ggsave(output_file, umap_dose_gg, dpi = 500, height = 6, width = 4)
}

for (plate in names(umap_cp_df)) {
    # Filter data for the two treatment types
    filtered_df <- umap_cp_df[[plate]] %>%
        dplyr::filter(Metadata_treatment_type %in% c("healthy + DMSO", "failing + DMSO"))
    
    # Generate output file path
    output_file <- output_umap_files[[plate]]
    output_file <- paste0(output_file, "_healthy_failing.png")
    
    # Create UMAP plot
    umap_gg <- (
        ggplot(filtered_df, aes(x = UMAP0, y = UMAP1))
        + geom_point(
            aes(color = Metadata_treatment_type), size = 0.4, alpha = 0.7
        )
        + theme_bw()
        + facet_wrap(~ Metadata_treatment_type, nrow = 1)
        + scale_color_manual(values = c("healthy + DMSO" = "#004400", "failing + DMSO" = "#a0004b"), name = "Treatment Type")
        + theme(legend.position = "none")
    )
    
    # Save the plot
    ggsave(output_file, umap_gg, dpi = 500, height = 4, width = 6)
}

custom_palette <- c(
  "Angiogenesis" = "#1b9e77",
  "Apoptosis" = "#d95f02",
  "DNA Damage" = "#7570b3",
  "Endocrinology &\nHormones" = "#e7298a",
  "Epigenetics" = "#66a61e",
  "MAPK" = "#e6ab02",
  "Metabolism" = "#a6761d",
  "Neuronal Signaling" = "#667665",
  "Others" = "#b3b3b3",
  "PI3K/Akt/mTOR" = "#8dd3c7",
  "Stem Cells & Wnt" = "#fb8072",
  "GPCR & G Protein" = "#984ea3",
  "healthy + DMSO" = "#004400",
  "failing + DMSO" = "#a0004b"
)

for (plate in names(umap_cp_df)) {
    # pathway UMAP
    output_file <- output_umap_files[[plate]]
    output_file <- paste0(output_file, "_pathway.png")

    # Move control facets to the front
    umap_cp_df[[plate]]$Metadata_Pathway <- factor(
        umap_cp_df[[plate]]$Metadata_Pathway,
        levels = c("healthy + DMSO", "failing + DMSO", setdiff(unique(umap_cp_df[[plate]]$Metadata_Pathway), c("healthy + DMSO", "failing + DMSO")))
        )

    umap_dose_gg <- (
        ggplot(umap_cp_df[[plate]], aes(x = UMAP0, y = UMAP1))
        + geom_point(
            aes(color = Metadata_Pathway), size = 0.4, alpha = 0.4
        )
        + theme_bw()
        + facet_wrap(~ Metadata_Pathway, nrow=2)
        + scale_color_manual(values = custom_palette) # Use the custom color palette
        + theme(legend.position = "none")
    )
    
    ggsave(output_file, umap_dose_gg, dpi = 500, height = 4, width = 10)
}
