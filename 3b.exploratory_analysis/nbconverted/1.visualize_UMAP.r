suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(arrow))
suppressPackageStartupMessages(library(ggExtra))
source("./plot_themes.r")

# Set up output directory for UMAP figures
output_fig_dir <- file.path("figures", "UMAPs")
dir.create(output_fig_dir, recursive = TRUE, showWarnings = FALSE)

# Set directory and file structure
umap_dir <- "results"
umap_files <- list.files(umap_dir, pattern = "\\.parquet$", recursive = TRUE, full.names = TRUE)

# Build a small table of UMAP files with their platemap directories
umap_info <- data.frame(
  umap_file = umap_files,
  umap_basename = basename(umap_files),
  platemap = basename(dirname(umap_files)),
  stringsAsFactors = FALSE
)

# Make sure each platemap output directory exists
for (platemap_name in unique(umap_info$platemap)) {
  dir.create(file.path(output_fig_dir, platemap_name), recursive = TRUE, showWarnings = FALSE)
}

# Define output figure roots for each UMAP parquet file
umap_info$output_root <- file.path(
  output_fig_dir,
  umap_info$platemap,
  gsub("\\.parquet$", "", umap_info$umap_basename)
)
output_umap_files <- setNames(umap_info$output_root, umap_info$umap_file)

# Print the mapping in a cleaner format
cat("Mapping of input files to output paths:\n")
formatted_output <- data.frame(
  UMAP_File = umap_info$umap_file,
  Platemap = umap_info$platemap,
  Output_Path = umap_info$output_root,
  stringsAsFactors = FALSE
)
print(formatted_output, row.names = FALSE)

# Load data
umap_cp_df <- list()

for (umap_file in umap_info$umap_file) {
    if (file.exists(umap_file)) {
        # Load the umap data directly from Parquet file
        df <- arrow::read_parquet(umap_file)
         
        # Group by Metadata_Well and count cells
        cell_count_df <- df %>%
            dplyr::group_by(Metadata_Well) %>%
            dplyr::count() %>%
            dplyr::rename(Metadata_Cell_Count = n)
        
        # Merge the cell count data with the original dataframe
        key <- umap_file
        umap_cp_df[[key]] <- df %>%
            dplyr::left_join(cell_count_df, by = "Metadata_Well")
        
        # Update 'Endocrinology & Hormones' in Metadata_Pathway
        umap_cp_df[[key]] <- umap_cp_df[[key]] %>%
            dplyr::mutate(Metadata_Pathway = dplyr::recode(Metadata_Pathway,
                                                           "Endocrinology & Hormones" = "Endocrinology &\nHormones"))
            
    } else {
        message(paste("No file found:", umap_file))
    }
}

# Inspect the first processed plate's data and print its dimensions
if (length(umap_cp_df) > 0) {
    file_to_inspect <- names(umap_cp_df)[1]
    df_to_inspect <- umap_cp_df[[file_to_inspect]]
    print(paste("Inspecting file:", file_to_inspect))
    print(paste("Dimensions:", dim(df_to_inspect)[1], "rows x", dim(df_to_inspect)[2], "columns"))
    head(df_to_inspect)
}


umap_combined_bulk_files <- grep(
    "UMAP_combined_bulk",
    names(umap_cp_df),
    value = TRUE,
    ignore.case = TRUE
)

selected_key <- if (length(umap_combined_bulk_files) > 0) {
    umap_combined_bulk_files[1]
} else {
    message("No UMAP_combined_bulk profile found; inspecting the first available UMAP profile instead.")
    names(umap_cp_df)[1]
}

selected_df <- umap_cp_df[[selected_key]]

cat("Inspecting UMAP profile:", selected_key, "\n")
cat("Dimensions:", nrow(selected_df), "rows x", ncol(selected_df), "columns\n\n")

cat("Column names:\n")
print(names(selected_df))
cat("\nFirst 6 rows:\n")
print(head(selected_df))

cat("\nCounts by Metadata_treatment_type:\n")
print(selected_df %>%
    dplyr::count(Metadata_treatment_type) %>%
    dplyr::arrange(desc(n)))

cat("\nUnique Metadata_Pathway values (up to 10):\n")
print(unique(selected_df$Metadata_Pathway)[1:min(10, length(unique(selected_df$Metadata_Pathway)))])

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
    
    ggsave(output_file, umap_dose_gg, dpi = 500, height = 6, width = 3)
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

highlight_ids <- c(
  "UCD-0001921","UCD-0001812","UCD-0159268","UCD-0159406",
  "UCD-0000841","UCD-0159442","UCD-0001419","UCD-0159486","UCD-0159487"
)

for (plate in names(umap_cp_df)) {

    output_file <- output_umap_files[[plate]]
    output_file <- paste0(output_file, "_pathway.png")

    df <- umap_cp_df[[plate]]

    # detect bulk UMAP (for dot size scaling)
    is_bulk <- any(grepl("UMAP_combined_bulk", colnames(df)))

    # create combined color logic (highlights override everything else)
    df$color_group <- ifelse(
        df$Metadata_treatment %in% highlight_ids,
        df$Metadata_treatment,
        df$Metadata_Pathway
    )

    # determine which highlighted compounds are actually present
    present_highlights <- intersect(
        highlight_ids,
        unique(df$Metadata_treatment)
    )

    # Reorder Metadata_Pathway to put DMSO controls first
    df$Metadata_Pathway <- factor(
        df$Metadata_Pathway,
        levels = c(
            "healthy + DMSO",
            "failing + DMSO",
            setdiff(unique(df$Metadata_Pathway), c("healthy + DMSO", "failing + DMSO"))
        )
    )

    # dynamic point sizing
    point_size <- if (is_bulk) 2.8 else 0.3

    # Dynamic figure sizing based on number of facets
    n_facets <- length(unique(df$Metadata_Pathway))

    plot_height <- if (n_facets <= 13) {
        10
    } else if (n_facets <= 18) {
        12
    } else {
        14
    }

    plot_width <- 10

    umap_dose_gg <- (
        ggplot(df, aes(x = UMAP0, y = UMAP1))

        + geom_point(
            aes(color = color_group),
            size = point_size,
            alpha = 0.7
        )

        + facet_wrap(~Metadata_Pathway, ncol = 4)

        + theme_bw()

        + scale_color_manual(
            values = highlight_colors,
            breaks = present_highlights,
            name = "Hit compound(s)"
        )

        + guides(
            color = guide_legend(
                override.aes = list(
                    size = 4,
                    alpha = 1
                )
            )
        )

        + theme(
            legend.position = "right",
            legend.title = element_text(size = 10),
            legend.text = element_text(size = 8),
            strip.text = element_text(size = 10),
            strip.background = element_rect(fill = "grey90", color = NA),
            panel.spacing = unit(0.6, "lines")
        )
    )

    ggsave(
        output_file,
        umap_dose_gg,
        dpi = 600,
        height = plot_height,
        width = plot_width
    )
}
