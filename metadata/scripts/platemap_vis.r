suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(platetools))

# directory for fig to be outputted
output_dir <- file.path("./platemap_fig")

# Create the output directory if it doesn't exist
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# path to one platemap file since all have same layout but different compounds
platemap_file <- file.path("./target_screen_plate1.csv") 

# Directory for example platemap figure
output_fig <- file.path(output_dir, "example_platemap.png")

platemap_df <- readr::read_csv(
    platemap_file,
    col_types = readr::cols(.default = "c")
)

# Add new column for "condition" to plot on platemap for DMSO versus compound
platemap_df <- platemap_df %>%
    mutate(condition = ifelse(treatment == "DMSO", "DMSO", 
                              ifelse(grepl("^UCD", treatment), "compound", NA)))

print(dim(platemap_df))
head(platemap_df)

plate_replicate_gg <-
    platetools::raw_map(
        data = platemap_df$condition, # nolint
        well = platemap_df$well_position,
        plate = 96,
        size = 8
    ) +
    ggtitle("Platemap layout for all plates") +
    theme(plot.title = element_text(hjust = 0.75, size = 15, face = "bold")) +
    ggplot2::geom_point(aes(shape = platemap_df$cell_type)) +
    ggplot2::scale_shape_discrete(name = "Cell Type") +
    ggplot2::scale_fill_discrete(name = "Treatment")

ggsave(
    output_fig,
    plate_replicate_gg,
    dpi = 500,
    height = 3.5,
    width = 6
)

plate_replicate_gg
