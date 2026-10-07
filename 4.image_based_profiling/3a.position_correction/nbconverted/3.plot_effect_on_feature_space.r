library(ggplot2)
library(dplyr)
library(readr)
library(cowplot)

figure_dir <- "figures"
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

read_result <- function(name) {
    readr::read_csv(name, show_col_types = FALSE)
}

tilt_sizes <- read_result("tilt_size_by_feature.csv")

family_panel <- function(column, title) {
    data <- tilt_sizes |>
        group_by(group = .data[[column]]) |>
        mutate(label = paste0(group, " (", n(), ")")) |>
        ungroup()
    data$label <- reorder(data$label, data$tilt_rms, FUN = median)
    (
        ggplot(data, aes(x = tilt_rms, y = label))
        + geom_boxplot(fill = "lightsteelblue", outlier.size = 0.6, outlier.alpha = 0.5)
        + labs(title = title, x = "Tilt size (RMS over plate positions)", y = NULL)
        + theme_bw(base_size = 13)
        + theme(plot.title.position = "plot")
    )
}

families <- c("compartment", "measurement", "channel")
titles <- c("By compartment", "By measurement type", "By channel")
n_groups <- sapply(families, function(column) n_distinct(tilt_sizes[[column]]))
panels <- Map(family_panel, families, titles)
figure <- cowplot::plot_grid(
    cowplot::ggdraw() + cowplot::draw_label(
        "Size of the plate-position tilt by feature family (number of features)",
        size = 15, x = 0.02, hjust = 0
    ),
    cowplot::plot_grid(plotlist = panels, ncol = 1, rel_heights = n_groups + 2, align = "v", axis = "lr"),
    ncol = 1, rel_heights = c(0.05, 1)
)
print(figure)
ggsave(
    file.path(figure_dir, "tilt_size_by_feature_family.png"), figure,
    dpi = 300, width = 7.5, height = 0.32 * sum(n_groups) + 4.5, bg = "white"
)

row_letters <- c("B", "C", "D", "E", "F", "G")
row_colors <- setNames(scales::viridis_pal(direction = -1)(6), row_letters)

scores <- read_result("well_pca_scores.csv") |>
    mutate(
        version = factor(version, levels = c("Before correction", "After correction")),
        well_row = factor(well_row, levels = row_letters)
    )
variance_by_row <- read_result("pc_variance_explained_by_row.csv") |>
    mutate(number = as.integer(sub("PC", "", component)))

row_components <- variance_by_row |>
    slice_max(row_explains_before, n = 2) |>
    arrange(number)
pc_x <- row_components$component[1]
pc_y <- row_components$component[2]

axis_label <- function(component) {
    sprintf("%s (%.0f%% of variance)", component, 100 * variance_by_row$variance_explained[variance_by_row$component == component])
}
row_share <- function(component, column) {
    100 * variance_by_row[[column]][variance_by_row$component == component]
}
panel_titles <- c(
    "Before correction" = sprintf(
        "Before correction\nplate row explains %.0f%% of %s and %.0f%% of %s",
        row_share(pc_x, "row_explains_before"), pc_x, row_share(pc_y, "row_explains_before"), pc_y
    ),
    "After correction" = sprintf(
        "After correction\nplate row explains %.0f%% of %s and %.0f%% of %s",
        row_share(pc_x, "row_explains_after"), pc_x, row_share(pc_y, "row_explains_after"), pc_y
    )
)
n_examples <- 100
displacement_title <- sprintf(
    "Displacement from before to after (zoomed)\nbold: mean of each plate row, thin: %d example wells", n_examples
)
panel_levels <- c(panel_titles, displacement_title)
scores <- scores |> mutate(panel = factor(panel_titles[as.character(version)], levels = panel_levels))

before <- scores |>
    filter(version == "Before correction") |>
    select(plate, well, well_row, bx = all_of(pc_x), by = all_of(pc_y))
after <- scores |>
    filter(version == "After correction") |>
    select(plate, well, ax = all_of(pc_x), ay = all_of(pc_y))
displacement <- inner_join(before, after, by = c("plate", "well")) |>
    mutate(panel = factor(displacement_title, levels = panel_levels))
row_means <- displacement |>
    group_by(panel, well_row) |>
    summarize(across(c(bx, by, ax, ay), mean), .groups = "drop")
set.seed(0)
examples <- slice_sample(displacement, n = n_examples)

view_limits <- function(values) {
    q <- quantile(values, c(0.01, 0.99))
    pad <- 0.08 * diff(q)
    c(q[1] - pad, q[2] + pad)
}
x_limits <- view_limits(scores[[pc_x]])
y_limits <- view_limits(scores[[pc_y]])
n_outside <- scores |>
    group_by(plate, well) |>
    summarize(
        outside = any(.data[[pc_x]] < x_limits[1] | .data[[pc_x]] > x_limits[2] |
            .data[[pc_y]] < y_limits[1] | .data[[pc_y]] > y_limits[2]),
        .groups = "drop"
    )

# the displacement panel has its own view, zoomed on the mean arrows
x_range <- range(c(row_means$bx, row_means$ax))
y_range <- range(c(row_means$by, row_means$ay))
half <- 0.65 * max(diff(x_range), diff(y_range))
x_center <- mean(x_range)
y_center <- mean(y_range)

control_label <- "DMSO control\n(outline)"
wells_plot <- (
    ggplot(scores, aes(x = .data[[pc_x]], y = .data[[pc_y]]))
    + geom_point(data = filter(scores, !is_control), aes(color = well_row), size = 1, alpha = 0.85)
    + geom_point(
        data = filter(scores, is_control), aes(fill = well_row, shape = control_label),
        color = "black", size = 2.2, stroke = 0.4
    )
    + facet_wrap(vars(panel))
    + scale_color_manual(values = row_colors, name = "Plate row")
    + scale_fill_manual(values = row_colors, guide = "none")
    + scale_shape_manual(values = setNames(21, control_label), name = NULL)
    + guides(
        color = guide_legend(order = 1, override.aes = list(size = 3)),
        shape = guide_legend(order = 2, override.aes = list(fill = "white", color = "black", size = 3))
    )
    + coord_cartesian(xlim = x_limits, ylim = y_limits)
    + labs(
        x = axis_label(pc_x), y = axis_label(pc_y),
        caption = sprintf("%d of %d wells are outside the view", sum(n_outside$outside), nrow(n_outside))
    )
    + theme_bw(base_size = 13)
    + theme(legend.position = "right")
)
displacement_plot <- (
    ggplot()
    + geom_segment(
        data = examples, aes(x = bx, y = by, xend = ax, yend = ay, color = well_row),
        arrow = arrow(length = unit(0.08, "cm"), type = "closed"), alpha = 0.3, linewidth = 0.3
    )
    + geom_segment(
        data = row_means, aes(x = bx, y = by, xend = ax, yend = ay, color = well_row),
        arrow = arrow(length = unit(0.25, "cm"), type = "closed"), linewidth = 1.3
    )
    + geom_label(
        data = row_means, aes(x = ax, y = ay, label = well_row),
        nudge_x = -0.06 * half, nudge_y = 0.06 * half, fontface = "bold", size = 4,
        label.padding = unit(0.12, "lines"), alpha = 0.85
    )
    + facet_wrap(vars(panel))
    + scale_color_manual(values = row_colors, name = "Plate row")
    + coord_cartesian(xlim = x_center + c(-half, half), ylim = y_center + c(-half, half))
    + labs(x = axis_label(pc_x), y = NULL)
    + theme_bw(base_size = 13)
    + theme(legend.position = "none")
)

legend <- cowplot::get_legend(wells_plot)
figure <- cowplot::plot_grid(
    cowplot::ggdraw() + cowplot::draw_label(
        sprintf("Wells in %s and %s of the uncorrected wells, before and after the correction", pc_x, pc_y),
        size = 15, x = 0.01, hjust = 0
    ),
    cowplot::plot_grid(
        wells_plot + theme(legend.position = "none"), displacement_plot, legend,
        nrow = 1, rel_widths = c(2, 1, 0.45), align = "h", axis = "tb"
    ),
    ncol = 1, rel_heights = c(0.06, 1)
) + theme(plot.background = element_rect(fill = "white", color = NA))
print(figure)
ggsave(file.path(figure_dir, "well_pca_before_after_by_row.png"), figure, dpi = 300, width = 20, height = 6.5)


top_component <- variance_by_row |> slice_max(row_explains_before, n = 1) |> pull(component)
platemap_levels <- c(paste("Platemap", sort(unique(scores$platemap_number))), "All platemaps")

plate_grid <- bind_rows(
    scores |> mutate(platemap_label = paste("Platemap", platemap_number)),
    scores |> mutate(platemap_label = "All platemaps")
) |>
    group_by(version, platemap_label, well_row, well_col) |>
    summarize(score = mean(.data[[top_component]]), .groups = "drop") |>
    mutate(
        platemap_label = factor(platemap_label, levels = platemap_levels),
        well_row = factor(well_row, levels = rev(row_letters)),
        well_col = factor(well_col, levels = 2:11)
    )
limit <- quantile(abs(filter(plate_grid, version == "Before correction")$score), 0.98, na.rm = TRUE)

figure <- (
    ggplot(plate_grid, aes(x = well_col, y = well_row, fill = score))
    + geom_tile()
    + facet_grid(rows = vars(version), cols = vars(platemap_label))
    + scale_fill_gradient2(
        low = "#2166AC", mid = "white", high = "#B2182B", limits = c(-limit, limit),
        oob = scales::squish, na.value = "gray85", name = paste("Mean", top_component)
    )
    + scale_x_discrete(breaks = c(2, 11))
    + labs(
        title = sprintf("%s at each plate position, before and after the correction", top_component),
        x = "Plate column", y = "Plate row"
    )
    + theme_minimal(base_size = 12)
    + theme(panel.grid = element_blank(), axis.text = element_text(size = 8), plot.title.position = "plot")
)
print(figure)
ggsave(file.path(figure_dir, "plate_layout_top_row_pc_before_after.png"), figure, dpi = 300, width = 16, height = 4.6)
