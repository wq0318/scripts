#!/usr/bin/env Rscript
# Plot 5 spatial tile SVGs from generate_nfrags_grid.py output CSV.
# Coordinate range fixed at xmin=0, xmax=23520, ymin=0, ymax=23520.
# Usage: Rscript plot_spatial_tiles.R <grid.csv> <outdir> <prefix>

suppressPackageStartupMessages({
    library(data.table)
    library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
    stop("Usage: Rscript plot_spatial_tiles.R <grid.csv> <outdir> <prefix>")
}

csv_path <- args[1]
outdir <- args[2]
prefix <- args[3]

if (!file.exists(csv_path)) {
    stop(sprintf("Input CSV not found: %s", csv_path))
}
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

dt <- fread(csv_path)
required_cols <- c("CB", "x", "y", "in_tissue",
                   "nFrags_filtered", "nFrags_chrM",
                   "nFrags_blacklist", "nFrags_raw")
missing_cols <- setdiff(required_cols, colnames(dt))
if (length(missing_cols) > 0) {
    stop(sprintf("Missing columns in %s: %s",
                 csv_path, paste(missing_cols, collapse = ", ")))
}

dt[, in_tissue := tolower(as.character(in_tissue)) == "true"]

XMIN <- 0
XMAX <- 23520
YMIN <- 0
YMAX <- 23520

# Palette aligned to generate_nfrags_grid.py custom_colorscale (16 stops)
fill_palette <- c(
    "#0E458F", "#0F5298", "#0E6BA8", "#0C86B8", "#3399A1", "#3B9C9C",
    "#B2C061", "#F2CE38", "#F2AB38", "#F2AB38", "#EB7232", "#E65B2E",
    "#E14428", "#DC2E22", "#DB2921", "#CC2623"
)

plot_tile <- function(data, value_col, title, outfile) {
    if (nrow(data) == 0) {
        message(sprintf("[plot_tile] No data for %s -- writing blank SVG", outfile))
        blank <- ggplot() +
            annotate("text", x = 0.5, y = 0.5,
                     label = paste0("No data: ", title), colour = "white") +
            theme_void() +
            theme(plot.background  = element_rect(fill = "black", colour = NA),
                  panel.background = element_rect(fill = "black", colour = NA))
        ggsave(outfile, blank, width = 12, height = 8, device = "svg")
        return(invisible(NULL))
    }
    p <- ggplot(data, aes(x = x, y = y, fill = .data[[value_col]])) +
        geom_tile() +
        scale_fill_gradientn(
            colours = fill_palette,
            name = "nFrags",
            na.value = "black"
        ) +
        coord_fixed(xlim = c(XMIN, XMAX), ylim = c(YMIN, YMAX), expand = FALSE) +
        labs(title = title, x = NULL, y = NULL) +
        theme_void() +
        theme(
            plot.background   = element_rect(fill = "black", colour = NA),
            panel.background  = element_rect(fill = "black", colour = NA),
            legend.background = element_rect(fill = "black", colour = NA),
            legend.key        = element_rect(fill = "black", colour = NA),
            text              = element_text(colour = "white", size = 14),
            plot.title        = element_text(hjust = 0.5, colour = "white"),
            legend.text       = element_text(colour = "white"),
            legend.title      = element_text(colour = "white"),
            legend.position   = "right"
        )
    # width=12, height=8 keeps the 23520x23520 square plot panel (constrained
    # by coord_fixed) and reserves ~4 inches of width for the legend outside
    # the data area. Browser scales preserve aspect, producing a ~2/3 visual
    # height vs the previous 8x8 (square).
    ggsave(outfile, p, width = 12, height = 8, device = "svg")
    message(sprintf("[plot_tile] saved %s (n=%d)", outfile, nrow(data)))
}

out <- function(suffix) file.path(outdir, paste0(prefix, suffix))

plot_tile(dt[nFrags_filtered > 0],
          "nFrags_filtered",
          "All bins - nFrags (filtered)",
          out("_all_bins_nFrags_filtered.svg"))

plot_tile(dt[nFrags_raw > 0],
          "nFrags_raw",
          "All bins - nFrags (raw)",
          out("_all_bins_nFrags_raw.svg"))

plot_tile(dt[nFrags_chrM > 0],
          "nFrags_chrM",
          "All bins - nFrags (chrM)",
          out("_all_bins_nFrags_chrM.svg"))

plot_tile(dt[nFrags_blacklist > 0],
          "nFrags_blacklist",
          "All bins - nFrags (blacklist)",
          out("_all_bins_nFrags_blacklist.svg"))

plot_tile(dt[in_tissue == TRUE & nFrags_filtered > 0],
          "nFrags_filtered",
          "Bins under tissue - nFrags (filtered)",
          out("_bins_under_tissue_nFrags_filtered.svg"))

message("[plot_spatial_tiles] All 5 SVGs written to ", outdir)
