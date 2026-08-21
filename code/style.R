# =============================================================================
# field-fnirs-soundscape — shared publication plot style
# Sourced by every figure script. PNG 600 dpi, opaque white.
# Paths resolve relative to the repository root (run everything from there).
# =============================================================================
library(ggplot2)

WIDTH_MM <- c(single_column = 85, double_column = 178)
mm2in <- function(mm) mm / 25.4

okabe_ito <- c("#0072B2", "#E69F00", "#009E73", "#D55E00", "#56B4E9", "#CC79A7")

# Consistent side mapping across every figure
COL_NATURAL   <- "#009E73"
COL_COMPOSITE <- "#E69F00"

theme_pub <- function(base_size = 9, axis_title_size = 10) {
  theme_classic(base_size = base_size, base_family = "Helvetica") %+replace%
    theme(
      axis.line         = element_line(colour = "black", linewidth = 0.5),
      axis.ticks        = element_line(colour = "black", linewidth = 0.4),
      axis.ticks.length = unit(0.10, "cm"),
      axis.title        = element_text(size = axis_title_size),
      axis.title.x      = element_text(size = axis_title_size,
                                       margin = margin(t = 5)),
      axis.title.y      = element_text(size = axis_title_size, angle = 90,
                                       margin = margin(r = 6)),
      axis.text         = element_text(size = base_size, colour = "black"),
      legend.text       = element_text(size = base_size),
      legend.title      = element_blank(),
      legend.position   = "inside",
      legend.position.inside = c(0.98, 0.98),
      legend.justification.inside = c(1, 1),
      legend.background = element_blank(),
      legend.key        = element_blank(),
      panel.grid        = element_blank(),
      plot.title        = element_blank(),
      plot.background   = element_rect(fill = "white", colour = NA),
      panel.background  = element_rect(fill = "white", colour = NA)
    )
}

save_fig <- function(plot, path, width_mm, height_mm) {
  ggsave(path, plot, width = mm2in(width_mm), height = mm2in(height_mm),
         dpi = 600, bg = "white", device = ragg::agg_png)
}

FIGDIR  <- file.path("output", "figures")
LOCKDIR <- file.path(FIGDIR, "data_lock")
# Analysis modules write to output/aNN_<slug>/; display scripts read the same
# files. The mapping keeps the working pipeline's exploration folder names
# ("analysis_05_fnirs_signal_anatomy" -> "output/a05_fnirs_signal_anatomy").
OUTC    <- function(a, f) file.path("output", sub("^analysis_", "a", a), f)
dir.create(LOCKDIR, showWarnings = FALSE, recursive = TRUE)
