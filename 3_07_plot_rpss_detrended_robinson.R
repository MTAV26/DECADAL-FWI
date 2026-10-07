rm(list = ls())
graphics.off()

# ============================================================
# Figure S4 — Detrended probabilistic skill (RPSS)
# Verification period: 1963–2022
#
# a) INIT vs climatology
# b) INIT vs HIST+245
#
# Reads pre-calculated detrended results from 2_03.
# NO RPSS calculation is performed here.
# ============================================================

PROJECT_ROOT <- "/Users/marco/Library/CloudStorage/Dropbox/estcena/scripts/ONFIRE/decadal"

source(
  file.path(
    PROJECT_ROOT,
    "scripts",
    "3_figures",
    "plot_functions_robinson.R"
  )
)

input_file <- file.path(
  PROJECT_ROOT,
  "data",
  "verification_new",
  "results",
  "skill_detrended_commontrend_corrected.rds"
)

out_dir <- file.path(
  PROJECT_ROOT,
  "plots",
  "revision",
  "FigureS4_RPSS_detrended"
)

dir.create(
  out_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

# ============================================================
# Read results
# ============================================================

x <- readRDS(input_file)

stopifnot(
  all(
    c(
      "lon",
      "lat",
      "years",
      "clim_years",
      "rpss_clim",
      "rpss_hist"
    ) %in% names(x)
  ),
  all(x$years == 1963:2022)
)

lon <- as.numeric(x$lon)
lat <- as.numeric(x$lat)

rpss_clim <- as.matrix(x$rpss_clim$rpss)
sig_clim  <- as.matrix(x$rpss_clim$sign)

rpss_hist <- as.matrix(x$rpss_hist$rpss)
sig_hist  <- as.matrix(x$rpss_hist$sign)

stopifnot(
  all(dim(rpss_clim) == c(length(lon), length(lat))),
  all(dim(rpss_hist) == c(length(lon), length(lat))),
  all(dim(sig_clim) == dim(rpss_clim)),
  all(dim(sig_hist) == dim(rpss_hist))
)

# ============================================================
# Geographic layers and palette
# ============================================================

layers <- load_robinson_base_layers()
pal <- rpss_palette()

# ============================================================
# Panel a — INIT vs climatology
# ============================================================

p1 <- plot_europe_robinson_tiles(
  data_mat = rpss_clim,
  sig_mask = sig_clim,
  title_str =
    "a) RPSS (INIT vs Climatology)",
  lon = lon,
  lat = lat,
  layers = layers,
  z_limits = c(-0.5, 0.5),
  step_breaks = seq(-0.5, 0.5, by = 0.1),
  palette = pal,
  title_size = 23,
  axis_text_size = 17,
  legend_text_size = 16
)

# ============================================================
# Panel b — INIT vs HIST+245
# ============================================================

p2 <- plot_europe_robinson_tiles(
  data_mat = rpss_hist,
  sig_mask = sig_hist,
  title_str =
    "b) RPSS (INIT vs HIST+245)",
  lon = lon,
  lat = lat,
  layers = layers,
  z_limits = c(-0.5, 0.5),
  step_breaks = seq(-0.5, 0.5, by = 0.1),
  palette = pal,
  title_size = 23,
  axis_text_size = 17,
  legend_text_size = 16
)

# ============================================================
# Save individual panels
# ============================================================

save_single_panel(
  p1,
  file.path(
    out_dir,
    "FigS4a_RPSS_INIT_vs_climatology_detrended.pdf"
  ),
  width = 8.2,
  height = 9.4
)

save_single_panel(
  p2,
  file.path(
    out_dir,
    "FigS4b_RPSS_INIT_vs_HIST245_detrended.pdf"
  ),
  width = 8.2,
  height = 9.4
)

# ============================================================
# Combined Figure S4
# ============================================================

save_two_panel_figure(
  p1,
  p2,
  file.path(
    out_dir,
    "FigureS4_RPSS_detrended.pdf"
  ),
  width = 15.5,
  height = 8.2
)

# ============================================================
# Brief numerical summary
# ============================================================

land_clim <- is.finite(rpss_clim)
land_hist <- is.finite(rpss_hist)

cat("\n====================================================\n")
cat("FIGURE S4 — DETRENDED RPSS SUMMARY\n")
cat("====================================================\n")

cat("\nINIT vs climatology\n")
cat("Mean RPSS   :",
    mean(rpss_clim[land_clim], na.rm = TRUE), "\n")
cat("Median RPSS :",
    median(rpss_clim[land_clim], na.rm = TRUE), "\n")
cat("RPSS > 0    :",
    mean(rpss_clim[land_clim] > 0, na.rm = TRUE), "\n")
cat("Significant :",
    mean(sig_clim[land_clim], na.rm = TRUE), "\n")

cat("\nINIT vs HIST+245\n")
cat("Mean RPSS   :",
    mean(rpss_hist[land_hist], na.rm = TRUE), "\n")
cat("Median RPSS :",
    median(rpss_hist[land_hist], na.rm = TRUE), "\n")
cat("RPSS > 0    :",
    mean(rpss_hist[land_hist] > 0, na.rm = TRUE), "\n")
cat("Significant :",
    mean(sig_hist[land_hist], na.rm = TRUE), "\n")

cat("\n====================================================\n")
cat("DONE — FIGURE S4\n")
cat("====================================================\n")
cat("Input       :", input_file, "\n")
cat("Years       :", paste(range(x$years), collapse = "-"), "\n")
cat("Climatology :", paste(range(x$clim_years), collapse = "-"), "\n")
cat("Output      :", out_dir, "\n")
cat("====================================================\n")
