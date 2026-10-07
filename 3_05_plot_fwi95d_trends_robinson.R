
rm(list = ls())
graphics.off()

# ============================================================
# Supplementary Figure — FWI95d linear trends
#
# a) ERA5 (OBS)
# b) INIT ensemble mean
# c) HIST+245 ensemble mean
#
# Verification period: 1963–2022
# Units: FWI95d days per decade
#
# Reads the already-calculated trend RDS.
# NO trend calculation is performed here.
# ============================================================

PROJECT_ROOT <-
  "/Users/marco/Library/CloudStorage/Dropbox/estcena/scripts/ONFIRE/decadal"

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
  "FWI95d_trends_1963-2022_corrected.rds"
)

out_dir <- file.path(
  PROJECT_ROOT,
  "plots",
  "revision",
  "Supplementary_FWI95d_trends_corrected"
)

dir.create(
  out_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

# ============================================================
# Read trends
# ============================================================

x <- readRDS(input_file)

stopifnot(
  all(
    c(
      "lon",
      "lat",
      "years",
      "trend_obs",
      "trend_init",
      "trend_hist"
    ) %in% names(x)
  ),
  all(x$years == 1963:2022)
)

lon <- as.numeric(x$lon)
lat <- as.numeric(x$lat)

trend_obs  <- as.matrix(x$trend_obs)
trend_init <- as.matrix(x$trend_init)
trend_hist <- as.matrix(x$trend_hist)

stopifnot(
  all(dim(trend_obs)  == c(length(lon), length(lat))),
  all(dim(trend_init) == c(length(lon), length(lat))),
  all(dim(trend_hist) == c(length(lon), length(lat)))
)

# ============================================================
# Common colour scale
#
# Use exactly the same symmetric scale for all three panels.
# The limit is derived from the largest absolute trend and
# rounded upward to the nearest integer.
# ============================================================

all_trends <- c(
  trend_obs,
  trend_init,
  trend_hist
)

zmax <- ceiling(
  max(abs(all_trends), na.rm = TRUE)
)

# Avoid an unnecessarily narrow scale
zmax <- max(zmax, 4)

z_limits <- c(-zmax, zmax)

step_breaks <- seq(
  -zmax,
  zmax,
  length.out = 9
)

cat("\nCommon trend scale:", -zmax, "to", zmax,
    "FWI95d days per decade\n")

# ============================================================
# Geographic layers
# ============================================================

layers <- load_robinson_base_layers()

# Same blue-white-red family as deterministic correlations:
# negative trends = blue
# positive trends = red
pal <- correlation_palette()

# ============================================================
# Panel a — ERA5
# ============================================================

p1 <- plot_europe_robinson_tiles(
  data_mat = trend_obs,
  sig_mask = NULL,

  title_str =
    "a) ERA5",

  lon = lon,
  lat = lat,
  layers = layers,

  z_limits = z_limits,
  step_breaks = step_breaks,
  palette = pal,

  title_size = 24,
  axis_text_size = 17,
  legend_text_size = 15
)

# ============================================================
# Panel b — INIT
# ============================================================

p2 <- plot_europe_robinson_tiles(
  data_mat = trend_init,
  sig_mask = NULL,

  title_str =
    "b) INIT",

  lon = lon,
  lat = lat,
  layers = layers,

  z_limits = z_limits,
  step_breaks = step_breaks,
  palette = pal,

  title_size = 24,
  axis_text_size = 17,
  legend_text_size = 15
)

# ============================================================
# Panel c — HIST+245
# ============================================================

p3 <- plot_europe_robinson_tiles(
  data_mat = trend_hist,
  sig_mask = NULL,

  title_str =
    "c) HIST+245",

  lon = lon,
  lat = lat,
  layers = layers,

  z_limits = z_limits,
  step_breaks = step_breaks,
  palette = pal,

  title_size = 24,
  axis_text_size = 17,
  legend_text_size = 15
)

# ============================================================
# Save individual panels
# ============================================================

save_single_panel(
  p1,
  file.path(
    out_dir,
    "FigS_FWI95d_trend_ERA5.pdf"
  ),
  width = 8.2,
  height = 9.4
)

save_single_panel(
  p2,
  file.path(
    out_dir,
    "FigS_FWI95d_trend_INIT.pdf"
  ),
  width = 8.2,
  height = 9.4
)

save_single_panel(
  p3,
  file.path(
    out_dir,
    "FigS_FWI95d_trend_HIST245.pdf"
  ),
  width = 8.2,
  height = 9.4
)

# ============================================================
# Combined supplementary figure
# ============================================================

save_three_panel_figure(
  p1,
  p2,
  p3,

  file.path(
    out_dir,
    "FigureS_FWI95d_trends_1963-2022_corrected.pdf"
  ),

  width = 15.5,
  height = 13.5
)

# ============================================================
# Summary
# ============================================================

cat("\n====================================================\n")
cat("FWI95d TREND FIGURE\n")
cat("====================================================\n")

if (!is.null(x$summary)) {
  print(x$summary)
}

cat("\nUnits  : FWI95d days per decade\n")
cat("Period :", paste(range(x$years), collapse = "-"), "\n")
cat("Scale  :", paste(z_limits, collapse = " to "), "\n")
cat("Output :", out_dir, "\n")

cat("\nDONE\n")

