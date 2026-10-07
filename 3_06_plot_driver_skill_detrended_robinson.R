
rm(list = ls())
graphics.off()

# ============================================================
# Supplementary Figure — Detrended driver prediction skill
#
# a) Temperature
# b) Relative Humidity
# c) Precipitation
# d) Wind Speed
#
# INIT vs ERA5, FY1–FY5
# Verification period: 1963–2022
#
# Reads pre-calculated detrended driver skill.
# NO recalculation is performed here.
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
  "verification_drivers_new",
  "results",
  "driver_skill_detrended_commontrend_FY1-FY5_corrected.rds"
)

out_dir <- file.path(
  PROJECT_ROOT,
  "plots",
  "revision",
  "Supplementary_driver_skill_detrended_commontrend_corrected"
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

vars <- c(
  "tas",
  "hurs",
  "pr",
  "sfcWind"
)

stopifnot(
  all(vars %in% names(x))
)

lon <- as.numeric(x$tas$lon)
lat <- as.numeric(x$tas$lat)

# ============================================================
# Geographic layers and palette
# ============================================================

layers <- load_robinson_base_layers()

pal <- correlation_palette()

z_limits <- c(-1, 1)
step_breaks <- seq(-1, 1, by = 0.2)

# ============================================================
# Panel a — Temperature
# ============================================================

p1 <- plot_europe_robinson_tiles(
  data_mat = as.matrix(
    x$tas$r_init_land
  ),

  sig_mask = as.matrix(
    x$tas$sig_init_land
  ),

  title_str =
    "a) Temperature",

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
# Panel b — Relative Humidity
# ============================================================

p2 <- plot_europe_robinson_tiles(
  data_mat = as.matrix(
    x$hurs$r_init_land
  ),

  sig_mask = as.matrix(
    x$hurs$sig_init_land
  ),

  title_str =
    "b) Relative Humidity",

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
# Panel c — Precipitation
# ============================================================

p3 <- plot_europe_robinson_tiles(
  data_mat = as.matrix(
    x$pr$r_init_land
  ),

  sig_mask = as.matrix(
    x$pr$sig_init_land
  ),

  title_str =
    "c) Precipitation",

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
# Panel d — Wind Speed
# ============================================================

p4 <- plot_europe_robinson_tiles(
  data_mat = as.matrix(
    x$sfcWind$r_init_land
  ),

  sig_mask = as.matrix(
    x$sfcWind$sig_init_land
  ),

  title_str =
    "d) Wind Speed",

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
    "FigS_detrended_driver_temperature.pdf"
  ),
  width = 8.2,
  height = 9.4
)

save_single_panel(
  p2,
  file.path(
    out_dir,
    "FigS_detrended_driver_humidity.pdf"
  ),
  width = 8.2,
  height = 9.4
)

save_single_panel(
  p3,
  file.path(
    out_dir,
    "FigS_detrended_driver_precipitation.pdf"
  ),
  width = 8.2,
  height = 9.4
)

save_single_panel(
  p4,
  file.path(
    out_dir,
    "FigS_detrended_driver_wind.pdf"
  ),
  width = 8.2,
  height = 9.4
)

# ============================================================
# Combined supplementary figure
# ============================================================

save_four_panel_figure(
  p1,
  p2,
  p3,
  p4,

  file.path(
    out_dir,
    "FigureS_detrended_driver_skill_INIT_vs_ERA5_commontrend_corrected.pdf"
  ),

  width = 15.5,
  height = 14.5
)

# ============================================================
# Numerical summary
# ============================================================

cat("\n====================================================\n")
cat("DETREND DRIVER SKILL FIGURE\n")
cat("====================================================\n")

for (var in vars) {

  r <- x[[var]]$r_init_land
  s <- x[[var]]$sig_init_land

  valid <- is.finite(r)

  cat("\n", var, "\n", sep = "")
  cat(
    "mean r      :",
    mean(r[valid], na.rm = TRUE),
    "\n"
  )
  cat(
    "median r    :",
    median(r[valid], na.rm = TRUE),
    "\n"
  )
  cat(
    "fraction >0 :",
    mean(r[valid] > 0, na.rm = TRUE),
    "\n"
  )
  cat(
    "significant :",
    mean(s[valid], na.rm = TRUE),
    "\n"
  )
}

cat("\nOutput:", out_dir, "\n")
cat("\nDONE\n")

