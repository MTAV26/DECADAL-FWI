rm(list = ls())
graphics.off()

# ============================================================
# Figure 3 — Detrended deterministic skill
# Exact analogue of Figure 2 after linear detrending
# Verification period: 1963–2022
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
  "Figure3_detrended_deterministic_skill_commontrend_corrected"
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
      "corr",
      "residual",
      "diffcorr"
    ) %in% names(x)
  ),
  all(x$years == 1963:2022)
)

lon <- as.numeric(x$lon)
lat <- as.numeric(x$lat)

corr_mat <- as.matrix(x$corr$r)
corr_sig <- as.matrix(x$corr$sign)

res_mat <- as.matrix(x$residual$res.corr)
res_sig <- as.matrix(x$residual$sign)

diff_mat <- as.matrix(x$diffcorr$diff.corr)
diff_sig <- as.matrix(x$diffcorr$sign)

stopifnot(
  all(dim(corr_mat) == c(length(lon), length(lat))),
  all(dim(res_mat)  == c(length(lon), length(lat))),
  all(dim(diff_mat) == c(length(lon), length(lat)))
)

# ============================================================
# Geographic layers and palette
# ============================================================

layers <- load_robinson_base_layers()
pal <- correlation_palette()

# ============================================================
# Panel a
# ============================================================

p1 <- plot_europe_robinson_tiles(
  data_mat = corr_mat,
  sig_mask = corr_sig,

  title_str =
    "a) Detrended Correlation (INIT vs OBS)",

  lon = lon,
  lat = lat,
  layers = layers,

  z_limits = c(-1, 1),
  step_breaks = seq(-1, 1, by = 0.2),
  palette = pal,

  title_size = 22,
  axis_text_size = 17,
  legend_text_size = 16
)

# ============================================================
# Panel b
# ============================================================

p2 <- plot_europe_robinson_tiles(
  data_mat = res_mat,
  sig_mask = res_sig,

  title_str =
    "b) Detrended Residual Correlation",

  lon = lon,
  lat = lat,
  layers = layers,

  z_limits = c(-1, 1),
  step_breaks = seq(-1, 1, by = 0.2),
  palette = pal,

  title_size = 22,
  axis_text_size = 17,
  legend_text_size = 16
)

# ============================================================
# Panel c
# Slightly smaller title because this is the longest panel title
# ============================================================

p3 <- plot_europe_robinson_tiles(
  data_mat = diff_mat,
  sig_mask = diff_sig,

  title_str =
    "c) Detrended Difference in Correlation\n(INIT - HIST+245)",

  lon = lon,
  lat = lat,
  layers = layers,

  z_limits = c(-1, 1),
  step_breaks = seq(-1, 1, by = 0.2),
  palette = pal,

  title_size = 20,
  axis_text_size = 17,
  legend_text_size = 16
)

# ============================================================
# Individual panels
# Slightly wider than before
# ============================================================

save_single_panel(
  p1,
  file.path(
    out_dir,
    "Fig3a_detrended_INIT_vs_OBS.pdf"
  ),
  width = 8.2,
  height = 9.4
)

save_single_panel(
  p2,
  file.path(
    out_dir,
    "Fig3b_detrended_residual_correlation.pdf"
  ),
  width = 8.2,
  height = 9.4
)

save_single_panel(
  p3,
  file.path(
    out_dir,
    "Fig3c_detrended_difference_INIT_minus_HIST245.pdf"
  ),
  width = 8.2,
  height = 9.4
)

# ============================================================
# Combined Figure 3
#
# Wider canvas than Figure 2 original so that panel c title
# is fully visible without enlarging the artboard in Illustrator.
# ============================================================

save_three_panel_figure(
  p1,
  p2,
  p3,

  file.path(
    out_dir,
    "Figure3_detrended_deterministic_skill_commontrend_corrected.pdf"
  ),

  width = 15.5,
  height = 13.5
)

cat("\n====================================================\n")
cat("DONE — FIGURE 3\n")
cat("====================================================\n")
cat("Input  :", input_file, "\n")
cat("Years  :", paste(range(x$years), collapse = "-"), "\n")
cat("Output :", out_dir, "\n")
cat("====================================================\n")
