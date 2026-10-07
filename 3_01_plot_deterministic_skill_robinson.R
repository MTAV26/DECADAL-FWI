rm(list = ls())
graphics.off()

# Figure 2 — Raw deterministic skill
# Reads the already-calculated corrected RDS.
# NO skill calculation is performed here.

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
  "deterministic_skill_corrected.rds"
)

out_dir <- file.path(
  PROJECT_ROOT,
  "plots",
  "revision",
  "Figure2_raw_deterministic_skill_corrected"
)

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

x <- readRDS(input_file)

stopifnot(
  all(c("lon", "lat", "years", "corr", "residual", "diffcorr") %in% names(x)),
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

layers <- load_robinson_base_layers()

pal <- correlation_palette()

p1 <- plot_europe_robinson_tiles(
  data_mat = corr_mat,
  sig_mask = corr_sig,
  title_str = "a) Correlation (INIT vs OBS)",
  lon = lon,
  lat = lat,
  layers = layers,
  z_limits = c(-1, 1),
  step_breaks = seq(-1, 1, by = 0.2),
  palette = pal,
  title_size = 24,
  axis_text_size = 17,
  legend_text_size = 16
)

p2 <- plot_europe_robinson_tiles(
  data_mat = res_mat,
  sig_mask = res_sig,
  title_str = "b) Residual Correlation",
  lon = lon,
  lat = lat,
  layers = layers,
  z_limits = c(-1, 1),
  step_breaks = seq(-1, 1, by = 0.2),
  palette = pal,
  title_size = 24,
  axis_text_size = 17,
  legend_text_size = 16
)

p3 <- plot_europe_robinson_tiles(
  data_mat = diff_mat,
  sig_mask = diff_sig,
  title_str = "c) INIT - HIST+245 Correlation Difference",
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

save_single_panel(
  p1,
  file.path(out_dir, "Fig2a_INIT_vs_OBS.pdf")
)

save_single_panel(
  p2,
  file.path(out_dir, "Fig2b_residual_correlation.pdf")
)

save_single_panel(
  p3,
  file.path(out_dir, "Fig2c_difference_INIT_minus_HIST245.pdf")
)

save_three_panel_figure(
  p1, p2, p3,
  file.path(out_dir, "Figure2_deterministic_skill_raw_corrected.pdf")
)

cat("\nDONE\n")
cat("Input :", input_file, "\n")
cat("Years :", paste(range(x$years), collapse = "-"), "\n")
cat("Output:", out_dir, "\n")
