rm(list = ls())
graphics.off()

# Supplementary figure — deterministic prediction skill of FWI inputs
# Four panels: temperature, relative humidity, precipitation, wind speed.
#
# This uses INIT-vs-ERA5 correlation skill for the same FY1-FY5,
# 1963-2022 verification windows used for FWI95d.

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
  "verification_drivers_new",
  "results",
  "driver_skill_FY1-FY5_corrected.rds"
)

out_dir <- file.path(
  PROJECT_ROOT,
  "plots",
  "revision",
  "Supplementary_driver_skill_corrected"
)

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

x <- readRDS(input_file)

vars <- c("tas", "hurs", "pr", "sfcWind")

stopifnot(all(vars %in% names(x)))

lon <- as.numeric(x$tas$lon)
lat <- as.numeric(x$tas$lat)

for (v in vars) {
  stopifnot(
    identical(as.numeric(x[[v]]$lon), lon),
    identical(as.numeric(x[[v]]$lat), lat)
  )
}

layers <- load_robinson_base_layers()

pal <- correlation_palette()

make_driver_panel <- function(var, title) {

  mat <- as.matrix(x[[var]]$r_init_land)
  sig <- as.matrix(x[[var]]$sig_init_land)

  stopifnot(
    all(dim(mat) == c(length(lon), length(lat))),
    all(dim(sig) == dim(mat))
  )

  plot_europe_robinson_tiles(
    data_mat = mat,
    sig_mask = sig,
    title_str = title,
    lon = lon,
    lat = lat,
    layers = layers,
    z_limits = c(-1, 1),
    step_breaks = seq(-1, 1, by = 0.2),
    palette = pal,
    title_size = 23,
    axis_text_size = 16,
    legend_text_size = 15
  )
}

p1 <- make_driver_panel(
  "tas",
  "a) Temperature"
)

p2 <- make_driver_panel(
  "hurs",
  "b) Relative Humidity"
)

p3 <- make_driver_panel(
  "pr",
  "c) Precipitation"
)

p4 <- make_driver_panel(
  "sfcWind",
  "d) Wind Speed"
)

save_single_panel(
  p1,
  file.path(out_dir, "FigS_driver_tas.pdf")
)

save_single_panel(
  p2,
  file.path(out_dir, "FigS_driver_hurs.pdf")
)

save_single_panel(
  p3,
  file.path(out_dir, "FigS_driver_pr.pdf")
)

save_single_panel(
  p4,
  file.path(out_dir, "FigS_driver_sfcWind.pdf")
)

save_four_panel_figure(
  p1, p2, p3, p4,
  file.path(out_dir, "FigureS_driver_skill_INIT_vs_ERA5_corrected.pdf")
)

cat("\nDONE\n")
cat("Input :", input_file, "\n")
cat("Output:", out_dir, "\n")
