###############################################################################
#  SCRIPT COMPLETO DESDE CERO (4 MAPAS) — ROBINSON + TILES + PALETA ESCALONADA
#  - OBS / INIT / NO-INIT (runmean5)
#  - CorrEno INIT vs OBS
#  - CorrEno HIST+245 vs OBS
#  - ResidualCorr
#  - DiffCorr
#  - Plot final: graticule, puntos significativos, leyenda “a todo el ancho”
#  - El grid de correlación SOLO se pinta en celdas cuyo centroide cae en tierra
#  - La capa de land NO tapa la correlación (fill = NA)
#  - Diagnóstico completo FULL vs BBOX
#  - Resumen final con porcentajes listos para el abstract
###############################################################################

rm(list = ls())
gc()
while (dev.cur() > 1) dev.off()
try(closeAllConnections(), silent = TRUE)

#—————————————————————————————————————————————
# 0) Paquetes
#—————————————————————————————————————————————
suppressPackageStartupMessages({
  library(ncdf4)
  library(abind)
  library(s2dv)
  library(fields)
  library(sf)
  library(ggplot2)
  library(reshape2)
  library(RColorBrewer)
  library(scales)
  library(rnaturalearth)
  library(grid)
})
# install.packages("s2dv")
sf::sf_use_s2(FALSE)

#—————————————————————————————————————————————
# 0.1) Rutas
#—————————————————————————————————————————————
out_dir  <- "/mnt/disco6tb/Dropbox/decadal/plots/"
obs_dir  <- "/mnt/disco6tb/Dropbox/decadal/data/runmean5/"
source("CorrEno.R")

if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

#—————————————————————————————————————————————
# 0.2) Capas base
#—————————————————————————————————————————————
ocean     <- rnaturalearth::ne_download(scale = 50, type = "ocean", category = "physical", returnclass = "sf")
land      <- rnaturalearth::ne_download(scale = 50, type = "land",  category = "physical", returnclass = "sf")
countries <- rnaturalearth::ne_countries(scale = 50, returnclass = "sf")

ocean     <- st_transform(ocean, 4326)
land      <- st_transform(land, 4326)
countries <- st_transform(countries, 4326)

ocean     <- st_make_valid(ocean)
land      <- st_make_valid(land)
countries <- st_make_valid(countries)

land_union <- st_union(land)
land_union <- st_make_valid(land_union)

#—————————————————————————————————————————————
# 0.3) Función de plot
#—————————————————————————————————————————————
plot_europe_robinson_tiles <- function(data_mat, sig_mask,
                                       title_str, filename,
                                       lon, lat, out_dir,
                                       ocean, land, countries, land_union,
                                       xlim_deg = c(-15, 45),
                                       ylim_deg = c(30, 72),
                                       lon_breaks = seq(-15, 45, 10),
                                       lat_breaks = seq(30, 72, 10),
                                       r_limits = c(-1, 1),
                                       step_breaks = seq(-1, 1, by = 0.2),
                                       tile_alpha = 1,
                                       ocean_fill = "lightgrey",
                                       graticule_color = "black",
                                       graticule_lwd = 0.4,
                                       graticule_lty = "solid",
                                       title_size = 18,
                                       axis_text_size = 12,
                                       cbar_height = unit(1.0, "cm"),
                                       cbar_width  = unit(0.95, "npc"),
                                       legend_title_size = 14,
                                       legend_text_size  = 13,
                                       sig_shape = 16,
                                       sig_size  = 1.1,
                                       sig_alpha = 0.85,
                                       reverse_prgn = TRUE) {
  
  stopifnot(is.matrix(data_mat))
  stopifnot(length(lon) == nrow(data_mat))
  stopifnot(length(lat) == ncol(data_mat))
  
  wrap_lon <- function(x) ((x + 180) %% 360) - 180
  lon <- wrap_lon(lon)
  o <- order(lon)
  
  lon <- lon[o]
  data_mat <- data_mat[o, , drop = FALSE]
  if (!is.null(sig_mask)) {
    stopifnot(all(dim(sig_mask) == dim(data_mat)))
    sig_mask <- sig_mask[o, , drop = FALSE]
  }
  
  lon_u <- sort(unique(lon))
  lat_u <- sort(unique(lat))
  dlon <- if (length(lon_u) > 1) median(diff(lon_u), na.rm = TRUE) else 1
  dlat <- if (length(lat_u) > 1) median(diff(lat_u), na.rm = TRUE) else 1
  
  colnames(data_mat) <- lat
  rownames(data_mat) <- lon
  
  df <- reshape2::melt(data_mat, varnames = c("lon", "lat"), value.name = "r")
  df$lon <- as.numeric(as.character(df$lon))
  df$lat <- as.numeric(as.character(df$lat))
  df <- df[is.finite(df$r), , drop = FALSE]
  
  # bbox del grid
  bb <- sf::st_bbox(c(
    xmin = min(lon_u) - dlon/2, xmax = max(lon_u) + dlon/2,
    ymin = min(lat_u) - dlat/2, ymax = max(lat_u) + dlat/2
  ), crs = sf::st_crs(4326))
  
  # crear grid
  grid_polys <- sf::st_make_grid(
    sf::st_as_sfc(bb),
    cellsize = c(dlon, dlat),
    what = "polygons",
    square = TRUE
  )
  
  grid_sf <- sf::st_sf(geometry = grid_polys)
  grid_sf <- st_set_crs(grid_sf, 4326)
  grid_sf <- st_make_valid(grid_sf)
  
  # centroides
  cent <- sf::st_coordinates(sf::st_centroid(grid_sf))
  grid_sf$lon_c <- cent[, 1]
  grid_sf$lat_c <- cent[, 2]
  
  # recorte a bbox de la figura
  grid_sf <- grid_sf[
    grid_sf$lon_c >= xlim_deg[1] & grid_sf$lon_c <= xlim_deg[2] &
      grid_sf$lat_c >= ylim_deg[1] & grid_sf$lat_c <= ylim_deg[2], ,
    drop = FALSE
  ]
  
  # asignar valor a cada tile
  grid_sf$lon_key <- lon_u[pmax(1, pmin(length(lon_u), round((grid_sf$lon_c - min(lon_u))/dlon) + 1))]
  grid_sf$lat_key <- lat_u[pmax(1, pmin(length(lat_u), round((grid_sf$lat_c - min(lat_u))/dlat) + 1))]
  
  df$key <- paste(df$lon, df$lat, sep = "_")
  grid_sf$key <- paste(grid_sf$lon_key, grid_sf$lat_key, sep = "_")
  grid_sf$r <- df$r[match(grid_sf$key, df$key)]
  grid_sf <- grid_sf[is.finite(grid_sf$r), , drop = FALSE]
  
  # quedarse SOLO con celdas cuyo centroide cae en tierra
  cent_sf <- st_as_sf(
    data.frame(lon = grid_sf$lon_c, lat = grid_sf$lat_c),
    coords = c("lon", "lat"),
    crs = 4326,
    remove = FALSE
  )
  
  idx_land <- lengths(st_within(cent_sf, land_union)) > 0
  grid_land <- grid_sf[idx_land, , drop = FALSE]
  
  cat("\n", title_str, "\n")
  cat("Tiles bbox:", nrow(grid_sf), "\n")
  cat("Tiles sobre tierra (centroide):", nrow(grid_land), "\n")
  
  # significancia
  sig_pts <- NULL
  if (!is.null(sig_mask)) {
    colnames(sig_mask) <- lat
    rownames(sig_mask) <- lon
    
    s <- reshape2::melt(sig_mask, varnames = c("lon", "lat"), value.name = "sig")
    s$lon <- as.numeric(as.character(s$lon))
    s$lat <- as.numeric(as.character(s$lat))
    s <- s[!is.na(s$sig), , drop = FALSE]
    s <- subset(s, sig == TRUE)
    
    if (nrow(s) > 0) {
      sig_pts <- st_as_sf(s, coords = c("lon", "lat"), crs = 4326, remove = FALSE)
      idx_sig_land <- lengths(st_within(sig_pts, land_union)) > 0
      sig_pts <- sig_pts[idx_sig_land, , drop = FALSE]
      if (nrow(sig_pts) == 0) sig_pts <- NULL
    }
  }
  
  # proyección Robinson
  rob_crs <- st_crs("+proj=robin +lon_0=0 +datum=WGS84 +units=m +no_defs")
  
  # graticule
  bb_crop <- st_bbox(c(
    xmin = xlim_deg[1], xmax = xlim_deg[2],
    ymin = ylim_deg[1], ymax = ylim_deg[2]
  ), crs = st_crs(4326))
  grat <- st_graticule(lat = lat_breaks, lon = lon_breaks, bbox = bb_crop)
  
  # paleta
  pal <- c(
    "#a50026",
    "#d73027",
    "#f46d43",
    "#fdae61",
    "#fee090",
    "#ffffff",
    "#e0f3f8",
    "#abd9e9",
    "#74add1",
    "#4575b4",
    "#313695"
  )
  if (isTRUE(reverse_prgn)) pal <- rev(pal)
  
  step_breaks <- step_breaks[step_breaks >= r_limits[1] & step_breaks <= r_limits[2]]
  step_vals   <- scales::rescale(step_breaks, from = r_limits)
  
  p <- ggplot() +
    geom_sf(data = ocean, fill = ocean_fill, color = NA) +
    geom_sf(data = grid_land, aes(fill = r), color = NA, alpha = tile_alpha) +
    geom_sf(data = land, fill = NA, color = "grey35", linewidth = 0.20) +
    geom_sf(data = countries, fill = NA, color = "grey20", linewidth = 0.25) +
    geom_sf(data = grat, color = graticule_color, linewidth = graticule_lwd, linetype = graticule_lty) +
    { if (!is.null(sig_pts)) geom_sf(data = sig_pts, shape = sig_shape, size = sig_size, color = "black", alpha = sig_alpha) else NULL } +
    scale_fill_stepsn(
      colors = pal,
      values = step_vals,
      breaks = step_breaks,
      limits = r_limits,
      oob = scales::squish,
      na.value = "transparent",
      name = "",
      guide = guide_colorsteps(
        title.position = "top",
        barwidth  = cbar_width,
        barheight = cbar_height,
        show.limits = TRUE,
        ticks = TRUE,
        even.steps = TRUE
      )
    ) +
    coord_sf(
      crs = rob_crs,
      default_crs = st_crs(4326),
      xlim = xlim_deg, ylim = ylim_deg,
      expand = FALSE
    ) +
    scale_x_continuous(breaks = lon_breaks, labels = function(x) paste0(x, "°")) +
    scale_y_continuous(breaks = lat_breaks, labels = function(y) paste0(y, "°")) +
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = title_size),
      legend.position = "bottom",
      legend.direction = "horizontal",
      legend.box = "horizontal",
      legend.box.just = "center",
      legend.key.width = unit(1, "npc"),
      legend.margin = margin(t = 2, r = 10, b = 2, l = 10),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_rect(fill = NA, color = NA),
      plot.background  = element_rect(fill = "white", color = NA),
      axis.title = element_blank(),
      axis.text  = element_text(color = "grey20", size = axis_text_size),
      axis.ticks = element_line(color = "grey30"),
      legend.title = element_text(size = legend_title_size),
      legend.text  = element_text(size = legend_text_size),
      plot.margin = margin(t = 8, r = 8, b = 14, l = 8)
    ) +
    labs(title = title_str)
  
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  ggsave(file.path(out_dir, filename), p, width = 7.6, height = 9.4, dpi = 300)
  
  return(p)
}

#—————————————————————————————————————————————
# 1) Observaciones (runmean5)
#—————————————————————————————————————————————
obs_file <- file.path(obs_dir, "FWI_obs_FWI95d_1961-2024_fixed_runmean5.nc")

nc       <- nc_open(obs_file)
lon      <- ncvar_get(nc, "lon")
lat      <- ncvar_get(nc, "lat")
time_obs <- ncvar_get(nc, "time")
obs      <- ncvar_get(nc, "FWI")
nc_close(nc)

dimnames(obs) <- list(lon = lon, lat = lat, sdate = as.character(time_obs))
names(dimnames(obs)) <- c("lon", "lat", "sdate")
.dobs <- dim(obs); names(.dobs) <- c("lon", "lat", "sdate"); dim(obs) <- .dobs

#—————————————————————————————————————————————
# 2) Predicciones INIT (runmean5)
#—————————————————————————————————————————————
dcpp_files <- list.files(
  obs_dir,
  pattern = "FWI95d_reconstructed_.*_1961-2024_runmean5\\.nc$",
  full.names = TRUE
)

nmem_dcpp <- length(dcpp_files)
ntime     <- length(time_obs)
nlon      <- length(lon)
nlat      <- length(lat)

dcpp <- array(NA_real_, dim = c(nlon, nlat, ntime, nmem_dcpp))
member_names_dcpp <- character(nmem_dcpp)

for (i in seq_along(dcpp_files)) {
  nc  <- nc_open(dcpp_files[i])
  tmp <- ncvar_get(nc, "FWI")
  nc_close(nc)
  
  if (!all(dim(tmp) == c(nlon, nlat, ntime))) {
    stop("Un fichero INIT no tiene las mismas dimensiones que obs.")
  }
  
  dcpp[,,,i] <- tmp
  member_names_dcpp[i] <- sub(".*reconstructed_(r[0-9]+i1p1f1)_.*", "\\1", basename(dcpp_files[i]))
}

dimnames(dcpp) <- list(lon = lon, lat = lat, sdate = as.character(time_obs), member = member_names_dcpp)
names(dimnames(dcpp)) <- c("lon", "lat", "sdate", "member")
.dcpp <- dim(dcpp); names(.dcpp) <- c("lon", "lat", "sdate", "member"); dim(dcpp) <- .dcpp

#—————————————————————————————————————————————
# 3) Histórico NO-INIT (runmean5)
#—————————————————————————————————————————————
hist_files <- list.files(
  obs_dir,
  pattern = "FWI_CMCC-CM2-SR5_.*_1960_2024_FWI95d_runmean5\\.nc$",
  full.names = TRUE
)
stopifnot(length(hist_files) == 10)

hist <- array(NA_real_, dim = c(nlon, nlat, ntime, length(hist_files)))
member_names_hist <- character(length(hist_files))

for (i in seq_along(hist_files)) {
  nc    <- nc_open(hist_files[i])
  tmp_h <- ncvar_get(nc, "FWI")
  nc_close(nc)
  
  if (!all(dim(tmp_h) == c(nlon, nlat, ntime))) {
    stop("Un fichero HIST+245 no tiene las mismas dimensiones que obs/dcpp.")
  }
  
  hist[,,,i] <- tmp_h
  member_names_hist[i] <- sub("FWI_CMCC-CM2-SR5_(r[0-9]+i1p2f1)_.*", "\\1", basename(hist_files[i]))
}

# Si tus ficheros runmean5 ya están alineados con obs, usa time_obs
dimnames(hist) <- list(lon = lon, lat = lat, sdate = as.character(time_obs), member = member_names_hist)
names(dimnames(hist)) <- c("lon", "lat", "sdate", "member")
.ndh <- dim(hist); names(.ndh) <- c("lon", "lat", "sdate", "member"); dim(hist) <- .ndh

# Comprobación rápida
if (!identical(dim(obs)[1:3], dim(dcpp)[1:3])) stop("obs y dcpp no tienen mismas dimensiones espaciales/temporales.")
if (!identical(dim(obs)[1:3], dim(hist)[1:3])) stop("obs y hist no tienen mismas dimensiones espaciales/temporales.")

#—————————————————————————————————————————————
# 4) Correlación Pearson (INIT vs OBS)
#—————————————————————————————————————————————
corr <- CorrEno(
  exp        = dcpp,
  obs        = obs,
  time_dim   = "sdate",
  member_dim = 4,
  method     = "pearson",
  alpha      = 0.05,
  test.type  = "two-sided",
  ncores     = 8,
  pval       = TRUE
)

#—————————————————————————————————————————————
# 5) Correlación Pearson (HIST+245 vs OBS)
#—————————————————————————————————————————————
corr_hist <- CorrEno(
  exp        = hist,
  obs        = obs,
  time_dim   = "sdate",
  member_dim = 4,
  method     = "pearson",
  alpha      = 0.05,
  test.type  = "two-sided",
  ncores     = 8,
  pval       = TRUE
)

#—————————————————————————————————————————————
# 6) Residual correlation
#—————————————————————————————————————————————
resid_corr <- s2dv::ResidualCorr(
  exp        = dcpp,
  obs        = obs,
  ref        = hist,
  time_dim   = "sdate",
  memb_dim   = "member",
  method     = "pearson",
  handle.na  = "return.na",
  pval       = FALSE,
  sign       = TRUE,
  alpha      = 0.05,
  ncores     = 8
)

#—————————————————————————————————————————————
# 7) Diferencia de correlación
#—————————————————————————————————————————————
diff_corr <- s2dv::DiffCorr(
  exp       = dcpp,
  ref       = hist,
  obs       = obs,
  time_dim  = "sdate",
  memb_dim  = "member",
  method    = "pearson",
  handle.na = "return.na",
  pval      = FALSE,
  sign      = TRUE,
  alpha     = 0.05,
  ncores    = 8
)

#—————————————————————————————————————————————
# 8) Mapas finales
#—————————————————————————————————————————————
cat("Generando los 4 mapas finales...\n")

xlim_deg_common   <- c(-15, 45)
ylim_deg_common   <- c(30, 72)
lon_breaks_common <- seq(-15, 45, 10)
lat_breaks_common <- seq(30, 72, 10)

p1 <- plot_europe_robinson_tiles(
  data_mat  = corr$r,
  sig_mask  = corr$sign,
  title_str = "a) Correlation (INIT vs OBS)",
  filename  = "map_skill_init_robinson_steps_landonly.pdf",
  lon = lon, lat = lat, out_dir = out_dir,
  land = land, ocean = ocean, countries = countries, land_union = land_union,
  xlim_deg = xlim_deg_common,
  ylim_deg = ylim_deg_common,
  lon_breaks = lon_breaks_common,
  lat_breaks = lat_breaks_common,
  r_limits = c(-1, 1),
  step_breaks = seq(-1, 1, by = 0.2),
  title_size = 25,
  axis_text_size = 18,
  cbar_width  = unit(0.90, "npc"),
  cbar_height = unit(1.0, "cm"),
  legend_title_size = 14,
  legend_text_size  = 18,
  reverse_prgn = TRUE
)

p2 <- plot_europe_robinson_tiles(
  data_mat  = corr_hist$r,
  sig_mask  = corr_hist$sign,
  title_str = "b) Correlation (HIST+245 vs OBS)",
  filename  = "map_skill_hist245_robinson_steps_landonly.pdf",
  lon = lon, lat = lat, out_dir = out_dir,
  land = land, ocean = ocean, countries = countries, land_union = land_union,
  xlim_deg = xlim_deg_common,
  ylim_deg = ylim_deg_common,
  lon_breaks = lon_breaks_common,
  lat_breaks = lat_breaks_common,
  r_limits = c(-1, 1),
  step_breaks = seq(-1, 1, by = 0.2),
  title_size = 25,
  axis_text_size = 18,
  cbar_width  = unit(0.90, "npc"),
  cbar_height = unit(1.0, "cm"),
  legend_title_size = 14,
  legend_text_size  = 18,
  reverse_prgn = TRUE
)

p3 <- plot_europe_robinson_tiles(
  data_mat  = resid_corr$res.corr,
  sig_mask  = resid_corr$sign,
  title_str = "c) Residual Correlation",
  filename  = "map_residual_skill_robinson_steps_landonly.pdf",
  lon = lon, lat = lat, out_dir = out_dir,
  land = land, ocean = ocean, countries = countries, land_union = land_union,
  xlim_deg = xlim_deg_common,
  ylim_deg = ylim_deg_common,
  lon_breaks = lon_breaks_common,
  lat_breaks = lat_breaks_common,
  r_limits = c(-1, 1),
  step_breaks = seq(-1, 1, by = 0.2),
  title_size = 25,
  axis_text_size = 18,
  cbar_width  = unit(0.90, "npc"),
  cbar_height = unit(1.0, "cm"),
  legend_title_size = 14,
  legend_text_size  = 18,
  reverse_prgn = TRUE
)

p4 <- plot_europe_robinson_tiles(
  data_mat  = diff_corr$diff,
  sig_mask  = diff_corr$sign,
  title_str = "d) Difference in Correlation (INIT - HIST+245)",
  filename  = "map_added_value_robinson_steps_landonly.pdf",
  lon = lon, lat = lat, out_dir = out_dir,
  land = land, ocean = ocean, countries = countries, land_union = land_union,
  xlim_deg = xlim_deg_common,
  ylim_deg = ylim_deg_common,
  lon_breaks = lon_breaks_common,
  lat_breaks = lat_breaks_common,
  r_limits = c(-1, 1),
  step_breaks = seq(-1, 1, by = 0.2),
  title_size = 24,
  axis_text_size = 18,
  cbar_width  = unit(0.90, "npc"),
  cbar_height = unit(1.0, "cm"),
  legend_title_size = 14,
  legend_text_size  = 18,
  reverse_prgn = TRUE
)

cat("¡Listo! PDFs guardados en: ", out_dir, "\n")

###############################################################################
# DIAGNÓSTICO + ESTADÍSTICAS COMPLETAS (DOMINIO COMPLETO vs BBOX FIGURA)
###############################################################################

make_df <- function(mat, sig, lon, lat, name = "val") {
  stopifnot(is.matrix(mat))
  stopifnot(all(dim(mat) == c(length(lon), length(lat))))
  
  df <- reshape2::melt(mat, varnames = c("lon", "lat"), value.name = name)
  df$lon <- as.numeric(as.character(df$lon))
  df$lat <- as.numeric(as.character(df$lat))
  
  s <- reshape2::melt(sig, varnames = c("lon", "lat"), value.name = "sig")
  df$sig <- as.logical(s$sig)
  df$w <- cos(df$lat * pi / 180)
  df
}

crop_bbox <- function(df, xlim, ylim) {
  subset(df,
         lon >= xlim[1] & lon <= xlim[2] &
           lat >= ylim[1] & lat <= ylim[2])
}

summ_domain <- function(df, value_col = "val", sig_col = "sig") {
  n_total_rows <- nrow(df)
  ok <- is.finite(df[[value_col]])
  d  <- df[ok, , drop = FALSE]
  n  <- nrow(d)
  
  sig_avail <- !is.na(d[[sig_col]])
  n_sig_avail <- sum(sig_avail)
  
  pos <- d[[value_col]] > 0
  neg <- d[[value_col]] < 0
  zer <- d[[value_col]] == 0
  
  n_pos <- sum(pos)
  n_neg <- sum(neg)
  n_zer <- sum(zer)
  
  sig_true <- (d[[sig_col]] == TRUE)
  n_sig <- sum(sig_true, na.rm = TRUE)
  
  n_pos_sig <- sum(pos & sig_true, na.rm = TRUE)
  n_neg_sig <- sum(neg & sig_true, na.rm = TRUE)
  n_zer_sig <- sum(zer & sig_true, na.rm = TRUE)
  
  mean_all <- if (n > 0) mean(d[[value_col]]) else NA_real_
  mean_sig <- if (n_sig > 0) mean(d[[value_col]][sig_true]) else NA_real_
  mean_pos_sig <- if (n_pos_sig > 0) mean(d[[value_col]][pos & sig_true]) else NA_real_
  mean_neg_sig <- if (n_neg_sig > 0) mean(d[[value_col]][neg & sig_true]) else NA_real_
  
  pct_pos <- if (n > 0) 100 * n_pos / n else NA_real_
  pct_neg <- if (n > 0) 100 * n_neg / n else NA_real_
  pct_sig_total <- if (n > 0) 100 * n_sig / n else NA_real_
  pct_pos_sig_total <- if (n > 0) 100 * n_pos_sig / n else NA_real_
  pct_neg_sig_total <- if (n > 0) 100 * n_neg_sig / n else NA_real_
  pct_sig_in_pos <- if (n_pos > 0) 100 * n_pos_sig / n_pos else NA_real_
  pct_sig_in_neg <- if (n_neg > 0) 100 * n_neg_sig / n_neg else NA_real_
  
  sanity_ok <- (n_sig == (n_pos_sig + n_neg_sig + n_zer_sig))
  
  list(
    n_total_rows = n_total_rows,
    n_valid = n,
    n_sig_avail = n_sig_avail,
    n_sig = n_sig,
    mean_all = mean_all,
    mean_sig = mean_sig,
    pct_pos = pct_pos,
    pct_neg = pct_neg,
    pct_sig_total = pct_sig_total,
    pct_pos_sig_total = pct_pos_sig_total,
    pct_neg_sig_total = pct_neg_sig_total,
    pct_sig_in_pos = pct_sig_in_pos,
    pct_sig_in_neg = pct_sig_in_neg,
    mean_pos_sig = mean_pos_sig,
    mean_neg_sig = mean_neg_sig,
    sanity = c(
      n_sig = n_sig,
      n_pos_sig = n_pos_sig,
      n_neg_sig = n_neg_sig,
      n_zer_sig = n_zer_sig,
      sum_parts = n_pos_sig + n_neg_sig + n_zer_sig,
      sanity_ok = sanity_ok
    )
  )
}

print_summary <- function(title, res) {
  cat("\n============================\n")
  cat(title, "\n")
  cat("============================\n")
  cat(sprintf("Total filas (antes NA): %d\n", res$n_total_rows))
  cat(sprintf("Celdas válidas (finito): %d\n", res$n_valid))
  cat(sprintf("Celdas significativas: %d (%.1f%% de válidas)\n", res$n_sig, res$pct_sig_total))
  cat(sprintf("Media (todas): %.3f | Media (solo sig): %.3f\n", res$mean_all, res$mean_sig))
  cat(sprintf("%%pos = %.1f; %%neg = %.1f\n", res$pct_pos, res$pct_neg))
  cat(sprintf("%%pos&sig = %.1f; %%neg&sig = %.1f  (sobre válidas)\n",
              res$pct_pos_sig_total, res$pct_neg_sig_total))
  cat(sprintf("sig dentro de pos = %.1f; sig dentro de neg = %.1f\n",
              res$pct_sig_in_pos, res$pct_sig_in_neg))
  cat(sprintf("Media (pos&sig) = %.3f; Media (neg&sig) = %.3f\n",
              res$mean_pos_sig, res$mean_neg_sig))
  cat("Sanity check:\n")
  print(res$sanity)
}

check_matrix <- function(mat, name) {
  cat("\n---", name, "---\n")
  cat("Dim:", paste(dim(mat), collapse = " x "), "\n")
  cat("Total celdas:", prod(dim(mat)), "\n")
  cat("Finitas:", sum(is.finite(mat)), "\n")
  cat("No finitas (NA/Inf):", sum(!is.finite(mat)), "\n")
}

check_matrix(corr$r, "INIT (corr$r)")
check_matrix(corr_hist$r, "HIST+245 (corr_hist$r)")
check_matrix(resid_corr$res.corr, "Residual (resid_corr$res.corr)")
check_matrix(diff_corr$diff, "Diff (diff_corr$diff)")

cat("\nDominio teórico (lon x lat) =", length(lon) * length(lat), "celdas\n")

df_init_full  <- make_df(corr$r,              corr$sign,       lon, lat, name = "val")
df_hist_full  <- make_df(corr_hist$r,         corr_hist$sign,  lon, lat, name = "val")
df_resid_full <- make_df(resid_corr$res.corr, resid_corr$sign, lon, lat, name = "val")
df_diff_full  <- make_df(diff_corr$diff,      diff_corr$sign,  lon, lat, name = "val")

df_init_bbox  <- crop_bbox(df_init_full,  xlim_deg_common, ylim_deg_common)
df_hist_bbox  <- crop_bbox(df_hist_full,  xlim_deg_common, ylim_deg_common)
df_resid_bbox <- crop_bbox(df_resid_full, xlim_deg_common, ylim_deg_common)
df_diff_bbox  <- crop_bbox(df_diff_full,  xlim_deg_common, ylim_deg_common)

r_init_full  <- summ_domain(df_init_full)
r_hist_full  <- summ_domain(df_hist_full)
r_resid_full <- summ_domain(df_resid_full)
r_diff_full  <- summ_domain(df_diff_full)

r_init_bbox  <- summ_domain(df_init_bbox)
r_hist_bbox  <- summ_domain(df_hist_bbox)
r_resid_bbox <- summ_domain(df_resid_bbox)
r_diff_bbox  <- summ_domain(df_diff_bbox)

print_summary("INIT vs OBS — DOMINIO COMPLETO",     r_init_full)
print_summary("INIT vs OBS — BBOX FIGURA",          r_init_bbox)

print_summary("HIST+245 vs OBS — DOMINIO COMPLETO", r_hist_full)
print_summary("HIST+245 vs OBS — BBOX FIGURA",      r_hist_bbox)

print_summary("ResidualCorr — DOMINIO COMPLETO",    r_resid_full)
print_summary("ResidualCorr — BBOX FIGURA",         r_resid_bbox)

print_summary("DiffCorr — DOMINIO COMPLETO",        r_diff_full)
print_summary("DiffCorr — BBOX FIGURA",             r_diff_bbox)

to_row <- function(map_name, dom_name, r) {
  data.frame(
    map = map_name,
    domain = dom_name,
    n_total_rows = r$n_total_rows,
    n_valid = r$n_valid,
    n_sig = r$n_sig,
    pct_sig_total = r$pct_sig_total,
    mean_all = r$mean_all,
    mean_sig = r$mean_sig,
    pct_pos = r$pct_pos,
    pct_neg = r$pct_neg,
    pct_pos_sig_total = r$pct_pos_sig_total,
    pct_neg_sig_total = r$pct_neg_sig_total,
    pct_sig_in_pos = r$pct_sig_in_pos,
    pct_sig_in_neg = r$pct_sig_in_neg,
    mean_pos_sig = r$mean_pos_sig,
    mean_neg_sig = r$mean_neg_sig,
    sanity_ok = as.logical(r$sanity["sanity_ok"]),
    stringsAsFactors = FALSE
  )
}

out_compare <- rbind(
  to_row("INIT_vs_OBS",    "FULL", r_init_full),
  to_row("INIT_vs_OBS",    "BBOX", r_init_bbox),
  to_row("HIST245_vs_OBS", "FULL", r_hist_full),
  to_row("HIST245_vs_OBS", "BBOX", r_hist_bbox),
  to_row("ResidualCorr",   "FULL", r_resid_full),
  to_row("ResidualCorr",   "BBOX", r_resid_bbox),
  to_row("DiffCorr",       "FULL", r_diff_full),
  to_row("DiffCorr",       "BBOX", r_diff_bbox)
)

print(out_compare)

write.csv(out_compare, file.path(out_dir, "debug_domain_full_vs_bbox.csv"), row.names = FALSE)

###############################################################################
# RESUMEN FINAL PARA ABSTRACT / RESPUESTA A COAUTORES
###############################################################################

init_pos_sig  <- r_init_full$pct_pos_sig_total
init_neg_sig  <- r_init_full$pct_neg_sig_total
hist_pos_sig  <- r_hist_full$pct_pos_sig_total
hist_neg_sig  <- r_hist_full$pct_neg_sig_total
resid_pos_sig <- r_resid_full$pct_pos_sig_total
resid_neg_sig <- r_resid_full$pct_neg_sig_total
gain_sig      <- r_diff_full$pct_pos_sig_total
loss_sig      <- r_diff_full$pct_neg_sig_total

init_pos      <- r_init_full$pct_pos
init_neg      <- r_init_full$pct_neg
hist_pos      <- r_hist_full$pct_pos
hist_neg      <- r_hist_full$pct_neg
resid_pos     <- r_resid_full$pct_pos
resid_neg     <- r_resid_full$pct_neg
diff_pos      <- r_diff_full$pct_pos
diff_neg      <- r_diff_full$pct_neg

cat("\n############################################################\n")
cat("RESUMEN FINAL — PORCENTAJES LISTOS PARA MANUSCRITO\n")
cat("############################################################\n")

cat(sprintf("\nINIT vs OBS:\n"))
cat(sprintf("  Positive correlations: %.1f%%\n", init_pos))
cat(sprintf("  Negative correlations: %.1f%%\n", init_neg))
cat(sprintf("  Positive & significant: %.1f%%\n", init_pos_sig))
cat(sprintf("  Negative & significant: %.1f%%\n", init_neg_sig))

cat(sprintf("\nHIST+245 vs OBS:\n"))
cat(sprintf("  Positive correlations: %.1f%%\n", hist_pos))
cat(sprintf("  Negative correlations: %.1f%%\n", hist_neg))
cat(sprintf("  Positive & significant: %.1f%%\n", hist_pos_sig))
cat(sprintf("  Negative & significant: %.1f%%\n", hist_neg_sig))

cat(sprintf("\nResidualCorr:\n"))
cat(sprintf("  Positive correlations: %.1f%%\n", resid_pos))
cat(sprintf("  Negative correlations: %.1f%%\n", resid_neg))
cat(sprintf("  Positive & significant: %.1f%%\n", resid_pos_sig))
cat(sprintf("  Negative & significant: %.1f%%\n", resid_neg_sig))

cat(sprintf("\nDiffCorr (INIT - HIST+245):\n"))
cat(sprintf("  Positive differences: %.1f%%\n", diff_pos))
cat(sprintf("  Negative differences: %.1f%%\n", diff_neg))
cat(sprintf("  Positive & significant added value: %.1f%%\n", gain_sig))
cat(sprintf("  Negative & significant differences: %.1f%%\n", loss_sig))

cat("\n------------------------------------------------------------\n")
cat("Texto orientativo para abstract:\n\n")
cat(sprintf(
  "We find moderate yet spatially robust deterministic skill, with %.1f%% of the domain showing statistically significant positive correlations, compared with %.1f%% in the non-initialised simulations. Initialization adds skill beyond the externally forced signal, with %.1f%% of the domain exhibiting significant positive residual correlations. Direct comparison with the non-initialised simulations shows a significant increase in skill over %.1f%% of the domain, compared with %.1f%% showing a significant decrease.\n",
  init_pos_sig, hist_pos_sig, resid_pos_sig, gain_sig, loss_sig
))
cat("------------------------------------------------------------\n")

###############################################################################
# CHEQUEO EXTRA
###############################################################################
check_domain <- function(mat, name) {
  cat("\n---", name, "---\n")
  cat("Dim:", paste(dim(mat), collapse = " x "), "\n")
  cat("Total cells:", prod(dim(mat)), "\n")
  cat("Finite cells:", sum(is.finite(mat)), "\n")
  cat("NA cells:", sum(!is.finite(mat)), "\n")
}

check_domain(corr$r, "INIT")
check_domain(corr_hist$r, "HIST+245")
check_domain(resid_corr$res.corr, "Residual")
check_domain(diff_corr$diff, "Diff")

cat("\nFIN DEL SCRIPT\n")

