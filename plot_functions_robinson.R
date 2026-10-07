# Common plotting functions for the decadal FWI paper
# Keeps the Robinson-projection style used in the submitted manuscript.
#
# Required packages:
#   sf, ggplot2, reshape2, scales, rnaturalearth, grid
#
# This file only contains plotting utilities. It does not calculate skill.

suppressPackageStartupMessages({
  library(sf)
  library(ggplot2)
  library(reshape2)
  library(scales)
  library(rnaturalearth)
  library(grid)
})

sf::sf_use_s2(FALSE)

# -------------------------------------------------------------------------
# Base geographic layers
# -------------------------------------------------------------------------

load_robinson_base_layers <- function() {

  ocean <- rnaturalearth::ne_download(
    scale = 50,
    type = "ocean",
    category = "physical",
    returnclass = "sf"
  )

  land <- rnaturalearth::ne_download(
    scale = 50,
    type = "land",
    category = "physical",
    returnclass = "sf"
  )

  countries <- rnaturalearth::ne_countries(
    scale = 50,
    returnclass = "sf"
  )

  ocean     <- sf::st_make_valid(sf::st_transform(ocean, 4326))
  land      <- sf::st_make_valid(sf::st_transform(land, 4326))
  countries <- sf::st_make_valid(sf::st_transform(countries, 4326))

  land_union <- sf::st_make_valid(sf::st_union(land))

  list(
    ocean = ocean,
    land = land,
    countries = countries,
    land_union = land_union
  )
}


# -------------------------------------------------------------------------
# Generic Robinson tile map
# -------------------------------------------------------------------------

plot_europe_robinson_tiles <- function(
    data_mat,
    sig_mask = NULL,
    title_str = "",
    lon,
    lat,
    layers,
    z_limits,
    step_breaks,
    palette,
    xlim_deg = c(-15, 45),
    ylim_deg = c(30, 72),
    lon_breaks = seq(-15, 45, 10),
    lat_breaks = seq(30, 72, 10),
    tile_alpha = 1,
    ocean_fill = "lightgrey",
    graticule_color = "black",
    graticule_lwd = 0.4,
    graticule_lty = "solid",
    title_size = 22,
    axis_text_size = 15,
    legend_title_size = 14,
    legend_text_size = 13,
    sig_shape = 16,
    sig_size = 1.0,
    sig_alpha = 0.85
) {

  data_mat <- as.matrix(data_mat)

  stopifnot(
    length(lon) == nrow(data_mat),
    length(lat) == ncol(data_mat)
  )

  if (!is.null(sig_mask)) {
    sig_mask <- as.matrix(sig_mask)
    stopifnot(all(dim(sig_mask) == dim(data_mat)))
  }

  ocean      <- layers$ocean
  land       <- layers$land
  countries  <- layers$countries
  land_union <- layers$land_union

  # Wrap longitude to [-180, 180)
  wrap_lon <- function(x) ((x + 180) %% 360) - 180

  lon <- wrap_lon(as.numeric(lon))
  lat <- as.numeric(lat)

  ord <- order(lon)
  lon <- lon[ord]
  data_mat <- data_mat[ord, , drop = FALSE]

  if (!is.null(sig_mask)) {
    sig_mask <- sig_mask[ord, , drop = FALSE]
  }

  lon_u <- sort(unique(lon))
  lat_u <- sort(unique(lat))

  dlon <- if (length(lon_u) > 1) median(diff(lon_u), na.rm = TRUE) else 1
  dlat <- if (length(lat_u) > 1) median(diff(lat_u), na.rm = TRUE) else 1

  rownames(data_mat) <- lon
  colnames(data_mat) <- lat

  df <- reshape2::melt(
    data_mat,
    varnames = c("lon", "lat"),
    value.name = "z"
  )

  df$lon <- as.numeric(as.character(df$lon))
  df$lat <- as.numeric(as.character(df$lat))
  df <- df[is.finite(df$z), , drop = FALSE]

  # Build grid polygons around the model grid centres
  bb <- sf::st_bbox(
    c(
      xmin = min(lon_u) - dlon / 2,
      xmax = max(lon_u) + dlon / 2,
      ymin = min(lat_u) - dlat / 2,
      ymax = max(lat_u) + dlat / 2
    ),
    crs = sf::st_crs(4326)
  )

  grid_polys <- sf::st_make_grid(
    sf::st_as_sfc(bb),
    cellsize = c(dlon, dlat),
    what = "polygons",
    square = TRUE
  )

  grid_sf <- sf::st_sf(geometry = grid_polys)
  grid_sf <- sf::st_make_valid(sf::st_set_crs(grid_sf, 4326))

  cent <- sf::st_coordinates(sf::st_centroid(grid_sf))
  grid_sf$lon_c <- cent[, 1]
  grid_sf$lat_c <- cent[, 2]

  grid_sf <- grid_sf[
    grid_sf$lon_c >= xlim_deg[1] &
      grid_sf$lon_c <= xlim_deg[2] &
      grid_sf$lat_c >= ylim_deg[1] &
      grid_sf$lat_c <= ylim_deg[2],
    ,
    drop = FALSE
  ]

  grid_sf$lon_key <- lon_u[
    pmax(
      1,
      pmin(
        length(lon_u),
        round((grid_sf$lon_c - min(lon_u)) / dlon) + 1
      )
    )
  ]

  grid_sf$lat_key <- lat_u[
    pmax(
      1,
      pmin(
        length(lat_u),
        round((grid_sf$lat_c - min(lat_u)) / dlat) + 1
      )
    )
  ]

  df$key <- paste(df$lon, df$lat, sep = "_")
  grid_sf$key <- paste(grid_sf$lon_key, grid_sf$lat_key, sep = "_")

  grid_sf$z <- df$z[match(grid_sf$key, df$key)]
  grid_sf <- grid_sf[is.finite(grid_sf$z), , drop = FALSE]

  # Keep grid cells whose centroid falls on land
  cent_sf <- sf::st_as_sf(
    data.frame(
      lon = grid_sf$lon_c,
      lat = grid_sf$lat_c
    ),
    coords = c("lon", "lat"),
    crs = 4326,
    remove = FALSE
  )

  idx_land <- lengths(sf::st_within(cent_sf, land_union)) > 0
  grid_land <- grid_sf[idx_land, , drop = FALSE]

  cat("\n", title_str, "\n", sep = "")
  cat("Tiles over land:", nrow(grid_land), "\n")

  # Significant points, also restricted to land
  sig_pts <- NULL

  if (!is.null(sig_mask)) {

    rownames(sig_mask) <- lon
    colnames(sig_mask) <- lat

    ss <- reshape2::melt(
      sig_mask,
      varnames = c("lon", "lat"),
      value.name = "sig"
    )

    ss$lon <- as.numeric(as.character(ss$lon))
    ss$lat <- as.numeric(as.character(ss$lat))
    ss <- ss[!is.na(ss$sig) & ss$sig == TRUE, , drop = FALSE]

    if (nrow(ss) > 0) {

      sig_pts <- sf::st_as_sf(
        ss,
        coords = c("lon", "lat"),
        crs = 4326,
        remove = FALSE
      )

      idx_sig_land <- lengths(sf::st_within(sig_pts, land_union)) > 0
      sig_pts <- sig_pts[idx_sig_land, , drop = FALSE]

      if (nrow(sig_pts) == 0) sig_pts <- NULL
    }
  }

  rob_crs <- sf::st_crs(
    "+proj=robin +lon_0=0 +datum=WGS84 +units=m +no_defs"
  )

  bb_crop <- sf::st_bbox(
    c(
      xmin = xlim_deg[1],
      xmax = xlim_deg[2],
      ymin = ylim_deg[1],
      ymax = ylim_deg[2]
    ),
    crs = sf::st_crs(4326)
  )

  grat <- sf::st_graticule(
    lat = lat_breaks,
    lon = lon_breaks,
    bbox = bb_crop
  )

  step_breaks <- step_breaks[
    step_breaks >= z_limits[1] &
      step_breaks <= z_limits[2]
  ]

  step_vals <- scales::rescale(
    step_breaks,
    from = z_limits
  )

  p <- ggplot2::ggplot() +
    ggplot2::geom_sf(
      data = ocean,
      fill = ocean_fill,
      color = NA
    ) +
    ggplot2::geom_sf(
      data = grid_land,
      ggplot2::aes(fill = z),
      color = NA,
      alpha = tile_alpha
    ) +
    ggplot2::geom_sf(
      data = land,
      fill = NA,
      color = "grey35",
      linewidth = 0.20
    ) +
    ggplot2::geom_sf(
      data = countries,
      fill = NA,
      color = "grey20",
      linewidth = 0.25
    ) +
    ggplot2::geom_sf(
      data = grat,
      color = graticule_color,
      linewidth = graticule_lwd,
      linetype = graticule_lty
    ) +
    {
      if (!is.null(sig_pts)) {
        ggplot2::geom_sf(
          data = sig_pts,
          shape = sig_shape,
          size = sig_size,
          color = "black",
          alpha = sig_alpha
        )
      } else {
        NULL
      }
    } +
    ggplot2::scale_fill_stepsn(
      colors = palette,
      values = step_vals,
      breaks = step_breaks,
      limits = z_limits,
      oob = scales::squish,
      na.value = "transparent",
      name = "",
      guide = ggplot2::guide_colorsteps(
        title.position = "top",
        barwidth = grid::unit(0.90, "npc"),
        barheight = grid::unit(0.8, "cm"),
        show.limits = TRUE,
        ticks = TRUE,
        even.steps = TRUE
      )
    ) +
    ggplot2::coord_sf(
      crs = rob_crs,
      default_crs = sf::st_crs(4326),
      xlim = xlim_deg,
      ylim = ylim_deg,
      expand = FALSE
    ) +
    ggplot2::scale_x_continuous(
      breaks = lon_breaks,
      labels = function(x) paste0(x, "°")
    ) +
    ggplot2::scale_y_continuous(
      breaks = lat_breaks,
      labels = function(y) paste0(y, "°")
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(
        hjust = 0.5,
        face = "bold",
        size = title_size
      ),
      legend.position = "bottom",
      legend.direction = "horizontal",
      legend.box = "horizontal",
      legend.box.just = "center",
      legend.key.width = grid::unit(1, "npc"),
      legend.margin = ggplot2::margin(t = 2, r = 10, b = 2, l = 10),
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      panel.background = ggplot2::element_rect(fill = NA, color = NA),
      plot.background = ggplot2::element_rect(fill = "white", color = NA),
      axis.title = ggplot2::element_blank(),
      axis.text = ggplot2::element_text(
        color = "grey20",
        size = axis_text_size
      ),
      axis.ticks = ggplot2::element_line(color = "grey30"),
      legend.title = ggplot2::element_text(size = legend_title_size),
      legend.text = ggplot2::element_text(size = legend_text_size),
      plot.margin = ggplot2::margin(t = 8, r = 8, b = 8, l = 8)
    ) +
    ggplot2::labs(title = title_str)

  p
}


# -------------------------------------------------------------------------
# Palettes from the manuscript
# -------------------------------------------------------------------------

correlation_palette <- function() {
  # Negative = blue; positive = red, as in the submitted Figure 2
  rev(c(
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
  ))
}

rpss_palette <- function() {
  # Negative = purple; positive = green, as in the submitted RPSS figure
  c(
    "#8e0152",
    "#c51b7d",
    "#de77ae",
    "#f1b6da",
    "#fde0ef",
    "#e6f5d0",
    "#b8e186",
    "#7fbc41",
    "#4d9221",
    "#276419"
  )
}


# -------------------------------------------------------------------------
# Legend extraction + multipanel output
# -------------------------------------------------------------------------

extract_bottom_legend <- function(p) {

  g <- ggplot2::ggplotGrob(
    p + ggplot2::theme(legend.position = "bottom")
  )

  nm <- vapply(g$grobs, function(x) x$name, character(1))
  ii <- grep("^guide-box", nm)

  if (length(ii) == 0) {
    stop("Could not find legend in ggplot object.")
  }

  g$grobs[[ii[1]]]
}


save_three_panel_figure <- function(
    p1, p2, p3, filename,
    width = 13,
    height = 13
) {

  dir.create(dirname(filename), recursive = TRUE, showWarnings = FALSE)

  legend <- extract_bottom_legend(p1)

  p1n <- p1 + ggplot2::theme(legend.position = "none")
  p2n <- p2 + ggplot2::theme(legend.position = "none")
  p3n <- p3 + ggplot2::theme(legend.position = "none")

  grDevices::pdf(
    filename,
    width = width,
    height = height,
    onefile = TRUE
  )

  grid::grid.newpage()

  lay <- grid::grid.layout(
    nrow = 3,
    ncol = 2,
    heights = grid::unit(c(1, 1, 0.16), "null")
  )

  vp <- grid::viewport(layout = lay)
  grid::pushViewport(vp)

  print(
    p1n,
    vp = grid::viewport(
      layout.pos.row = 1,
      layout.pos.col = 1:2
    )
  )

  print(
    p2n,
    vp = grid::viewport(
      layout.pos.row = 2,
      layout.pos.col = 1
    )
  )

  print(
    p3n,
    vp = grid::viewport(
      layout.pos.row = 2,
      layout.pos.col = 2
    )
  )

  grid::grid.draw(
    grid::editGrob(
      legend,
      vp = grid::viewport(
        layout.pos.row = 3,
        layout.pos.col = 1:2
      )
    )
  )

  grid::popViewport()
  grDevices::dev.off()

  invisible(filename)
}


save_two_panel_figure <- function(
    p1, p2, filename,
    width = 13,
    height = 7.6
) {

  dir.create(dirname(filename), recursive = TRUE, showWarnings = FALSE)

  legend <- extract_bottom_legend(p1)

  p1n <- p1 + ggplot2::theme(legend.position = "none")
  p2n <- p2 + ggplot2::theme(legend.position = "none")

  grDevices::pdf(
    filename,
    width = width,
    height = height,
    onefile = TRUE
  )

  grid::grid.newpage()

  lay <- grid::grid.layout(
    nrow = 2,
    ncol = 2,
    heights = grid::unit(c(1, 0.18), "null")
  )

  vp <- grid::viewport(layout = lay)
  grid::pushViewport(vp)

  print(
    p1n,
    vp = grid::viewport(
      layout.pos.row = 1,
      layout.pos.col = 1
    )
  )

  print(
    p2n,
    vp = grid::viewport(
      layout.pos.row = 1,
      layout.pos.col = 2
    )
  )

  grid::grid.draw(
    grid::editGrob(
      legend,
      vp = grid::viewport(
        layout.pos.row = 2,
        layout.pos.col = 1:2
      )
    )
  )

  grid::popViewport()
  grDevices::dev.off()

  invisible(filename)
}


save_four_panel_figure <- function(
    p1, p2, p3, p4, filename,
    width = 13,
    height = 12
) {

  dir.create(dirname(filename), recursive = TRUE, showWarnings = FALSE)

  legend <- extract_bottom_legend(p1)

  pp <- lapply(
    list(p1, p2, p3, p4),
    function(x) x + ggplot2::theme(legend.position = "none")
  )

  grDevices::pdf(
    filename,
    width = width,
    height = height,
    onefile = TRUE
  )

  grid::grid.newpage()

  lay <- grid::grid.layout(
    nrow = 3,
    ncol = 2,
    heights = grid::unit(c(1, 1, 0.16), "null")
  )

  vp <- grid::viewport(layout = lay)
  grid::pushViewport(vp)

  print(pp[[1]], vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 1))
  print(pp[[2]], vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 2))
  print(pp[[3]], vp = grid::viewport(layout.pos.row = 2, layout.pos.col = 1))
  print(pp[[4]], vp = grid::viewport(layout.pos.row = 2, layout.pos.col = 2))

  grid::grid.draw(
    grid::editGrob(
      legend,
      vp = grid::viewport(
        layout.pos.row = 3,
        layout.pos.col = 1:2
      )
    )
  )

  grid::popViewport()
  grDevices::dev.off()

  invisible(filename)
}


save_single_panel <- function(
    p, filename,
    width = 7.6,
    height = 9.4
) {

  dir.create(dirname(filename), recursive = TRUE, showWarnings = FALSE)

  ggplot2::ggsave(
    filename,
    p,
    width = width,
    height = height,
    dpi = 300
  )

  invisible(filename)
}
