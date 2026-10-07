rm(list = ls())
graphics.off()

library(ncdf4)
library(s2dv)
library(multiApply)
library(fields)
library(maps)

# ============================================================
# Paths
# ============================================================

data_root <- "/Users/marco/Library/CloudStorage/Dropbox/estcena/scripts/ONFIRE/decadal/data"
driver_dir <- file.path(data_root, "verification_drivers_new")
out_dir <- file.path(driver_dir, "results")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

vars <- c("tas", "hurs", "pr", "sfcWind")

years <- 1963:2022
ntime <- length(years)

# ============================================================
# CorrEno: same logic used for FWI95d
# ============================================================

CorrEno <- function(exp, obs,
                    time_dim = "sdate",
                    member_dim = NULL,
                    method = "pearson",
                    alpha = 0.05,
                    test.type = "two-sided",
                    pval = TRUE,
                    handle.na = "return.na",
                    ncores = 1) {

  if (!is.null(member_dim)) {
    exp <- multiApply::Apply(
      data = exp,
      target_dims = member_dim,
      fun = mean,
      na.rm = FALSE,
      ncores = ncores
    )$output1
  }

  output <- multiApply::Apply(
    data = list(exp = exp, obs = obs),
    target_dims = time_dim,
    fun = .CorrEno,
    time_dim = time_dim,
    method = method,
    alpha = alpha,
    test.type = test.type,
    pval = pval,
    handle.na = handle.na,
    ncores = ncores
  )

  return(output)
}


.CorrEno <- function(exp, obs, time_dim, method,
                     alpha, test.type, pval, handle.na) {

  calc <- function(exp, obs) {

    r <- cor(exp, obs, method = method)

    xobs <- obs
    dim(xobs) <- length(xobs)
    names(dim(xobs)) <- time_dim

    n_eff <- s2dv::Eno(
      data = xobs,
      time_dim = time_dim,
      na.action = na.pass,
      ncores = 1
    )

    if (test.type == "one-sided") {

      tcrit <- qt(alpha, df = n_eff - 2, lower.tail = FALSE)
      tt <- r * sqrt(n_eff - 2) / sqrt(1 - r^2)

      sign <- is.finite(tt) &&
              is.finite(tcrit) &&
              tt >= tcrit &&
              r > 0

      pv <- pt(tt, df = n_eff - 2, lower.tail = FALSE)

    } else {

      tcrit <- qt(alpha / 2, df = n_eff - 2, lower.tail = FALSE)
      tt <- abs(r) * sqrt(n_eff - 2) / sqrt(1 - r^2)

      sign <- is.finite(tt) &&
              is.finite(tcrit) &&
              tt >= tcrit

      pv <- 2 * pt(tt, df = n_eff - 2, lower.tail = FALSE)
    }

    list(r = r, sign = sign, pval = pv)
  }


  if (anyNA(exp) || anyNA(obs)) {

    if (handle.na == "only.complete.pairs") {

      ok <- complete.cases(exp, obs)

      if (sum(ok) < 3)
        return(list(r = NA_real_,
                    sign = FALSE,
                    pval = NA_real_))

      return(calc(exp[ok], obs[ok]))

    } else {

      return(list(r = NA_real_,
                  sign = FALSE,
                  pval = NA_real_))
    }

  }

  calc(exp, obs)
}


# ============================================================
# Helpers
# ============================================================

read_nc <- function(file, var) {

  nc <- nc_open(file)

  if (!(var %in% names(nc$var))) {
    stop("Variable ", var, " not found in ", file,
         "\nVariables: ", paste(names(nc$var), collapse = ", "))
  }

  x <- ncvar_get(nc, var)

  lon <- nc$dim$lon$vals
  lat <- nc$dim$lat$vals

  nc_close(nc)

  list(data = x, lon = lon, lat = lat)
}


set_dims3 <- function(x, lon, lat) {

  if (!all(dim(x) == c(length(lon), length(lat), ntime))) {
    stop("Unexpected dimensions: ",
         paste(dim(x), collapse = " x "))
  }

  dimnames(x) <- list(
    lon   = as.character(lon),
    lat   = as.character(lat),
    sdate = as.character(years)
  )

  d <- dim(x)
  names(d) <- c("lon", "lat", "sdate")
  dim(x) <- d

  x
}


set_dims4 <- function(x, lon, lat, members) {

  dimnames(x) <- list(
    lon    = as.character(lon),
    lat    = as.character(lat),
    sdate  = as.character(years),
    member = members
  )

  d <- dim(x)
  names(d) <- c("lon", "lat", "sdate", "member")
  dim(x) <- d

  x
}


# ============================================================
# Land mask from FWI95d OBS
# Ensures spatial summaries/maps use the same FWI domain
# ============================================================

fwi_obs_file <- file.path(
  data_root,
  "verification_new",
  "FWI95d_OBS_runmean5_1963-2022_corrected.nc"
)

if (!file.exists(fwi_obs_file))
  stop("FWI95d OBS file not found: ", fwi_obs_file)

nc <- nc_open(fwi_obs_file)

fwi_var <- if ("FWI95d" %in% names(nc$var)) {
  "FWI95d"
} else if ("FWI" %in% names(nc$var)) {
  "FWI"
} else {
  names(nc$var)[1]
}

fwi_obs <- ncvar_get(nc, fwi_var)
nc_close(nc)

land_mask <- is.finite(fwi_obs[,,1])

cat("FWI land grid cells:", sum(land_mask), "\n\n")


# ============================================================
# Main calculation
# ============================================================

results <- list()
summary_table <- data.frame()


for (var in vars) {

  cat("\n====================================================\n")
  cat("VARIABLE:", var, "\n")
  cat("====================================================\n")


  # ----------------------------------------------------------
  # OBS
  # ----------------------------------------------------------

  obs_file <- file.path(
    driver_dir, "OBS", var,
    paste0(var, "_OBS_FY1-FY5_1963-2022.nc")
  )

  stopifnot(file.exists(obs_file))

  oo <- read_nc(obs_file, var)

  lon <- oo$lon
  lat <- oo$lat

  obs <- set_dims3(oo$data, lon, lat)

  nlon <- length(lon)
  nlat <- length(lat)


  # ----------------------------------------------------------
  # INIT
  # ----------------------------------------------------------

  init_files <- list.files(
    file.path(driver_dir, "INIT", var),
    pattern = paste0("^", var,
                     "_INIT_FY1-FY5_r[0-9]+i1p1f1_1963-2022\\.nc$"),
    full.names = TRUE
  )

  init_files <- init_files[
    order(
      as.numeric(
        sub(".*_r([0-9]+)i1p1f1_.*", "\\1",
            basename(init_files))
      )
    )
  ]

  if (length(init_files) != 24)
    stop(var, ": expected 24 INIT files, found ",
         length(init_files))


  init <- array(
    NA_real_,
    dim = c(nlon, nlat, ntime, length(init_files))
  )

  init_members <- character(length(init_files))

  for (i in seq_along(init_files)) {

    nc <- nc_open(init_files[i])
    xx <- ncvar_get(nc, var)
    nc_close(nc)

    if (!all(dim(xx) == c(nlon, nlat, ntime)))
      stop("INIT dimensions mismatch: ", init_files[i])

    init[,,,i] <- xx

    init_members[i] <- sub(
      paste0(var, "_INIT_FY1-FY5_(r[0-9]+i1p1f1)_.*"),
      "\\1",
      basename(init_files[i])
    )
  }

  init <- set_dims4(init, lon, lat, init_members)


  # ----------------------------------------------------------
  # HIST
  # ----------------------------------------------------------

  hist_files <- list.files(
    file.path(driver_dir, "HIST", var),
    pattern = paste0("^", var,
                     "_HIST_FY1-FY5_r[0-9]+i1p2f1_1963-2022\\.nc$"),
    full.names = TRUE
  )

  hist_files <- hist_files[
    order(
      as.numeric(
        sub(".*_r([0-9]+)i1p2f1_.*", "\\1",
            basename(hist_files))
      )
    )
  ]

  if (length(hist_files) != 10)
    stop(var, ": expected 10 HIST files, found ",
         length(hist_files))


  hist <- array(
    NA_real_,
    dim = c(nlon, nlat, ntime, length(hist_files))
  )

  hist_members <- character(length(hist_files))

  for (i in seq_along(hist_files)) {

    nc <- nc_open(hist_files[i])
    xx <- ncvar_get(nc, var)
    nc_close(nc)

    if (!all(dim(xx) == c(nlon, nlat, ntime)))
      stop("HIST dimensions mismatch: ", hist_files[i])

    hist[,,,i] <- xx

    hist_members[i] <- sub(
      paste0(var, "_HIST_FY1-FY5_(r[0-9]+i1p2f1)_.*"),
      "\\1",
      basename(hist_files[i])
    )
  }

  hist <- set_dims4(hist, lon, lat, hist_members)


  cat("OBS :", paste(dim(obs), collapse = " x "), "\n")
  cat("INIT:", paste(dim(init), collapse = " x "), "\n")
  cat("HIST:", paste(dim(hist), collapse = " x "), "\n")


  # ----------------------------------------------------------
  # Correlation skill
  # ----------------------------------------------------------

  corr_init <- CorrEno(
    exp = init,
    obs = obs,
    time_dim = "sdate",
    member_dim = "member",
    method = "pearson",
    alpha = 0.05,
    test.type = "two-sided",
    pval = TRUE,
    ncores = 8
  )


  corr_hist <- CorrEno(
    exp = hist,
    obs = obs,
    time_dim = "sdate",
    member_dim = "member",
    method = "pearson",
    alpha = 0.05,
    test.type = "two-sided",
    pval = TRUE,
    ncores = 8
  )


  # Same land domain as FWI95d
  r_init <- corr_init$r
  r_hist <- corr_hist$r

  r_init[!land_mask] <- NA
  r_hist[!land_mask] <- NA

  sig_init <- corr_init$sign
  sig_hist <- corr_hist$sign

  sig_init[!land_mask] <- FALSE
  sig_hist[!land_mask] <- FALSE


  # ----------------------------------------------------------
  # Summary
  # ----------------------------------------------------------

  valid <- land_mask & is.finite(r_init)

  ss <- data.frame(
    variable = var,

    INIT_mean_r =
      mean(r_init[valid], na.rm = TRUE),

    INIT_median_r =
      median(r_init[valid], na.rm = TRUE),

    INIT_fraction_positive =
      mean(r_init[valid] > 0, na.rm = TRUE),

    INIT_fraction_significant =
      mean(sig_init[valid], na.rm = TRUE),

    HIST_mean_r =
      mean(r_hist[valid], na.rm = TRUE),

    HIST_median_r =
      median(r_hist[valid], na.rm = TRUE),

    HIST_fraction_positive =
      mean(r_hist[valid] > 0, na.rm = TRUE),

    HIST_fraction_significant =
      mean(sig_hist[valid], na.rm = TRUE)
  )

  print(ss)

  summary_table <- rbind(summary_table, ss)


  results[[var]] <- list(
    lon = lon,
    lat = lat,
    obs = obs,
    corr_init = corr_init,
    corr_hist = corr_hist,
    r_init_land = r_init,
    r_hist_land = r_hist,
    sig_init_land = sig_init,
    sig_hist_land = sig_hist
  )
}


# ============================================================
# Save numerical results
# ============================================================

saveRDS(
  results,
  file.path(out_dir, "driver_skill_FY1-FY5_corrected.rds")
)

write.csv(
  summary_table,
  file.path(out_dir, "driver_skill_summary_corrected.csv"),
  row.names = FALSE
)


# ============================================================
# Four main maps: INIT vs OBS
# ============================================================

pdf(
  file.path(out_dir, "driver_skill_INIT_vs_OBS_maps_corrected.pdf"),
  width = 8,
  height = 6
)

for (var in vars) {

  rr <- results[[var]]

  r <- rr$r_init_land
  sig <- rr$sig_init_land

  image.plot(
    rr$lon,
    rr$lat,
    r,
    zlim = c(-1, 1),
    main = paste0(var, ": INIT vs ERA5, FY1-FY5"),
    xlab = "Longitude",
    ylab = "Latitude",
    legend.lab = "r"
  )

  map("world", add = TRUE, lwd = 0.8)

  idx <- which(sig & is.finite(r), arr.ind = TRUE)

  if (nrow(idx) > 0) {
    points(
      rr$lon[idx[,1]],
      rr$lat[idx[,2]],
      pch = 1,
      cex = 0.45,
      lwd = 0.7
    )
  }
}

dev.off()


# ============================================================
# Print final table
# ============================================================

cat("\n\n====================================================\n")
cat("FINAL SUMMARY — LAND GRID POINTS\n")
cat("====================================================\n\n")

print(summary_table, digits = 3)

cat("\nSaved:\n")
cat(file.path(out_dir, "driver_skill_FY1-FY5_corrected.rds"), "\n")
cat(file.path(out_dir, "driver_skill_summary_corrected.csv"), "\n")
cat(file.path(out_dir, "driver_skill_INIT_vs_OBS_maps_corrected.pdf"), "\n")
