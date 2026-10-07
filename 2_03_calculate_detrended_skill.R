rm(list = ls())
gc()

library(ncdf4)
library(s2dv)

source("../functions/CorrEno.R")

# ============================================================
# CONFIGURACIÓN
# ============================================================

VER <- "/Users/marco/Library/CloudStorage/Dropbox/estcena/scripts/ONFIRE/decadal/data/verification_new"
OUT <- file.path(VER, "results")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

years <- 1963:2022
tt <- seq_along(years)

# Baseline climatológica 1991-2020
clim_idx <- which(years >= 1991 & years <= 2020)

cat("Climatology:", years[min(clim_idx)], "-", years[max(clim_idx)], "\n")
cat("Indices:", min(clim_idx), "-", max(clim_idx), "\n")
cat("N:", length(clim_idx), "\n\n")


# ============================================================
# FUNCIONES DE DETRENDING
# ============================================================

# Array [lon, lat, time]
detrend3 <- function(x) {

  out <- x

  for(i in seq_len(dim(x)[1])) {
    for(j in seq_len(dim(x)[2])) {

      y <- x[i, j, ]
      ok <- is.finite(y)

      if(sum(ok) >= 3) {

        fit <- lm.fit(
          cbind(1, tt[ok]),
          y[ok]
        )

        slope <- fit$coefficients[2]
        t0 <- mean(tt[ok])

        # Remove only the centered linear trend,
        # preserving the temporal mean
        out[i, j, ok] <-
          y[ok] - slope * (tt[ok] - t0)
      }
    }
  }

  out
}


# Array [lon, lat, time, member]
detrend4 <- function(x) {

  out <- x

  # Ensemble mean [lon, lat, time]
  ensmean <- apply(
    x,
    c(1, 2, 3),
    mean,
    na.rm = TRUE
  )

  ensmean[!is.finite(ensmean)] <- NA_real_

  for(i in seq_len(dim(x)[1])) {
    for(j in seq_len(dim(x)[2])) {

      ymean <- ensmean[i, j, ]
      okmean <- is.finite(ymean)

      if(sum(okmean) >= 3) {

        # Common trend estimated from ensemble mean
        fit <- lm.fit(
          cbind(1, tt[okmean]),
          ymean[okmean]
        )

        slope <- fit$coefficients[2]
        t0 <- mean(tt[okmean])

        common_trend <- slope * (tt - t0)

        # Remove the same trend from every member
        for(m in seq_len(dim(x)[4])) {

          y <- x[i, j, , m]
          ok <- is.finite(y)

          out[i, j, ok, m] <-
            y[ok] - common_trend[ok]
        }
      }
    }
  }

  out
}


# ============================================================
# 1. OBS
# ============================================================

obs_file <- file.path(
  VER,
  "FWI95d_OBS_runmean5_1963-2022_corrected.nc"
)

nc <- nc_open(obs_file)

lon <- ncvar_get(nc, "lon")
lat <- ncvar_get(nc, "lat")
obs <- ncvar_get(nc, "FWI95d")

nc_close(nc)

nlon <- length(lon)
nlat <- length(lat)
ntime <- length(years)

stopifnot(all(dim(obs) == c(nlon, nlat, ntime)))

dimnames(obs) <- list(
  lon   = lon,
  lat   = lat,
  sdate = as.character(years)
)

d <- dim(obs)
names(d) <- c("lon", "lat", "sdate")
dim(obs) <- d


# ============================================================
# 2. INIT
# ============================================================

init_files <- list.files(
  file.path(VER, "INIT"),
  pattern = "FWI95d_INIT_FY1-FY5_r[0-9]+i1p1f1_1963-2022\\.nc$",
  full.names = TRUE
)

stopifnot(length(init_files) == 24)

dcpp <- array(
  NA_real_,
  dim = c(nlon, nlat, ntime, length(init_files))
)

init_members <- character(length(init_files))

for(m in seq_along(init_files)) {

  nc <- nc_open(init_files[m])
  x <- ncvar_get(nc, "FWI95d")
  nc_close(nc)

  stopifnot(all(dim(x) == c(nlon, nlat, ntime)))

  dcpp[, , , m] <- x

  init_members[m] <- sub(
    ".*_(r[0-9]+i1p1f1)_1963.*",
    "\\1",
    basename(init_files[m])
  )
}

dimnames(dcpp) <- list(
  lon    = lon,
  lat    = lat,
  sdate  = as.character(years),
  member = init_members
)

d <- dim(dcpp)
names(d) <- c("lon", "lat", "sdate", "member")
dim(dcpp) <- d


# ============================================================
# 3. HIST+245
# ============================================================

hist_files <- list.files(
  file.path(VER, "HIST"),
  pattern = "FWI95d_HIST_runmean5_r[0-9]+i1p2f1_1963-2022\\.nc$",
  full.names = TRUE
)

stopifnot(length(hist_files) == 10)

hist <- array(
  NA_real_,
  dim = c(nlon, nlat, ntime, length(hist_files))
)

hist_members <- character(length(hist_files))

for(m in seq_along(hist_files)) {

  nc <- nc_open(hist_files[m])
  x <- ncvar_get(nc, "FWI95d")
  nc_close(nc)

  stopifnot(all(dim(x) == c(nlon, nlat, ntime)))

  hist[, , , m] <- x

  hist_members[m] <- sub(
    ".*_(r[0-9]+i1p2f1)_1963.*",
    "\\1",
    basename(hist_files[m])
  )
}

dimnames(hist) <- list(
  lon    = lon,
  lat    = lat,
  sdate  = as.character(years),
  member = hist_members
)

d <- dim(hist)
names(d) <- c("lon", "lat", "sdate", "member")
dim(hist) <- d


cat("OBS :", paste(dim(obs), collapse = " x "), "\n")
cat("INIT:", paste(dim(dcpp), collapse = " x "), "\n")
cat("HIST:", paste(dim(hist), collapse = " x "), "\n\n")


# ============================================================
# 4. DETREND
# ============================================================

cat("Detrending OBS...\n")
obs_dt <- detrend3(obs)

cat("\nDetrending INIT...\n")
dcpp_dt <- detrend4(dcpp)

cat("\nDetrending HIST...\n")
hist_dt <- detrend4(hist)


# Restaurar nombres de dimensiones

dimnames(obs_dt) <- dimnames(obs)
d <- dim(obs_dt)
names(d) <- c("lon", "lat", "sdate")
dim(obs_dt) <- d

dimnames(dcpp_dt) <- dimnames(dcpp)
d <- dim(dcpp_dt)
names(d) <- c("lon", "lat", "sdate", "member")
dim(dcpp_dt) <- d

dimnames(hist_dt) <- dimnames(hist)
d <- dim(hist_dt)
names(d) <- c("lon", "lat", "sdate", "member")
dim(hist_dt) <- d


# ============================================================
# 5. DETRENDED INIT vs OBS CORRELATION
# ============================================================

corr <- CorrEno(
  exp        = dcpp_dt,
  obs        = obs_dt,
  time_dim   = "sdate",
  member_dim = 4,
  method     = "pearson",
  alpha      = 0.05,
  test.type = "two-sided",
  ncores     = 8,
  pval       = TRUE
)

cat("\n========================================\n")
cat("DETRENDED INIT-OBS correlation\n")
cat("========================================\n")

cat("mean r =", mean(corr$r, na.rm = TRUE), "\n")
cat("median r =", median(corr$r, na.rm = TRUE), "\n")
cat("fraction r>0 =", mean(corr$r > 0, na.rm = TRUE), "\n")
cat("fraction significant =", mean(corr$sign, na.rm = TRUE), "\n")
cat("significant AND positive =",
    mean(corr$sign & corr$r > 0, na.rm = TRUE), "\n")
cat("significant AND negative =",
    mean(corr$sign & corr$r < 0, na.rm = TRUE), "\n")


# ============================================================
# 6. DETRENDED RESIDUAL CORRELATION
# ============================================================

resid_corr <- s2dv::ResidualCorr(
  exp       = dcpp_dt,
  obs       = obs_dt,
  ref       = hist_dt,
  time_dim  = "sdate",
  memb_dim  = "member",
  method    = "pearson",
  handle.na = "return.na",
  pval      = FALSE,
  sign      = TRUE,
  alpha     = 0.05,
  ncores    = 8
)

cat("\n========================================\n")
cat("DETRENDED residual correlation\n")
cat("========================================\n")

cat("mean r =", mean(resid_corr$res.corr, na.rm = TRUE), "\n")
cat("median r =", median(resid_corr$res.corr, na.rm = TRUE), "\n")
cat("fraction >0 =", mean(resid_corr$res.corr > 0, na.rm = TRUE), "\n")
cat("fraction significant =", mean(resid_corr$sign, na.rm = TRUE), "\n")
cat("significant AND positive =",
    mean(resid_corr$sign & resid_corr$res.corr > 0, na.rm = TRUE), "\n")
cat("significant AND negative =",
    mean(resid_corr$sign & resid_corr$res.corr < 0, na.rm = TRUE), "\n")


# ============================================================
# 7. DETRENDED DIFFERENCE IN CORRELATION
#    INIT - HIST
# ============================================================

diff_corr <- s2dv::DiffCorr(
  exp       = dcpp_dt,
  ref       = hist_dt,
  obs       = obs_dt,
  time_dim  = "sdate",
  memb_dim  = "member",
  method    = "pearson",
  handle.na = "return.na",
  pval      = FALSE,
  sign      = TRUE,
  alpha     = 0.05,
  test.type = "two-sided",
  ncores    = 8
)

cat("\n========================================\n")
cat("DETRENDED difference correlation INIT-HIST\n")
cat("========================================\n")

cat("mean difference =", mean(diff_corr$diff.corr, na.rm = TRUE), "\n")
cat("median difference =", median(diff_corr$diff.corr, na.rm = TRUE), "\n")
cat("fraction >0 =", mean(diff_corr$diff.corr > 0, na.rm = TRUE), "\n")
cat("fraction significant =", mean(diff_corr$sign, na.rm = TRUE), "\n")
cat("significant AND INIT better =",
    mean(diff_corr$sign & diff_corr$diff.corr > 0, na.rm = TRUE), "\n")
cat("significant AND HIST better =",
    mean(diff_corr$sign & diff_corr$diff.corr < 0, na.rm = TRUE), "\n")


# ============================================================
# 8. DETRENDED RPSS: INIT vs CLIMATOLOGY
# ============================================================

rpss_clim_dt <- s2dv::RPSS(
  exp              = dcpp_dt,
  obs              = obs_dt,
  ref              = NULL,
  prob_thresholds  = c(1/3, 2/3),
  indices_for_clim = clim_idx,
  Fair             = FALSE,
  sig_method.type  = "two.sided.approx",
  alpha            = 0.05,
  time_dim         = "sdate",
  memb_dim         = "member",
  dat_dim          = NULL,
  ncores           = 8
)

cat("\n========================================\n")
cat("DETRENDED RPSS INIT vs climatology\n")
cat("========================================\n")

cat("mean RPSS =", mean(rpss_clim_dt$rpss, na.rm = TRUE), "\n")
cat("median RPSS =", median(rpss_clim_dt$rpss, na.rm = TRUE), "\n")
cat("fraction >0 =", mean(rpss_clim_dt$rpss > 0, na.rm = TRUE), "\n")
cat("fraction significant =", mean(rpss_clim_dt$sign, na.rm = TRUE), "\n")
cat("significant AND positive =",
    mean(rpss_clim_dt$sign & rpss_clim_dt$rpss > 0, na.rm = TRUE), "\n")
cat("significant AND negative =",
    mean(rpss_clim_dt$sign & rpss_clim_dt$rpss < 0, na.rm = TRUE), "\n")


# ============================================================
# 9. DETRENDED RPSS: INIT vs HIST+245
# ============================================================

rpss_hist_dt <- s2dv::RPSS(
  exp              = dcpp_dt,
  obs              = obs_dt,
  ref              = hist_dt,
  prob_thresholds  = c(1/3, 2/3),
  indices_for_clim = clim_idx,
  Fair             = FALSE,
  sig_method.type  = "two.sided.approx",
  alpha            = 0.05,
  time_dim         = "sdate",
  memb_dim         = "member",
  dat_dim          = NULL,
  ncores           = 8
)

cat("\n========================================\n")
cat("DETRENDED RPSS INIT vs HIST+245\n")
cat("========================================\n")

cat("mean RPSS =", mean(rpss_hist_dt$rpss, na.rm = TRUE), "\n")
cat("median RPSS =", median(rpss_hist_dt$rpss, na.rm = TRUE), "\n")
cat("fraction >0 =", mean(rpss_hist_dt$rpss > 0, na.rm = TRUE), "\n")
cat("fraction significant =", mean(rpss_hist_dt$sign, na.rm = TRUE), "\n")
cat("significant AND positive =",
    mean(rpss_hist_dt$sign & rpss_hist_dt$rpss > 0, na.rm = TRUE), "\n")
cat("significant AND negative =",
    mean(rpss_hist_dt$sign & rpss_hist_dt$rpss < 0, na.rm = TRUE), "\n")


# ============================================================
# 10. GUARDAR RESULTADOS
# ============================================================

saveRDS(
  list(
    lon         = lon,
    lat         = lat,
    years       = years,
    clim_years  = years[clim_idx],
    corr        = corr,
    residual    = resid_corr,
    diffcorr    = diff_corr,
    rpss_clim   = rpss_clim_dt,
    rpss_hist   = rpss_hist_dt
  ),
  file = file.path(
    OUT,
    "skill_detrended_commontrend_corrected.rds"
  )
)

cat("\n========================================\n")
cat("DONE\n")
cat("Saved:", file.path(OUT, "skill_detrended_commontrend_corrected.rds"), "\n")
cat("========================================\n")
