rm(list=ls())
gc()

library(ncdf4)
library(s2dv)

VER <- "/Users/marco/Library/CloudStorage/Dropbox/estcena/scripts/ONFIRE/decadal/data/verification_new"
OUT <- file.path(VER, "results")
dir.create(OUT, showWarnings=FALSE, recursive=TRUE)

years <- 1963:2022
clim_idx <- which(years >= 1991 & years <= 2020)

cat("Climatology years:", years[min(clim_idx)], "-", years[max(clim_idx)], "\n")
cat("Indices:", min(clim_idx), "-", max(clim_idx), "\n")
cat("N =", length(clim_idx), "\n\n")

# OBS
obs_file <- file.path(VER, "FWI95d_OBS_runmean5_1963-2022_corrected.nc")

nc <- nc_open(obs_file)
lon <- ncvar_get(nc, "lon")
lat <- ncvar_get(nc, "lat")
obs <- ncvar_get(nc, "FWI95d")
nc_close(nc)

nlon <- length(lon)
nlat <- length(lat)
ntime <- length(years)

dimnames(obs) <- list(
  lon=lon, lat=lat, sdate=as.character(years)
)
d <- dim(obs)
names(d) <- c("lon","lat","sdate")
dim(obs) <- d


# INIT
init_files <- list.files(
  file.path(VER,"INIT"),
  pattern="FWI95d_INIT_FY1-FY5_r[0-9]+i1p1f1_1963-2022\\.nc$",
  full.names=TRUE
)

stopifnot(length(init_files)==24)

dcpp <- array(NA_real_,
              dim=c(nlon,nlat,ntime,length(init_files)))

members <- character(length(init_files))

for(i in seq_along(init_files)) {
  nc <- nc_open(init_files[i])
  dcpp[,,,i] <- ncvar_get(nc,"FWI95d")
  nc_close(nc)

  members[i] <- sub(
    ".*_(r[0-9]+i1p1f1)_1963.*",
    "\\1",
    basename(init_files[i])
  )
}

dimnames(dcpp) <- list(
  lon=lon, lat=lat, sdate=as.character(years), member=members
)
d <- dim(dcpp)
names(d) <- c("lon","lat","sdate","member")
dim(dcpp) <- d


# HIST+245
hist_files <- list.files(
  file.path(VER,"HIST"),
  pattern="FWI95d_HIST_runmean5_r[0-9]+i1p2f1_1963-2022\\.nc$",
  full.names=TRUE
)

stopifnot(length(hist_files)==10)

hist <- array(NA_real_,
              dim=c(nlon,nlat,ntime,length(hist_files)))

members_hist <- character(length(hist_files))

for(i in seq_along(hist_files)) {
  nc <- nc_open(hist_files[i])
  hist[,,,i] <- ncvar_get(nc,"FWI95d")
  nc_close(nc)

  members_hist[i] <- sub(
    ".*_(r[0-9]+i1p2f1)_1963.*",
    "\\1",
    basename(hist_files[i])
  )
}

dimnames(hist) <- list(
  lon=lon, lat=lat, sdate=as.character(years), member=members_hist
)
d <- dim(hist)
names(d) <- c("lon","lat","sdate","member")
dim(hist) <- d


cat("OBS :", paste(dim(obs), collapse=" x "), "\n")
cat("INIT:", paste(dim(dcpp), collapse=" x "), "\n")
cat("HIST:", paste(dim(hist), collapse=" x "), "\n\n")


# ============================================================
# RPSS INIT vs CLIMATOLOGY
# ============================================================

rpss_clim <- s2dv::RPSS(
  exp               = dcpp,
  obs               = obs,
  ref               = NULL,
  prob_thresholds   = c(1/3, 2/3),
  indices_for_clim  = which(years %in% 1991:2020),
  Fair              = FALSE,
  time_dim          = "sdate",
  memb_dim          = "member",
  dat_dim           = NULL,
  sig_method.type   = "two.sided.approx",
  alpha             = 0.05,
  ncores            = 8
)

cat("RPSS INIT vs climatology\n")
cat("mean RPSS =", mean(rpss_clim$rpss, na.rm=TRUE), "\n")
cat("median RPSS =", median(rpss_clim$rpss, na.rm=TRUE), "\n")
cat("fraction >0 =", mean(rpss_clim$rpss > 0, na.rm=TRUE), "\n")
cat("fraction significant =", mean(rpss_clim$sign, na.rm=TRUE), "\n\n")


# ============================================================
# RPSS INIT vs HIST+245
# ============================================================

# RPSS INIT vs HIST+245
rpss_hist <- s2dv::RPSS(
  exp               = dcpp,
  obs               = obs,
  ref               = hist,
  prob_thresholds   = c(1/3, 2/3),
  indices_for_clim  = which(years %in% 1991:2020),
  Fair              = FALSE,
  time_dim          = "sdate",
  memb_dim          = "member",
  dat_dim           = NULL,
  sig_method.type   = "two.sided.approx",
  alpha             = 0.05,
  ncores            = 8
)

cat("RPSS INIT vs HIST+245\n")
cat("mean RPSS =", mean(rpss_hist$rpss, na.rm=TRUE), "\n")
cat("median RPSS =", median(rpss_hist$rpss, na.rm=TRUE), "\n")
cat("fraction >0 =", mean(rpss_hist$rpss > 0, na.rm=TRUE), "\n")
cat("fraction significant =", mean(rpss_hist$sign, na.rm=TRUE), "\n\n")


saveRDS(
  list(
    lon=lon,
    lat=lat,
    years=years,
    clim_years=years[clim_idx],
    rpss_clim=rpss_clim,
    rpss_hist=rpss_hist
  ),
  file=file.path(OUT,"RPSS_corrected.rds")
)

cat("DONE\n")
