rm(list=ls())
gc()

library(ncdf4)
library(s2dv)

source("../functions/CorrEno.R")

VER <- "/Users/marco/Library/CloudStorage/Dropbox/estcena/scripts/ONFIRE/decadal/data/verification_new"
OUT <- file.path(VER, "results")
dir.create(OUT, showWarnings=FALSE, recursive=TRUE)

# --------------------------------------------------
# OBS
# --------------------------------------------------
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
ntime <- dim(obs)[3]

stopifnot(ntime == 60)

dimnames(obs) <- list(
  lon=lon,
  lat=lat,
  sdate=as.character(1963:2022)
)
d <- dim(obs)
names(d) <- c("lon","lat","sdate")
dim(obs) <- d


# --------------------------------------------------
# INIT
# --------------------------------------------------
init_files <- list.files(
  file.path(VER,"INIT"),
  pattern="FWI95d_INIT_FY1-FY5_r[0-9]+i1p1f1_1963-2022\\.nc$",
  full.names=TRUE
)

stopifnot(length(init_files) == 24)

dcpp <- array(
  NA_real_,
  dim=c(nlon,nlat,ntime,length(init_files))
)

init_members <- character(length(init_files))

for(i in seq_along(init_files)) {

  nc <- nc_open(init_files[i])
  x  <- ncvar_get(nc,"FWI95d")
  nc_close(nc)

  stopifnot(all(dim(x) == c(nlon,nlat,ntime)))

  dcpp[,,,i] <- x

  init_members[i] <- sub(
    ".*_(r[0-9]+i1p1f1)_1963.*",
    "\\1",
    basename(init_files[i])
  )
}

dimnames(dcpp) <- list(
  lon=lon,
  lat=lat,
  sdate=as.character(1963:2022),
  member=init_members
)

d <- dim(dcpp)
names(d) <- c("lon","lat","sdate","member")
dim(dcpp) <- d


# --------------------------------------------------
# HIST+245
# --------------------------------------------------
hist_files <- list.files(
  file.path(VER,"HIST"),
  pattern="FWI95d_HIST_runmean5_r[0-9]+i1p2f1_1963-2022\\.nc$",
  full.names=TRUE
)

stopifnot(length(hist_files) == 10)

hist <- array(
  NA_real_,
  dim=c(nlon,nlat,ntime,length(hist_files))
)

hist_members <- character(length(hist_files))

for(i in seq_along(hist_files)) {

  nc <- nc_open(hist_files[i])
  x  <- ncvar_get(nc,"FWI95d")
  nc_close(nc)

  stopifnot(all(dim(x) == c(nlon,nlat,ntime)))

  hist[,,,i] <- x

  hist_members[i] <- sub(
    ".*_(r[0-9]+i1p2f1)_1963.*",
    "\\1",
    basename(hist_files[i])
  )
}

dimnames(hist) <- list(
  lon=lon,
  lat=lat,
  sdate=as.character(1963:2022),
  member=hist_members
)

d <- dim(hist)
names(d) <- c("lon","lat","sdate","member")
dim(hist) <- d


cat("OBS :", paste(dim(obs), collapse=" x "), "\n")
cat("INIT:", paste(dim(dcpp), collapse=" x "), "\n")
cat("HIST:", paste(dim(hist), collapse=" x "), "\n")


# --------------------------------------------------
# 1. INIT vs OBS correlation
# --------------------------------------------------
corr <- CorrEno(
  exp=dcpp,
  obs=obs,
  time_dim="sdate",
  member_dim=4,
  method="pearson",
  alpha=0.05,
  test.type = "two-sided",
  ncores=8,
  pval=TRUE
)

cat("\nINIT-OBS correlation\n")
cat("mean r =", mean(corr$r, na.rm=TRUE), "\n")
cat("median r =", median(corr$r, na.rm=TRUE), "\n")
cat("fraction r>0 =", mean(corr$r > 0, na.rm=TRUE), "\n")
cat("fraction significant =", mean(corr$sign, na.rm=TRUE), "\n")
cat("significant AND positive =",
    mean(corr$sign & corr$r > 0, na.rm=TRUE), "\n")
cat("significant AND negative =",
    mean(corr$sign & corr$r < 0, na.rm=TRUE), "\n")


# --------------------------------------------------
# 2. HIST+245 vs OBS correlation
# --------------------------------------------------

corr_hist <- CorrEno(
  exp=hist,
  obs=obs,
  time_dim="sdate",
  member_dim=4,
  method="pearson",
  alpha=0.05,
  test.type = "two-sided",
  ncores=8,
  pval=TRUE
)

cat("\nHIST+245-OBS correlation\n")
cat("mean r =", mean(corr_hist$r, na.rm=TRUE), "\n")
cat("median r =", median(corr_hist$r, na.rm=TRUE), "\n")
cat("fraction r>0 =", mean(corr_hist$r > 0, na.rm=TRUE), "\n")
cat("fraction significant =", mean(corr_hist$sign, na.rm=TRUE), "\n")
cat("significant AND positive =",
    mean(corr_hist$sign & corr_hist$r > 0, na.rm=TRUE), "\n")
cat("significant AND negative =",
    mean(corr_hist$sign & corr_hist$r < 0, na.rm=TRUE), "\n")


# --------------------------------------------------
# 3. Residual correlation
# --------------------------------------------------
resid_corr <- s2dv::ResidualCorr(
  exp=dcpp,
  obs=obs,
  ref=hist,
  time_dim="sdate",
  memb_dim="member",
  method="pearson",
  handle.na="return.na",
  pval=FALSE,
  sign=TRUE,
  alpha=0.05,
  ncores=8
)

cat("\nResidual correlation\n")
cat("mean r =", mean(resid_corr$res.corr, na.rm=TRUE), "\n")
cat("fraction >0 =", mean(resid_corr$res.corr > 0, na.rm=TRUE), "\n")
cat("fraction significant =", mean(resid_corr$sign, na.rm=TRUE), "\n")
cat("significant AND positive =",
    mean(resid_corr$sign & resid_corr$res.corr > 0, na.rm=TRUE), "\n")
cat("significant AND negative =",
    mean(resid_corr$sign & resid_corr$res.corr < 0, na.rm=TRUE), "\n")


# --------------------------------------------------
# 4. Difference in correlation INIT - HIST
# --------------------------------------------------
diff_corr <- s2dv::DiffCorr(
  exp=dcpp,
  ref=hist,
  obs=obs,
  time_dim="sdate",
  memb_dim="member",
  method="pearson",
  handle.na="return.na",
  pval=FALSE,
  sign=TRUE,
  alpha=0.05,
  test.type = "two-sided",
  ncores=8
)

cat("\nDifference correlation INIT-HIST\n")
cat("mean difference =", mean(diff_corr$diff.corr, na.rm=TRUE), "\n")
cat("fraction >0 =", mean(diff_corr$diff.corr > 0, na.rm=TRUE), "\n")
cat("fraction significant =", mean(diff_corr$sign, na.rm=TRUE), "\n")
cat("significant AND INIT better =",
    mean(diff_corr$sign & diff_corr$diff.corr > 0, na.rm=TRUE), "\n")
cat("significant AND HIST better =",
    mean(diff_corr$sign & diff_corr$diff.corr < 0, na.rm=TRUE), "\n")


saveRDS(
  list(
    lon=lon,
    lat=lat,
    years=1963:2022,
    corr=corr,
    corr_hist=corr_hist,
    residual=resid_corr,
    diffcorr=diff_corr
  ),
  file=file.path(OUT,"deterministic_skill_corrected.rds")
)

cat("\nDONE\n")
