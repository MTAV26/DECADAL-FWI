
rm(list = ls())
graphics.off()

suppressPackageStartupMessages({
  library(ncdf4)
})

# ============================================================
# FWI95d linear trends
#
# OBS / INIT / HIST+245
# Verification period: 1963–2022
# Units: FWI95d days per decade
#
# IMPORTANT:
# - INIT trend is calculated from the ensemble-mean FY1–FY5 series
# - HIST trend is calculated from the ensemble-mean HIST+245 series
# - Member-specific trends are also saved for diagnostics
# ============================================================

PROJECT_ROOT <-
  "/Users/marco/Library/CloudStorage/Dropbox/estcena/scripts/ONFIRE/decadal"

DATA_ROOT <- file.path(
  PROJECT_ROOT,
  "data",
  "verification_new"
)

OUT <- file.path(
  DATA_ROOT,
  "results"
)

dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

YEARS <- 1963:2022
NTIME <- length(YEARS)

# ============================================================
# Helper: read FWI95d NetCDF
# ============================================================

read_fwi <- function(file) {

  nc <- nc_open(file)

  on.exit(nc_close(nc))

  lon <- ncvar_get(nc, "lon")
  lat <- ncvar_get(nc, "lat")

  vars <- names(nc$var)

  if ("FWI" %in% vars) {
    vname <- "FWI"
  } else if ("FWI95d" %in% vars) {
    vname <- "FWI95d"
  } else {
    stop(
      "Cannot identify FWI variable in: ",
      basename(file),
      "\nVariables: ",
      paste(vars, collapse = ", ")
    )
  }

  x <- ncvar_get(nc, vname)

  list(
    lon = lon,
    lat = lat,
    data = x
  )
}

# ============================================================
# Helper: linear trend
# slope returned in DAYS PER DECADE
# ============================================================

trend_one <- function(y, years = YEARS) {

  ok <- is.finite(y)

  if (sum(ok) < 10)
    return(NA_real_)

  fit <- lm(y[ok] ~ years[ok])

  as.numeric(coef(fit)[2]) * 10
}

trend_field <- function(x, years = YEARS) {

  # x = lon x lat x time

  stopifnot(
    length(dim(x)) == 3,
    dim(x)[3] == length(years)
  )

  apply(
    x,
    c(1, 2),
    trend_one,
    years = years
  )
}

# ============================================================
# 1. OBS
# ============================================================

obs_file <- file.path(
  DATA_ROOT,
  "FWI95d_OBS_runmean5_1963-2022_corrected.nc"
)

stopifnot(file.exists(obs_file))

cat("\nReading OBS...\n")

oo <- read_fwi(obs_file)

lon <- as.numeric(oo$lon)
lat <- as.numeric(oo$lat)
obs <- oo$data

stopifnot(
  all(dim(obs) == c(length(lon), length(lat), NTIME))
)

cat(
  "OBS dimensions:",
  paste(dim(obs), collapse = " x "),
  "\n"
)

# ============================================================
# 2. INIT
# ============================================================

init_files <- list.files(
  file.path(DATA_ROOT, "INIT"),
  pattern = "^FWI95d_INIT_FY1-FY5_.*_1963-2022\\.nc$",
  full.names = TRUE
)

init_files <- sort(init_files)

cat("\nINIT files:", length(init_files), "\n")

stopifnot(length(init_files) == 24)

init <- array(
  NA_real_,
  dim = c(
    length(lon),
    length(lat),
    NTIME,
    length(init_files)
  )
)

init_members <- character(length(init_files))

for (i in seq_along(init_files)) {

  cat(
    sprintf(
      "Reading INIT %02d/%02d: %s\n",
      i,
      length(init_files),
      basename(init_files[i])
    )
  )

  xx <- read_fwi(init_files[i])

  stopifnot(
    identical(as.numeric(xx$lon), lon),
    identical(as.numeric(xx$lat), lat),
    all(dim(xx$data) ==
          c(length(lon), length(lat), NTIME))
  )

  init[, , , i] <- xx$data

  init_members[i] <- sub(
    "^FWI95d_INIT_FY1-FY5_(.*)_1963-2022\\.nc$",
    "\\1",
    basename(init_files[i])
  )
}

# Ensemble mean at each time step
init_mean <- apply(
  init,
  c(1, 2, 3),
  mean,
  na.rm = TRUE
)

init_mean[!is.finite(init_mean)] <- NA_real_

# ============================================================
# 3. HIST+245
# ============================================================

hist_files <- list.files(
  file.path(DATA_ROOT, "HIST"),
  pattern = "^FWI95d_HIST_runmean5_.*_1963-2022\\.nc$",
  full.names = TRUE
)

hist_files <- sort(hist_files)

cat("\nHIST files:", length(hist_files), "\n")

stopifnot(length(hist_files) == 10)

hist <- array(
  NA_real_,
  dim = c(
    length(lon),
    length(lat),
    NTIME,
    length(hist_files)
  )
)

hist_members <- character(length(hist_files))

for (i in seq_along(hist_files)) {

  cat(
    sprintf(
      "Reading HIST %02d/%02d: %s\n",
      i,
      length(hist_files),
      basename(hist_files[i])
    )
  )

  xx <- read_fwi(hist_files[i])

  stopifnot(
    identical(as.numeric(xx$lon), lon),
    identical(as.numeric(xx$lat), lat),
    all(dim(xx$data) ==
          c(length(lon), length(lat), NTIME))
  )

  hist[, , , i] <- xx$data

  hist_members[i] <- sub(
    "^FWI95d_HIST_runmean5_(.*)_1963-2022\\.nc$",
    "\\1",
    basename(hist_files[i])
  )
}

hist_mean <- apply(
  hist,
  c(1, 2, 3),
  mean,
  na.rm = TRUE
)

hist_mean[!is.finite(hist_mean)] <- NA_real_

# ============================================================
# 4. Ensemble-mean trends
# ============================================================

cat("\nCalculating OBS trend...\n")
trend_obs <- trend_field(obs)

cat("Calculating INIT ensemble-mean trend...\n")
trend_init <- trend_field(init_mean)

cat("Calculating HIST ensemble-mean trend...\n")
trend_hist <- trend_field(hist_mean)

# ============================================================
# 5. Member-specific trends
# ============================================================

cat("\nCalculating INIT member trends...\n")

trend_init_members <- array(
  NA_real_,
  dim = c(
    length(lon),
    length(lat),
    length(init_files)
  )
)

for (m in seq_along(init_files)) {

  trend_init_members[, , m] <-
    trend_field(init[, , , m])
}

cat("Calculating HIST member trends...\n")

trend_hist_members <- array(
  NA_real_,
  dim = c(
    length(lon),
    length(lat),
    length(hist_files)
  )
)

for (m in seq_along(hist_files)) {

  trend_hist_members[, , m] <-
    trend_field(hist[, , , m])
}

# ============================================================
# 6. Ensemble spread of member trends
# ============================================================

qfun <- function(x, p) {

  if (all(!is.finite(x)))
    return(NA_real_)

  quantile(
    x,
    probs = p,
    na.rm = TRUE,
    names = FALSE
  )
}

init_trend_p05 <- apply(
  trend_init_members,
  c(1, 2),
  qfun,
  p = 0.05
)

init_trend_p95 <- apply(
  trend_init_members,
  c(1, 2),
  qfun,
  p = 0.95
)

hist_trend_p05 <- apply(
  trend_hist_members,
  c(1, 2),
  qfun,
  p = 0.05
)

hist_trend_p95 <- apply(
  trend_hist_members,
  c(1, 2),
  qfun,
  p = 0.95
)

# ============================================================
# 7. Differences in trends
# ============================================================

trend_diff_init_obs <- trend_init - trend_obs
trend_diff_hist_obs <- trend_hist - trend_obs
trend_diff_init_hist <- trend_init - trend_hist

# ============================================================
# 8. Summary
# ============================================================

summ <- function(x) {

  x <- x[is.finite(x)]

  c(
    mean = mean(x),
    median = median(x),
    pct_positive = mean(x > 0) * 100,
    pct_negative = mean(x < 0) * 100
  )
}

summary_table <- rbind(
  OBS = summ(trend_obs),
  INIT = summ(trend_init),
  HIST245 = summ(trend_hist),
  INIT_minus_OBS = summ(trend_diff_init_obs),
  HIST245_minus_OBS = summ(trend_diff_hist_obs),
  INIT_minus_HIST245 = summ(trend_diff_init_hist)
)

cat("\n====================================================\n")
cat("FWI95d LINEAR TRENDS — days per decade\n")
cat("====================================================\n")

print(summary_table)

# ============================================================
# 9. Save
# ============================================================

out_rds <- file.path(
  OUT,
  "FWI95d_trends_1963-2022_corrected.rds"
)

saveRDS(
  list(
    lon = lon,
    lat = lat,
    years = YEARS,

    trend_obs = trend_obs,
    trend_init = trend_init,
    trend_hist = trend_hist,

    trend_diff_init_obs = trend_diff_init_obs,
    trend_diff_hist_obs = trend_diff_hist_obs,
    trend_diff_init_hist = trend_diff_init_hist,

    trend_init_members = trend_init_members,
    trend_hist_members = trend_hist_members,

    init_members = init_members,
    hist_members = hist_members,

    init_trend_p05 = init_trend_p05,
    init_trend_p95 = init_trend_p95,
    hist_trend_p05 = hist_trend_p05,
    hist_trend_p95 = hist_trend_p95,

    summary = summary_table,

    units = "FWI95d days per decade"
  ),
  out_rds
)

write.csv(
  summary_table,
  file.path(
    OUT,
    "FWI95d_trends_summary_1963-2022_corrected.csv"
  )
)

cat("\nSaved:\n")
cat(out_rds, "\n")
cat(
  file.path(
    OUT,
    "FWI95d_trends_summary_1963-2022_corrected.csv"
  ),
  "\n"
)

cat("\nDONE\n")

