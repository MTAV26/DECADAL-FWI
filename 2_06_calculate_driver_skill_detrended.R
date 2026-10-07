
rm(list = ls())
graphics.off()

suppressPackageStartupMessages({
  library(ncdf4)
  library(s2dv)
  library(multiApply)
})

# ============================================================
# Detrended prediction skill of FWI meteorological drivers
#
# Variables:
#   tas, hurs, pr, sfcWind
#
# Verification period: 1963-2022
#
# For each dataset, its own linear trend is removed before
# calculating Pearson correlation.
#
# INIT and HIST use the trend of their ensemble mean, which is removed from all members before
# calculating the ensemble-mean correlation.
# ============================================================

PROJECT_ROOT <-
  "/Users/marco/Library/CloudStorage/Dropbox/estcena/scripts/ONFIRE/decadal"

driver_dir <- file.path(
  PROJECT_ROOT,
  "data",
  "verification_drivers_new"
)

verification_dir <- file.path(
  PROJECT_ROOT,
  "data",
  "verification_new"
)

out_dir <- file.path(
  driver_dir,
  "results"
)

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

source(
  file.path(
    PROJECT_ROOT,
    "scripts",
    "functions",
    "CorrEno.R"
  )
)

years <- 1963:2022
ntime <- length(years)

variables <- c(
  "tas",
  "hurs",
  "pr",
  "sfcWind"
)

# ============================================================
# Helpers
# ============================================================

read_nc_field <- function(file, preferred_var = NULL) {

  nc <- nc_open(file)
  on.exit(nc_close(nc))

  lon <- ncvar_get(nc, "lon")
  lat <- ncvar_get(nc, "lat")

  vars <- names(nc$var)

  if (!is.null(preferred_var) &&
      preferred_var %in% vars) {

    vname <- preferred_var

  } else {

    candidates <- setdiff(
      vars,
      c("lon", "lat", "time")
    )

    if (length(candidates) < 1)
      stop("No data variable found in ", file)

    vname <- candidates[1]
  }

  x <- ncvar_get(nc, vname)

  list(
    lon = as.numeric(lon),
    lat = as.numeric(lat),
    data = x,
    variable = vname
  )
}


detrend_vector <- function(y) {

  if (all(!is.finite(y)))
    return(rep(NA_real_, length(y)))

  ok <- is.finite(y)

  if (sum(ok) < 10)
    return(rep(NA_real_, length(y)))

  out <- rep(NA_real_, length(y))

  fit <- lm(
    y[ok] ~ years[ok]
  )

  # Keep original mean after removing linear trend
  out[ok] <-
    residuals(fit) +
    mean(y[ok], na.rm = TRUE)

  out
}


detrend_3d <- function(x) {

  # x = lon x lat x time

  stopifnot(
    length(dim(x)) == 3,
    dim(x)[3] == ntime
  )

  out <- apply(
    x,
    c(1, 2),
    detrend_vector
  )

  # apply() returns time x lon x lat
  aperm(
    out,
    c(2, 3, 1)
  )
}


detrend_4d <- function(x) {

  # x = lon x lat x time x member
  stopifnot(
    length(dim(x)) == 4,
    dim(x)[3] == ntime
  )

  out <- x

  # Ensemble mean at each grid point and time
  ensmean <- apply(
    x,
    c(1, 2, 3),
    mean,
    na.rm = TRUE
  )

  ensmean[!is.finite(ensmean)] <- NA_real_

  for (i in seq_len(dim(x)[1])) {
    for (j in seq_len(dim(x)[2])) {

      ymean <- ensmean[i, j, ]
      okmean <- is.finite(ymean)

      if (sum(okmean) >= 10) {

        # Estimate one common linear trend from the ensemble mean
        fit <- lm(
          ymean[okmean] ~ years[okmean]
        )

        slope <- coef(fit)[2]
        t0 <- mean(years[okmean])

        # Centered trend: removing it preserves the mean level
        common_trend <- slope * (years - t0)

        # Subtract exactly the same trend from every member
        for (m in seq_len(dim(x)[4])) {

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


summary_land <- function(r, sig, mask) {

  valid <- mask & is.finite(r)

  c(
    mean_r =
      mean(r[valid], na.rm = TRUE),

    median_r =
      median(r[valid], na.rm = TRUE),

    fraction_positive =
      mean(r[valid] > 0, na.rm = TRUE),

    fraction_significant =
      mean(sig[valid], na.rm = TRUE)
  )
}

# ============================================================
# Land mask — exactly from FWI95d OBS verification field
# ============================================================

fwi_obs_file <- file.path(
  verification_dir,
  "FWI95d_OBS_runmean5_1963-2022_corrected.nc"
)

stopifnot(file.exists(fwi_obs_file))

fwi <- read_nc_field(
  fwi_obs_file,
  preferred_var = "FWI"
)

fwi_mean <- apply(
  fwi$data,
  c(1, 2),
  mean,
  na.rm = TRUE
)

land_mask <- is.finite(fwi_mean)

cat(
  "\nLand grid cells:",
  sum(land_mask),
  "\n"
)

# ============================================================
# Main loop
# ============================================================

results <- list()
summary_rows <- list()

for (var in variables) {

  cat("\n")
  cat("====================================================\n")
  cat("VARIABLE:", var, "\n")
  cat("====================================================\n")

  # ----------------------------------------------------------
  # OBS
  # ----------------------------------------------------------

  obs_file <- file.path(
    driver_dir,
    "OBS",
    var,
    paste0(
      var,
      "_OBS_FY1-FY5_1963-2022.nc"
    )
  )

  stopifnot(file.exists(obs_file))

  oo <- read_nc_field(
    obs_file,
    preferred_var = var
  )

  lon <- oo$lon
  lat <- oo$lat
  obs <- oo$data

  stopifnot(
    all(
      dim(obs) ==
        c(
          length(lon),
          length(lat),
          ntime
        )
    )
  )

  cat(
    "OBS :",
    paste(dim(obs), collapse = " x "),
    "\n"
  )

  # ----------------------------------------------------------
  # INIT
  # ----------------------------------------------------------

  init_files <- list.files(
    file.path(
      driver_dir,
      "INIT",
      var
    ),
    pattern = paste0(
      "^",
      var,
      "_INIT_FY1-FY5_r[0-9]+i1p1f1_1963-2022\\.nc$"
    ),
    full.names = TRUE
  )

  init_files <- sort(init_files)

  if (length(init_files) != 24)
    stop(
      var,
      ": expected 24 INIT files, found ",
      length(init_files)
    )

  init <- array(
    NA_real_,
    dim = c(
      length(lon),
      length(lat),
      ntime,
      length(init_files)
    )
  )

  init_members <- character(
    length(init_files)
  )

  for (i in seq_along(init_files)) {

    xx <- read_nc_field(
      init_files[i],
      preferred_var = var
    )

    if (!all(
      dim(xx$data) ==
        c(
          length(lon),
          length(lat),
          ntime
        )
    ))
      stop(
        "INIT dimensions mismatch: ",
        init_files[i]
      )

    init[, , , i] <- xx$data

    init_members[i] <- sub(
      paste0(
        var,
        "_INIT_FY1-FY5_(r[0-9]+i1p1f1)_.*"
      ),
      "\\1",
      basename(init_files[i])
    )
  }

  cat(
    "INIT:",
    paste(dim(init), collapse = " x "),
    "\n"
  )

  # ----------------------------------------------------------
  # HIST
  # ----------------------------------------------------------

  hist_files <- list.files(
    file.path(
      driver_dir,
      "HIST",
      var
    ),
    pattern = paste0(
      "^",
      var,
      "_HIST_FY1-FY5_r[0-9]+i1p2f1_1963-2022\\.nc$"
    ),
    full.names = TRUE
  )

  hist_files <- sort(hist_files)

  if (length(hist_files) != 10)
    stop(
      var,
      ": expected 10 HIST files, found ",
      length(hist_files)
    )

  hist <- array(
    NA_real_,
    dim = c(
      length(lon),
      length(lat),
      ntime,
      length(hist_files)
    )
  )

  hist_members <- character(
    length(hist_files)
  )

  for (i in seq_along(hist_files)) {

    xx <- read_nc_field(
      hist_files[i],
      preferred_var = var
    )

    if (!all(
      dim(xx$data) ==
        c(
          length(lon),
          length(lat),
          ntime
        )
    ))
      stop(
        "HIST dimensions mismatch: ",
        hist_files[i]
      )

    hist[, , , i] <- xx$data

    hist_members[i] <- sub(
      paste0(
        var,
        "_HIST_FY1-FY5_(r[0-9]+i1p2f1)_.*"
      ),
      "\\1",
      basename(hist_files[i])
    )
  }

  cat(
    "HIST:",
    paste(dim(hist), collapse = " x "),
    "\n"
  )

  # ==========================================================
  # Detrend
  # ==========================================================

  cat("\nDetrending OBS...\n")
  obs_dt <- detrend_3d(obs)

  cat("Detrending INIT...\n")
  init_dt <- detrend_4d(init)

  cat("Detrending HIST...\n")
  hist_dt <- detrend_4d(hist)

  # ==========================================================
  # Add named dimensions required by CorrEno
  # ==========================================================

  dimnames(obs_dt) <- list(
    lon = lon,
    lat = lat,
    sdate = as.character(years)
  )

  dn <- dim(obs_dt)
  names(dn) <- c(
    "lon",
    "lat",
    "sdate"
  )
  dim(obs_dt) <- dn


  dimnames(init_dt) <- list(
    lon = lon,
    lat = lat,
    sdate = as.character(years),
    member = init_members
  )

  dn <- dim(init_dt)
  names(dn) <- c(
    "lon",
    "lat",
    "sdate",
    "member"
  )
  dim(init_dt) <- dn


  dimnames(hist_dt) <- list(
    lon = lon,
    lat = lat,
    sdate = as.character(years),
    member = hist_members
  )

  dn <- dim(hist_dt)
  names(dn) <- c(
    "lon",
    "lat",
    "sdate",
    "member"
  )
  dim(hist_dt) <- dn

  # ==========================================================
  # Correlation — same method as raw driver analysis
  # ==========================================================

  cat("\nCorrelation INIT vs OBS...\n")

  corr_init <- CorrEno(
    exp = init_dt,
    obs = obs_dt,
    time_dim = "sdate",
    member_dim = "member",
    method = "pearson",
    alpha = 0.05,
    test.type = "two-sided",
    pval = TRUE,
    handle.na = "return.na",
    ncores = 8
  )

  cat("Correlation HIST vs OBS...\n")

  corr_hist <- CorrEno(
    exp = hist_dt,
    obs = obs_dt,
    time_dim = "sdate",
    member_dim = "member",
    method = "pearson",
    alpha = 0.05,
    test.type = "two-sided",
    pval = TRUE,
    handle.na = "return.na",
    ncores = 8
  )

  r_init <- corr_init$r
  r_hist <- corr_hist$r

  sig_init <- corr_init$sign
  sig_hist <- corr_hist$sign

  # ==========================================================
  # Summary over same land mask
  # ==========================================================

  init_summary <-
    summary_land(
      r_init,
      sig_init,
      land_mask
    )

  hist_summary <-
    summary_land(
      r_hist,
      sig_hist,
      land_mask
    )

  summary_rows[[var]] <- data.frame(
    variable = var,

    INIT_mean_r =
      init_summary["mean_r"],

    INIT_median_r =
      init_summary["median_r"],

    INIT_fraction_positive =
      init_summary["fraction_positive"],

    INIT_fraction_significant =
      init_summary["fraction_significant"],

    HIST_mean_r =
      hist_summary["mean_r"],

    HIST_median_r =
      hist_summary["median_r"],

    HIST_fraction_positive =
      hist_summary["fraction_positive"],

    HIST_fraction_significant =
      hist_summary["fraction_significant"]
  )

  cat("\nDETREND SUMMARY\n")
  cat(
    "INIT mean r =",
    round(init_summary["mean_r"], 3),
    " | positive =",
    round(init_summary["fraction_positive"], 3),
    " | significant =",
    round(init_summary["fraction_significant"], 3),
    "\n"
  )

  cat(
    "HIST mean r =",
    round(hist_summary["mean_r"], 3),
    " | positive =",
    round(hist_summary["fraction_positive"], 3),
    " | significant =",
    round(hist_summary["fraction_significant"], 3),
    "\n"
  )

  # ==========================================================
  # Save variable result
  # ==========================================================

  results[[var]] <- list(
    lon = lon,
    lat = lat,
    years = years,

    corr_init = corr_init,
    corr_hist = corr_hist,

    r_init_land =
      ifelse(
        land_mask,
        r_init,
        NA_real_
      ),

    r_hist_land =
      ifelse(
        land_mask,
        r_hist,
        NA_real_
      ),

    sig_init_land =
      ifelse(
        land_mask,
        sig_init,
        NA
      ),

    sig_hist_land =
      ifelse(
        land_mask,
        sig_hist,
        NA
      )
  )
}

# ============================================================
# Final summary
# ============================================================

summary_table <- do.call(
  rbind,
  summary_rows
)

rownames(summary_table) <- NULL

cat("\n")
cat("============================================================\n")
cat("DETREND DRIVER SKILL — FINAL SUMMARY\n")
cat("============================================================\n")

print(
  summary_table,
  row.names = FALSE
)

# ============================================================
# Compare directly against RAW results
# ============================================================

raw_file <- file.path(
  out_dir,
  "driver_skill_FY1-FY5_corrected.rds"
)

if (file.exists(raw_file)) {

  raw <- readRDS(raw_file)

  cat("\n")
  cat("============================================================\n")
  cat("RAW -> DETRENDED INIT MEAN CORRELATION\n")
  cat("============================================================\n")

  for (var in variables) {

    raw_mean <- mean(
      raw[[var]]$r_init_land,
      na.rm = TRUE
    )

    dt_mean <- mean(
      results[[var]]$r_init_land,
      na.rm = TRUE
    )

    cat(
      sprintf(
        "%-8s  %+.3f  ->  %+.3f\n",
        var,
        raw_mean,
        dt_mean
      )
    )
  }
}

# ============================================================
# Save
# ============================================================

out_rds <- file.path(
  out_dir,
  "driver_skill_detrended_commontrend_FY1-FY5_corrected.rds"
)

out_csv <- file.path(
  out_dir,
  "driver_skill_detrended_commontrend_summary_corrected.csv"
)

saveRDS(
  results,
  out_rds
)

write.csv(
  summary_table,
  out_csv,
  row.names = FALSE
)

cat("\nSaved:\n")
cat(out_rds, "\n")
cat(out_csv, "\n")
cat("\nDONE\n")

