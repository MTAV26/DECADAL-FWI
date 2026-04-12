# Clear environment and free memory
rm(list = ls())
graphics.off()
gc()

library(loadeR)
library(transformeR)
library(visualizeR)
library(loadeR.2nc)
library(fireDanger)
library(RColorBrewer)
library(fields)

# Define file paths
dir_data <- "/diskonfire/Decadal/"
output_dir <- "/diskonfire/Decadal/FWI/"

# Load land-sea mask
print("Loading land-sea mask")
sftlf <- loadGridData(dataset = paste0(dir_data, "landsea_mask.nc"), "FWI", dictionary = FALSE)

# List all tas files and extract the matching parts
tas_files <- list.files(path = paste0(dir_data, "tas_day"), pattern = "*_europe.nc", full.names = TRUE)

# Initialize a log file for missing files
log_file <- paste0(output_dir, "missing_files.log")
write("Missing files log:", file = log_file)

# Loop over each tas file and process corresponding variables
for (tas_file in tas_files) {
  
  # Extract the common identifier for the current run (e.g., 's2019-r40i1p1f1')
  run_id <- sub("tas_day_|_europe.nc", "", basename(tas_file))
  
  # Define corresponding files for hurs, pr, and sfcWind
  hurs_file <- paste0(dir_data, "hurs_day/hurs_day_", run_id)
  pr_file <- paste0(dir_data, "pr_day/pr_day_", run_id)
  wind_file <- paste0(dir_data, "sfcWind_day/sfcWind_day_", run_id)
  
  # Check if each file exists and log missing files separately
  missing <- FALSE
  
  if (!file.exists(hurs_file)) {
    message("Missing hurs file for run: ", hurs_file)
    write(paste("Missing hurs file for run:", hurs_file), file = log_file, append = TRUE)
    missing <- TRUE
  }
  
  if (!file.exists(pr_file)) {
    message("Missing pr file for run: ", pr_file)
    write(paste("Missing pr file for run:", pr_file), file = log_file, append = TRUE)
    missing <- TRUE
  }
  
  if (!file.exists(wind_file)) {
    message("Missing sfcWind file for run: ", wind_file)
    write(paste("Missing sfcWind file for run:", wind_file), file = log_file, append = TRUE)
    missing <- TRUE
  }
  
  # Skip the current run if any of the files are missing
  if (missing) {
    next
  }
  
  # Inside the loop (after checking if files exist)
  if (!file.exists(hurs_file) || !file.exists(pr_file) || !file.exists(wind_file)) {
    message("Missing files for run: ", run_id, " - Skipping this run.")
    write(paste("Missing files for run:", run_id), file = log_file, append = TRUE)
    next
  }
  
  print(paste("Processing run:", run_id))
  
  # Load data
  tas.y <- loadGridData(dataset = tas_file, var = "tas", dictionary = FALSE)
  hurs.y <- loadGridData(dataset = hurs_file, var = "hurs", dictionary = FALSE)
  pr.y <- loadGridData(dataset = pr_file, var = "pr", dictionary = FALSE)
  wss.y <- loadGridData(dataset = wind_file, var = "sfcWind", dictionary = FALSE)
  
  # Convert units
  pr.y$Data <- pr.y$Data * 86400
  attr(pr.y$Variable, "units") <- "mm/day"
  
  tas.y$Data <- tas.y$Data - 273.15
  attr(tas.y$Variable, "units") <- "celsius"
  
  wss.y$Data <- wss.y$Data * 3.6
  attr(wss.y$Variable, "units") <- "km/h"
  
  # Set date formatting
  pr.y$Dates$start <- as.Date(pr.y$Dates$start)
  pr.y$Dates$end <- as.Date(pr.y$Dates$end)
  tas.y$Dates$start <- as.Date(tas.y$Dates$start)
  tas.y$Dates$end <- as.Date(tas.y$Dates$end)
  hurs.y$Dates$start <- as.Date(hurs.y$Dates$start)
  hurs.y$Dates$end <- as.Date(hurs.y$Dates$end)
  wss.y$Dates$start <- as.Date(wss.y$Dates$start)
  wss.y$Dates$end <- as.Date(wss.y$Dates$end)
  
  # Set variable names
  tas.y$Variable$varName <- "tas"
  pr.y$Variable$varName <- "tp"
  hurs.y$Variable$varName <- "hurs"
  wss.y$Variable$varName <- "wss"
  
  # Calculate FWI
  print("FWI calculation start")
  multigrid_fwi <- makeMultiGrid(tasmax = tas.y, hurs = hurs.y, pr = pr.y, wind = wss.y, skip.temporal.check = TRUE)
  fwi.y <- fwiGrid(multigrid = multigrid_fwi, mask = sftlf, restart.annual = FALSE)
  print("FWI calculation end")
  
  # Extract the first time slice
  fwi.y$Data <- fwi.y$Data[1,,,]
  
  # Define output file name and export to NetCDF
  ncFile <- paste0(output_dir, "FWI_", run_id)
  grid2nc(data = fwi.y, NetCDFOutFile = ncFile, missval = 1e20, prec = "float", globalAttributes = list())
  
  print(paste("FWI saved to:", ncFile))
  
  # Clean up
  rm(list = c("pr.y", "hurs.y", "tas.y", "wss.y", "fwi.y"))
  gc()
}

print("All runs processed successfully.")
