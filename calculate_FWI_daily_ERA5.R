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
library(maps)

# PARAMETER SETTING FOR DATA LOADING, INDEX CALCULATION, AND EXPORT--------------------------------
dir_proxy <- "/diskonfire/ERA5/daily_values/"
dir_grid <- "/diskonfire/Decadal/"

# Load land-sea mask
print("Loading land-sea mask")
sftlf <-
  loadGridData(dataset = paste0(dir_grid, "landsea_mask.nc"),
               "FWI",
               dictionary = FALSE)
refGrid <- getGrid(sftlf)

dim(sftlf)
dim(refGrid)

# Define corresponding files for hurs, pr, and sfcWind
tas_file <-
  paste0(dir_proxy, "tasmean/tasmean_remapped.nc")
hurs_file <-
  paste0(dir_proxy, "hurs/rhmean_remapped.nc")
pr_file <- paste0(dir_proxy, "pr/precip_remapped.nc")
wind_file <-

# Load data
tas.y <-
  loadGridData(dataset = tas_file,
               var = "tas",
               dictionary = FALSE)
hurs.y <-
  loadGridData(dataset = hurs_file,
               var = "rhmean",
               dictionary = FALSE)
pr.y <-
  loadGridData(dataset = pr_file,
               var = "tp",
               dictionary = FALSE)
wss.y <-
  loadGridData(dataset = wind_file,
               var = "10u",
               dictionary = FALSE)


tas_name = "tt"
hurs_name = "rh"
pr_name = "pr"
ws_name = "ws"



# Unit conversions for tasmean/tasmax and wind
tas.y$Data <- tas.y$Data - 273.15  # Kelvin to Celsius
attr(tas.y$Variable, "units") <- "celsius"

wss.y$Data <- wss.y$Data * 3.6  # m/s to km/h
attr(wss.y$Variable, "units") <- "km/h"


pr.y$Data <-
  pr.y$Data * 1000  # m to mm only when pr is the proxy
attr(pr.y$Variable, "units") <- "mm"


tas.y$Variable$varName <- "tas"
pr.y$Variable$varName <- "tp"
hurs.y$Variable$varName <- "hurs"
wss.y$Variable$varName <- "wss"

check=sftlf$Data
#  library(fields)
lon = pr.y$xyCoords$x
lat = pr.y$xyCoords$y

# 
#  # check <- t(apply(wss.y$Data,c(2,3),mean))
# # check <- t(apply(tas.y$Data,c(2,3),mean))
# check <- t(apply(pr.y$Data,c(2,3),mean))
# # 
# check <- t(apply(hurs.y$Data,c(2,3),mean))
#  
image.plot(lon,lat,t(check))
world(add = TRUE, col = "black")

# Convertir las fechas a objetos Date (por si acaso no lo son ya)
# Asegurarse de que las fechas estén en formato Date
tas_dates <- as.Date(tas.y$Dates$start)
hurs_dates <- as.Date(hurs.y$Dates$start)
pr_dates <- as.Date(pr.y$Dates$start)
wss_dates <- as.Date(wss.y$Dates$start)

# Extraer el rango (mínimo y máximo) de fechas para cada variable
tas_range <- range(tas_dates)
hurs_range <- range(hurs_dates)
pr_range <- range(pr_dates)
wss_range <- range(wss_dates)

# Crear secuencias de fechas a partir de los rangos
dates_tas <- seq(from = tas_range[1], to = tas_range[2], by = "day")
dates_hurs <- seq(from = hurs_range[1], to = hurs_range[2], by = "day")
dates_pr <- seq(from = pr_range[1], to = pr_range[2], by = "day")
dates_wss <- seq(from = wss_range[1], to = wss_range[2], by = "day")

print(paste("Número de días en tas:", length(dates_tas)))
print(paste("Número de días en hurs:", length(dates_hurs)))
print(paste("Número de días en pr:", length(dates_pr)))
print(paste("Número de días en wss:", length(dates_wss)))




# Convert dates to Date objects
tas.y$Dates$start <- as.Date(tas.y$Dates$start)
tas.y$Dates$end <- as.Date(tas.y$Dates$end)
hurs.y$Dates$start <- as.Date(hurs.y$Dates$start)
hurs.y$Dates$end <- as.Date(hurs.y$Dates$end)
pr.y$Dates$start <- as.Date(pr.y$Dates$start)
pr.y$Dates$end <- as.Date(pr.y$Dates$end)
wss.y$Dates$start <- as.Date(wss.y$Dates$start)
wss.y$Dates$end <- as.Date(wss.y$Dates$end)





# Calculate FWI
cat("Calculating FWI")
multigrid_fwi <-
  makeMultiGrid(
    tasmax = tas.y,
    hurs = hurs.y,
    pr = pr.y,
    wind = wss.y,
    skip.temporal.check = TRUE
  )
fwi.y <-
  fwiGrid(multigrid = multigrid_fwi,
          mask = sftlf,
          restart.annual = FALSE)
print("FWI calculation end")

fwi.y$Data <-
  fwi.y$Data[1, , , ]  # This will index out the first dimension
dim(fwi.y$Data)

print("Creating output file")
# Export FWI to NetCDF
ncFile <-
  paste0(dir_proxy, "FWI_1961_2024_1degree_europe.nc")
cat("Saving FWI to", ncFile, "\n")
grid2nc(
  data = fwi.y,
  NetCDFOutFile = ncFile,
  missval = 1e20,
  prec = "float",
  globalAttributes = list()
)

# Clean up
rm(list = c("tas.y", "hurs.y", "pr.y", "wss.y", "fwi.y"))
gc()

