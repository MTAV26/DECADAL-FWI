rm(list = ls())
graphics.off()
gc()

library(loadeR)
library(transformeR)
library(visualizeR)
library(loadeR.2nc)
library(fireDanger)
library(fields)

# Directorios
dir_model <- "/diskonfire/Decadal/processed_NO-INIT/"
dir_grid <- "/diskonfire/Decadal/"

# Cargar máscara de tierra
cat("Loading land-sea mask...\n")
sftlf <- loadGridData(dataset = paste0(dir_grid, "landsea_mask.nc"),
                      var = "FWI", dictionary = FALSE)
refGrid <- getGrid(sftlf)

# Lista de miembros
members <- paste0("r", 2:11, "i1p2f1")

for (member in members) {
  cat("\n==== Procesando miembro:", member, "====\n")
  
  # Rutas a los archivos
  tas_file <- paste0(dir_model, "tas_merged_CMCC-CM2-SR5_", member, "_19600101-20241231_regridded.nc")
  hurs_file <- paste0(dir_model, "hurs_merged_CMCC-CM2-SR5_", member, "_19600101-20241231_regridded.nc")
  pr_file <- paste0(dir_model, "pr_merged_CMCC-CM2-SR5_", member, "_19600101-20241231_regridded.nc")
  wind_file <- paste0(dir_model, "sfcWind_merged_CMCC-CM2-SR5_", member, "_19600101-20241231_regridded.nc")
  
  # Cargar datos
  tas.y <- loadGridData(tas_file, var = "tas", dictionary = FALSE)
  hurs.y <- loadGridData(hurs_file, var = "hurs", dictionary = FALSE)
  pr.y <- loadGridData(pr_file, var = "pr", dictionary = FALSE)
  wss.y <- loadGridData(wind_file, var = "sfcWind", dictionary = FALSE)
  
  # Conversiones de unidades
  tas.y$Data <- tas.y$Data - 273.15
  attr(tas.y$Variable, "units") <- "celsius"
  
  wss.y$Data <- wss.y$Data * 3.6
  attr(wss.y$Variable, "units") <- "km/h"
  
  pr.y$Data <- pr.y$Data * 86400
  attr(pr.y$Variable, "units") <- "mm"
  
  # Renombrar variables
  tas.y$Variable$varName <- "tas"
  pr.y$Variable$varName <- "tp"
  hurs.y$Variable$varName <- "hurs"
  wss.y$Variable$varName <- "wss"
  
  # Fechas como objetos Date
  for (var in list(tas.y, hurs.y, pr.y, wss.y)) {
    var$Dates$start <- as.Date(var$Dates$start)
    var$Dates$end <- as.Date(var$Dates$end)
  }
  
  # Calcular FWI
  cat("Calculando FWI para ", member, "\n")
  multigrid_fwi <- makeMultiGrid(
    tasmax = tas.y,
    hurs = hurs.y,
    pr = pr.y,
    wind = wss.y,
    skip.temporal.check = TRUE
  )
  
  fwi.y <- fwiGrid(multigrid = multigrid_fwi,
                   mask = sftlf,
                   restart.annual = FALSE)
  
  # Eliminar primera dimensión
  fwi.y$Data <- fwi.y$Data[1, , , ]
  
  # Guardar resultado
  ncFile <- paste0(dir_model, "FWI_CMCC-CM2-SR5_", member, "_1960_2024.nc")
  cat("Guardando FWI en", ncFile, "\n")
  grid2nc(data = fwi.y,
          NetCDFOutFile = ncFile,
          missval = 1e20,
          prec = "float",
          globalAttributes = list())
  
  # Limpieza
  rm(tas.y, hurs.y, pr.y, wss.y, fwi.y)
  gc()
}

