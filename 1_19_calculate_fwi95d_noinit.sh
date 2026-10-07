#!/bin/bash

# Directorio donde están los FWI de cada miembro
DIR="/diskonfire/Decadal/processed_NO-INIT"
cd "$DIR"

# Baseline para percentil 95
BASELINE_START=1991
BASELINE_END=2020

# Loop sobre cada archivo de FWI
for file in FWI_CMCC-CM2-SR5_r*i1p2f1_1960_2024.nc; do
  echo "Procesando $file"

  base=$(basename "$file" .nc)

  # Paso 1: extraer periodo 1991-2020 para baseline
  cdo selyear,${BASELINE_START}/${BASELINE_END} "$file" ${base}_baseline.nc

  # Paso 2: calcular percentil 95 (FWI95)
  cdo timpctl,95 ${base}_baseline.nc -timmin ${base}_baseline.nc -timmax ${base}_baseline.nc ${base}_pctl95.nc

  # Paso 3: contar días por año en que FWI > FWI95
  cdo -yearsum -gt "$file" ${base}_pctl95.nc ${base}_FWI95d.nc

  # Limpieza intermedia (opcional)
  rm ${base}_baseline.nc ${base}_pctl95.nc ${base}_pctl95_const.nc

  echo "→ Guardado: ${base}_FWI95d.nc"
done

echo "✅ ¡Proceso completo para todos los miembros!"
