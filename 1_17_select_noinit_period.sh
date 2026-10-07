#!/bin/bash

# Directorio donde están los archivos
DIR="/diskonfire/Decadal/processed_NO-INIT"
cd "$DIR"

# Procesar todos los archivos *_regridded.nc
for file in *_regridded.nc; do
    echo "Recortando años 1960–2024 en $file"
    
    # Guardar en archivo temporal
    tmpfile="tmp_${file}"
    cdo selyear,1960/2024 "$file" "$tmpfile"
    
    # Sobrescribir archivo original
    mv "$tmpfile" "$file"
done

echo "¡Todos los archivos han sido recortados y sobrescritos!"
