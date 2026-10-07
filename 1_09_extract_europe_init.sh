#!/bin/bash

set -euo pipefail

BASE_DIR="/diskonfire/Decadal"
OUTPUT_SUFFIX="_europe.nc"
REGION="-10,40,30,72"

cut_files_in_folder() {
    local var_dir=$1

    cd "${BASE_DIR}/${var_dir}" || exit 1
    echo "Procesando archivos en: ${var_dir}"

    for file in *.nc; do

        # No volver a procesar archivos europeos
        if [[ "$file" == *_europe.nc ]]; then
            continue
        fi

        output_file="${file%.nc}${OUTPUT_SUFFIX}"

        if [[ ! -f "$output_file" ]]; then
            echo "Cortando: $file"
            cdo sellonlatbox,${REGION} "$file" "$output_file"
        else
            echo "Ya existe, omitiendo: $output_file"
        fi
    done
}

cut_files_in_folder "pr_day"
cut_files_in_folder "tas_day"
cut_files_in_folder "hurs_day"
cut_files_in_folder "sfcWind_day"

echo "Proceso completado."