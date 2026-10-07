#!/bin/bash
# remap_data.sh

# Directories (adjust if necessary)
dir_proxy="/diskonfire/ERA5/daily_values/"
dir_grid="/diskonfire/Decadal/"

# Define the land‐sea mask file (target grid)
mask_file="${dir_grid}landsea_mask.nc"

# Input files
tas_file="${dir_proxy}tasmean/tasmean_1961_2024_1degree_europe.nc"
hurs_file="${dir_proxy}hurs/rhmean_1961_2024_1degree_europe.nc"
pr_file="${dir_proxy}pr/precip_1961_2024_1degree_europe.nc"
wind_file="${dir_proxy}sfcWind/sfcWind_1961_2024_1degree_europe.nc"

# Output files (you may change the names as required)
tas_out="${dir_proxy}tasmean/tasmean_remapped.nc"
hurs_out="${dir_proxy}hurs/rhmean_remapped.nc"
pr_out="${dir_proxy}pr/precip_remapped.nc"
wind_out="${dir_proxy}sfcWind/sfcWind_remapped.nc"

# echo "Remapping tasmean using bilinear interpolation..."
# cdo remapbil,$mask_file $tas_file $tas_out
# if [ $? -ne 0 ]; then
#     echo "Error remapping tasmean."
#     exit 1
# fi

# echo "Remapping rhmean using bilinear interpolation..."
# cdo remapbil,$mask_file $hurs_file $hurs_out
# if [ $? -ne 0 ]; then
#     echo "Error remapping rhmean."
#     exit 1
# fi

# echo "Remapping precip using conservative remapping..."
# cdo remapcon,$mask_file $pr_file $pr_out
# if [ $? -ne 0 ]; then
#     echo "Error remapping precip."
#     exit 1
# fi

# echo "Remapping sfcWind using bilinear interpolation..."
# cdo remapbil,$mask_file $wind_file $wind_out
# if [ $? -ne 0 ]; then
#     echo "Error remapping sfcWind."
#     exit 1
# fi

cdo -O -b F32 remapbil,/diskonfire/Decadal/landsea_mask.nc \
    -sellonlatbox,-10,40,30,72 \
    -mergetime /diskonfire/ERA5/daily_values/sfcWind_v2/yearly/sfcWind_daily_global_*.nc \
    /diskonfire/ERA5/daily_values/sfcWind_v2/sfcWind_remapped.nc

echo "Remapping completed successfully."
