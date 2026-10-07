#!/bin/bash 

# Save original working directory
orig_dir=$(pwd)

cp ~/.cdsapirc_cds ~/.cdsapirc

y1=1976
y2=2024
m1=1
m2=12
d1=1
d2=31

# Land-sea mask file for regridding
LSM_FILE="/diskonfire/ERA5/daily_values/grids/land_sea_mask_1degree.nc4"

dirwrk="/home/marco/Dropbox/estcena/scripts/ONFIRE/decadal"

# Define output directory for processed files
output_dir="/diskonfire/ERA5/daily_values/tasmean"

# Variables to be downloaded
vars=("tas")
varlongnames=('2m_temperature')

for year in $(seq $y1 $y2); do
    for month in $(seq -w $m1 $m2); do
        for day in $(seq -w $d1 $d2); do
            # Convert month and day removing leading zeros then reformat them
            month=$(printf "%02d" $(echo $month | sed 's/^0*//'))
            day=$(printf "%02d" $(echo $day | sed 's/^0*//'))
            
            # Check if the date is valid (e.g., skips 1963-04-31)
            if ! date -d "${year}-${month}-${day}" >/dev/null 2>&1; then
                echo "Invalid date: ${year}-${month}-${day}. Skipping."
                continue
            fi

            for index in ${!vars[@]}; do
                var=${vars[$index]}
                varlongname=${varlongnames[$index]}

                # Final output filename
                final_nc_file="${output_dir}/tasmean_${year}${month}${day}_1degree.nc"

                # Check if the file already exists
                if [ -f "$final_nc_file" ]; then
                    echo "File already exists: $final_nc_file. Skipping download."
                    continue
                fi

                # Return to original working directory before starting a new iteration
                cd "$orig_dir"

                # Create a temporary directory for the download
                temp_dir=$(mktemp -d)
                echo "Temporary directory: $temp_dir"

                # Create the sed file in the temporary directory
                sed_file="${temp_dir}/kk4sed"
                cat > "$sed_file" <<EOF
s#SSSVARIABLE#${varlongname}#g
s#SSSYEAR#${year}#g
s#SSSMONTH#${month}#g
s#SSSDAY#${day}#g
EOF

                # Ensure the download script template exists
                if [ ! -f "$dirwrk/download-cds-era5-4sed.py" ]; then
                    echo "Error: download-cds-era5-4sed.py not found."
                    exit 1
                fi

                # Apply the sed substitutions and create the Python download script in the temp directory
                sed -f "$sed_file" "$dirwrk/download-cds-era5-4sed.py" > "$temp_dir/download_cds_era5_${var}.py"
                rm "$sed_file"

                chmod u+x "$temp_dir/download_cds_era5_${var}.py"

                # Change into the temporary directory to run the Python script
                cd "$temp_dir"

                # Run the Python script
                python3 "download_cds_era5_${var}.py"

                # Find the most recently downloaded .grib file
                latest_file=$(find . -name "*.grib" -type f | tail -n 1)
                if [ -z "$latest_file" ]; then
                    echo "No GRIB file found for ${var} on ${year}-${month}-${day}."
                    rm -rf "$temp_dir"
                    continue
                fi
                echo "Processing file: $latest_file"

                # Convert GRIB to NetCDF
                cdo -f nc copy "$latest_file" "${var}_${year}${month}${day}.nc"

                # Rename variable in NetCDF file
                cdo chname,2t,tas "${var}_${year}${month}${day}.nc" ofile.nc
                mv ofile.nc "${var}_${year}${month}${day}.nc"

                # Calculate daily mean of TAS
                cdo -b F32 -daymean "${var}_${year}${month}${day}.nc" "tasmean_${year}${month}${day}.nc"

                # Regrid the daily mean TAS to 1-degree resolution 
                cdo -b F32 remapbil,"$LSM_FILE" "tasmean_${year}${month}${day}.nc" "tasmean_${year}${month}${day}_1degree.nc"

                # Return to original directory before moving the final file
                cd "$orig_dir"

                # Move the final file to output directory
                mv "$temp_dir/tasmean_${year}${month}${day}_1degree.nc" "$final_nc_file"

                # Clean up temporary files and directory
                rm -rf "$temp_dir"
            done
        done
    done
done

