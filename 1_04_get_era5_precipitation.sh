#!/bin/bash

cp ~/.cdsapirc_cds ~/.cdsapirc

# Land-sea mask file for regridding
LSM_FILE="/diskonfire/ERA5/daily_values/grids/land_sea_mask_1degree.nc4"

dirwrk="/home/marco/Dropbox/estcena/scripts/ONFIRE/decadal"
dirout="/diskonfire/ERA5/daily_values/pr"

# Function to create the sed file and download GRIB for a given date.
sed_process_and_download() {
    local y=$1
    local m=$2
    local d=$3
    local file_var_name=$4

    # Create the sed file with substitutions.
    cat > "$dirwrk/kk4sed" <<EOF
s#SSSVARIABLE#${varlongname}#g
s#SSSYEAR#${y}#g
s#SSSMONTH#${m}#g
s#SSSDAY#${d}#g
EOF

    # Create a temporary directory for the download.
    temp_dir=$(mktemp -d)
    echo "Temporary directory: $temp_dir"

    cd "$temp_dir"
    sed -f "$dirwrk/kk4sed" "$dirwrk/download-cds-era5-4sed.py" > "$temp_dir/download_cds_era5_${var}.py"
    chmod u+x "$temp_dir/download_cds_era5_${var}.py"

    # Run the Python download script.
    python3 "$temp_dir/download_cds_era5_${var}.py"

    # Find the most recently downloaded .grib file.
    latest_file=$(find "$temp_dir" -name "*.grib" -type f | tail -n 1)
    if [ -z "$latest_file" ]; then
        echo "No GRIB file found for ${var} on ${y}-${m}-${d}."
        rm -rf "$temp_dir"
        return
    fi

    local file_name="${dirout}/precip_${y}${m}${d}_download.grib"
    mv "$latest_file" "$file_name"
    chmod u+x "$file_name"
    # Safely assign the file name to the variable passed as argument.
    printf -v "$file_var_name" '%s' "$file_name"
    echo "Set variable '$file_var_name' to '$file_name'"

    rm -rf "$temp_dir"
}

vars='tp'
y1=1961
y2=2024
m1=1
m2=12
d2=31

for variable in $vars; do
    case $variable in
        'tp') varlongname='total_precipitation' ;;
    esac

    next_file=""

    for year in $(seq $y1 $y2); do
        for month in $(seq -w $m1 $m2); do
            # month is already two-digit from seq -w.
            for d_val in $(seq -w 1 $d2); do

                echo "Processing date: ${year}-${month}-${d_val}"

                # Validate that the date is valid.
                if ! date -d "${year}-${month}-${d_val}" >/dev/null 2>&1; then
                    echo "Invalid date: ${year}-${month}-${d_val}. Skipping."
                    continue
                fi

                # Convert month and day to numbers to avoid octal issues.
                month_num=$((10#$month))
                d_val_num=$((10#$d_val))
                final_file="${dirout}/precip_$(printf "%04d" $year)$(printf "%02d" $month_num)$(printf "%02d" $d_val_num)_1degree.nc"
                echo "Final file will be: $final_file"
                if [ -f "$final_file" ]; then
                    echo "File already exists: $final_file. Skipping precipitation calculation for ${year}-${month}-${d_val}."
                    continue
                fi

                # Reset window variables for each iteration.
                curr_file=""
                next_file=""

                # Download for the current day.
                sed_process_and_download $year $month $d_val "curr_file"

                # Compute the next day using ISO format.
                next_day=$(date -I -d "$year-$month-$d_val + 1 day")
                echo "Next day computed: $next_day"
                sed_process_and_download ${next_day:0:4} ${next_day:5:2} ${next_day:8:2} "next_file"

                # If both files are available, combine them.
                if [ -f "$curr_file" ] && [ -f "$next_file" ]; then
                    cdo -b F32 timselsum,23,1 "$curr_file" "$dirout/tmp1.nc"
                    cdo -b F32 selhour,0 "$next_file" "$dirout/tmp2.nc"
                    cdo -b F32 add "$dirout/tmp1.nc" "$dirout/tmp2.nc" "$dirout/ofile.nc"

                    cdo -b F32 -f nc remapcon,$LSM_FILE "$dirout/ofile.nc" "$final_file"

                    rm "$dirout/tmp1.nc" "$dirout/tmp2.nc" "$curr_file" "$dirout/ofile.nc"
                else
                    echo "Missing files for ${year}-${month}-${d_val}, skipping precipitation calculation."
                fi

            done
        done
    done
done
