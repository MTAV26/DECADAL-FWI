#!/bin/bash

set -euo pipefail

indir="/diskonfire/Decadal/processed_INIT"
outdir="${indir}/FY1-FY5_mean"

mkdir -p "$outdir"

members=$(seq -f "r%gi1p1f1" 17 40)

for member in $members; do

    echo "=== $member ==="

    tmpdir=$(mktemp -d)

    # 60 initializations: 1960-2019
    for start in $(seq 1960 2019); do

        echo "  Initialization $start"

        selected=()

        # FY1 ... FY5
        for f in 1 2 3 4 5; do

            target_year=$((start + f))

            # Existing FY series:
            # fy1 = 1961-2020
            # fy2 = 1962-2021
            # ...
            # fy5 = 1965-2024
            first_year=$((1960 + f))
            last_year=$((2019 + f))

            infile="${indir}/FWI_${member}_fy${f}_FWI95d_${first_year}-${last_year}.nc"

            tmp="${tmpdir}/${member}_init${start}_fy${f}.nc"

            if [[ ! -f "$infile" ]]; then
                echo "ERROR: missing $infile"
                exit 1
            fi

            # Select the appropriate forecast year.
            # Standardize timestamp to the initialization year so that
            # all five FY fields have exactly the same time coordinate.
            cdo -O \
                settaxis,${start}-07-01,00:00:00,1year \
                -selyear,${target_year} \
                "$infile" \
                "$tmp"

            selected+=("$tmp")
        done


        # Direct FY1-FY5 mean for this initialization
        rawmean="${tmpdir}/mean_raw_${start}.nc"
        outfile_tmp="${tmpdir}/FWI95d_INIT_${member}_${start}.nc"

        cdo -O ensmean "${selected[@]}" "$rawmean"

        # Give variable an explicit name
        cdo -O setname,FWI95d \
            "$rawmean" \
            "$outfile_tmp"

    done


    # Merge the 60 initialization dates
    outfile="${outdir}/FWI95d_INIT_FY1-FY5_${member}_init1960-2019.nc"

    cdo -O mergetime \
        "${tmpdir}"/FWI95d_INIT_${member}_*.nc \
        "$outfile"

    rm -rf "$tmpdir"

    echo "Saved: $outfile"
done

echo "DONE."