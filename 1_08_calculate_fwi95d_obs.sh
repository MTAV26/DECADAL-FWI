#!/bin/bash

set -euo pipefail

# ============================================================
# ERA5 FWI95d
# Daily FWI: 1961-2024
# Reference period for p95: 1991-2020
# ============================================================

OBS_FILE="/diskonfire/ERA5/daily_values/FWI_1961_2024_1degree_europe_corrected.nc"
WORKDIR="/diskonfire/ERA5/daily_values/processed_europe"

mkdir -p "$WORKDIR"

REF_START="1991-01-01"
REF_END="2020-12-31"

REF_FILE="${WORKDIR}/FWI_ref_1991_2020_corrected.nc"
PCTL95_FILE="${WORKDIR}/FWI_obs_pctl95_1991_2020_corrected.nc"
FWI95D_FILE="${WORKDIR}/FWI_obs_FWI95d_1961_2024_ref1991_2020_corrected.nc"


echo "1) Extracting reference period 1991-2020..."
cdo -O seldate,${REF_START},${REF_END} \
    "$OBS_FILE" \
    "$REF_FILE"


echo "2) Calculating grid-point 95th percentile..."
cdo -O timpctl,95 \
    "$REF_FILE" \
    -timmin "$REF_FILE" \
    -timmax "$REF_FILE" \
    "$PCTL95_FILE"


echo "3) Counting annual days with FWI > p95..."
cdo -O -yearsum -gt \
    "$OBS_FILE" \
    "$PCTL95_FILE" \
    "$FWI95D_FILE"


echo "4) Renaming variable FWI -> FWI95d..."
TMP_FILE="${FWI95D_FILE%.nc}_tmp.nc"

cdo -O chname,FWI,FWI95d \
    "$FWI95D_FILE" \
    "$TMP_FILE"

mv "$TMP_FILE" "$FWI95D_FILE"


echo "DONE."
echo "FWI95d saved in:"
echo "$FWI95D_FILE"