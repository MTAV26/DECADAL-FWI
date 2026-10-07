#!/bin/bash
set -euo pipefail

IN="/diskonfire/ERA5/daily_values/processed_europe/FWI_obs_FWI95d_1961_2024_ref1991_2020_corrected.nc"

OUTDIR="/diskonfire/Decadal/verification_new"
TMP="${OUTDIR}/FWI95d_OBS_runmean5_1963-2022_corrected_tmp.nc"
OUT="${OUTDIR}/FWI95d_OBS_runmean5_1963-2022_corrected.nc"

echo "Calculating centered 5-year running mean..."

cdo -O runmean,5 "$IN" "$TMP"

# Force exactly the same annual time axis as the original verification file
cdo -O settaxis,1963-07-01,00:00:00,1year \
    "$TMP" \
    "$OUT"

rm -f "$TMP"

echo
echo "NTIME:"
cdo -s ntime "$OUT"

echo "DATES:"
cdo -s showdate "$OUT" | awk '{print $1, $NF}'

echo "VARIABLE:"
cdo -s showname "$OUT"

echo "OUTPUT:"
ls -lh "$OUT"
