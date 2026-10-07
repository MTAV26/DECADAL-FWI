#!/bin/bash
set -euo pipefail

IN="/diskonfire/ERA5/daily_values/hurs_corrected/rhmean_corrected_1961_2024_1degree.nc"

# Existing ERA5 RH file is used only as target grid definition
TARGET="/diskonfire/ERA5/daily_values/hurs/rhmean_remapped.nc"

OUT="/diskonfire/ERA5/daily_values/hurs_corrected/rhmean_remapped_corrected.nc"

echo "Remapping corrected RH to decadal European grid..."

time cdo -O remapbil,"$TARGET" "$IN" "$OUT"

echo
echo "GRID:"
cdo -s griddes "$OUT" | \
grep -E 'gridtype|xsize|ysize|xfirst|xinc|yfirst|yinc'

echo
echo "NTIME:"
cdo -s ntime "$OUT"

echo
echo "DATES:"
cdo -s showdate "$OUT" | awk '{print $1, $NF}'

echo
echo "MEAN RH:"
cdo -s output -fldmean -timmean "$OUT"

echo
echo "OUTPUT:"
ls -lh "$OUT"
