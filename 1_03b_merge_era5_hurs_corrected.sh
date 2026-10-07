#!/bin/bash
set -euo pipefail

DIR="/diskonfire/ERA5/daily_values/hurs_corrected"

OUT="${DIR}/rhmean_corrected_1961_2024_1degree.nc"
TMP="${DIR}/rhmean_corrected_1961_2024_1degree.tmp.nc"

EXPECTED_MONTHS=768
EXPECTED_DAYS=23376

mapfile -t FILES < <(
    find "$DIR" -maxdepth 1 -type f \
      -name 'rhmean_??????_1degree.nc' \
      | sort
)

NFILES=${#FILES[@]}

echo "Monthly files found: $NFILES"

if [ "$NFILES" -ne "$EXPECTED_MONTHS" ]; then
    echo "ERROR: expected $EXPECTED_MONTHS monthly files, found $NFILES"
    exit 1
fi

echo "First file: ${FILES[0]}"
echo "Last file : ${FILES[$((NFILES-1))]}"

rm -f "$TMP"

echo
echo "Merging RH files..."
time cdo -O mergetime "${FILES[@]}" "$TMP"

NTIME=$(cdo -s ntime "$TMP")
FIRST=$(cdo -s showdate "$TMP" | awk '{print $1}')
LAST=$(cdo -s showdate "$TMP" | awk '{print $NF}')

echo
echo "Number of days: $NTIME"
echo "First date: $FIRST"
echo "Last date : $LAST"

if [ "$NTIME" -ne "$EXPECTED_DAYS" ]; then
    echo "ERROR: expected $EXPECTED_DAYS days, found $NTIME"
    exit 1
fi

if [ "$FIRST" != "1961-01-01" ] || [ "$LAST" != "2024-12-31" ]; then
    echo "ERROR: unexpected temporal coverage"
    exit 1
fi

mv "$TMP" "$OUT"

echo
echo "Mean RH:"
cdo -s output -fldmean -timmean "$OUT"

echo
echo "Output:"
ls -lh "$OUT"

echo
echo "MERGE COMPLETED SUCCESSFULLY"
