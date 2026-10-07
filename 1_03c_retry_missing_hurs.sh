#!/bin/bash
set -euo pipefail

REBUILD="/diskonfire/ERA5/daily_values/hurs/rebuild_hurs_corrected.sh"
OUT="/diskonfire/ERA5/daily_values/hurs_corrected"

MONTHS=(
  1999-06
  2016-06
  2016-07
  2016-08
  2016-09
  2016-10
  2016-11
  2016-12
  2017-01
  2017-02
  2017-03
  2017-04
  2017-05
  2017-06
  2017-07
  2017-08
  2017-09
)

for YM in "${MONTHS[@]}"; do

    YEAR=${YM%-*}
    MONTH=${YM#*-}
    MNUM=$((10#$MONTH))

    FILE="${OUT}/rhmean_${YEAR}${MONTH}_1degree.nc"

    if [ -s "$FILE" ]; then
        echo "$YM already exists -> skipping"
        continue
    fi

    echo
    echo "=========================================="
    echo "RETRY $YM"
    echo "=========================================="

    "$REBUILD" "$YEAR" "$YEAR" "$MNUM" "$MNUM"

done

echo
echo "Files now present:"
find "$OUT" -maxdepth 1 -type f \
  -name 'rhmean_??????_1degree.nc' | wc -l
