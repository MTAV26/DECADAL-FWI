#!/usr/bin/env bash
set -euo pipefail

# ==============================================================================
# ERA5 GLOBAL daily-mean 10 m wind speed — block-download version
#
# Correct calculation:
#   1) hourly scalar wind speed = sqrt(u10^2 + v10^2)
#   2) daily mean of the hourly scalar wind speed
#
# Data are downloaded globally from ERA5 in multi-day blocks to reduce the
# number of CDS requests. No spatial subset and no remapping are applied.
#
# Defaults:
#   period      = 1961-01-01 to 2024-12-31
#   block size  = 10 days per CDS request
#
# Products:
#   monthly/sfcWind_daily_global_YYYYMM.nc
#   yearly/sfcWind_daily_global_YYYY.nc
#
# Restart behaviour:
#   - valid monthly files are skipped;
#   - valid processed block files are reused;
#   - raw downloads use a .part suffix and are renamed only after completion.
# ==============================================================================

ROOT="/diskonfire/ERA5/daily_values/sfcWind_v2"

START_YEAR="${START_YEAR:-1961}"
END_YEAR="${END_YEAR:-2024}"
START_MONTH="${START_MONTH:-1}"
END_MONTH="${END_MONTH:-12}"
BLOCK_DAYS="${BLOCK_DAYS:-10}"

KEEP_RAW="${KEEP_RAW:-false}"
KEEP_BLOCKS="${KEEP_BLOCKS:-false}"
BUILD_YEARLY="${BUILD_YEARLY:-true}"

RAW_DIR="${ROOT}/raw_blocks"
BLOCK_DIR="${ROOT}/block_cache"
MONTHLY_DIR="${ROOT}/monthly"
YEARLY_DIR="${ROOT}/yearly"
TMP_ROOT="${ROOT}/tmp_blocks"
LOG_DIR="${ROOT}/logs"

mkdir -p \
    "$RAW_DIR" \
    "$BLOCK_DIR" \
    "$MONTHLY_DIR" \
    "$YEARLY_DIR" \
    "$TMP_ROOT" \
    "$LOG_DIR"

if ! [[ "$BLOCK_DAYS" =~ ^[0-9]+$ ]] || (( BLOCK_DAYS < 1 || BLOCK_DAYS > 15 )); then
    echo "ERROR: BLOCK_DAYS must be an integer between 1 and 15." >&2
    exit 1
fi

# Use the CDS credentials already used in the previous workflow.
if [[ -f "${HOME}/.cdsapirc_cds" ]]; then
    cp "${HOME}/.cdsapirc_cds" "${HOME}/.cdsapirc"
elif [[ ! -f "${HOME}/.cdsapirc" ]]; then
    echo "ERROR: neither ~/.cdsapirc_cds nor ~/.cdsapirc exists." >&2
    exit 1
fi

for cmd in python3 cdo; do
    command -v "$cmd" >/dev/null 2>&1 || {
        echo "ERROR: required command not found: $cmd" >&2
        exit 1
    }
done

CURRENT_TMP=""

cleanup() {
    if [[ -n "${CURRENT_TMP}" && -d "${CURRENT_TMP}" ]]; then
        rm -rf "${CURRENT_TMP}"
    fi
}
trap cleanup EXIT INT TERM

days_in_month() {
    local year="$1"
    local month="$2"
    date -d "${year}-$(printf '%02d' "$month")-01 +1 month -1 day" +%d
}

download_block() {
    local year="$1"
    local month="$2"
    local first_day="$3"
    local last_day="$4"
    local target="$5"

    local part="${target}.part"
    rm -f "$part"

    python3 - "$year" "$month" "$first_day" "$last_day" "$part" <<'PY'
import sys
from pathlib import Path

import cdsapi

year = int(sys.argv[1])
month = int(sys.argv[2])
first_day = int(sys.argv[3])
last_day = int(sys.argv[4])
target = Path(sys.argv[5])

request = {
    "product_type": ["reanalysis"],
    "variable": [
        "10m_u_component_of_wind",
        "10m_v_component_of_wind",
    ],
    "year": [f"{year:04d}"],
    "month": [f"{month:02d}"],
    "day": [f"{day:02d}" for day in range(first_day, last_day + 1)],
    "time": [f"{hour:02d}:00" for hour in range(24)],
    "data_format": "grib",
    "download_format": "unarchived",
}

target.parent.mkdir(parents=True, exist_ok=True)
client = cdsapi.Client()
client.retrieve("reanalysis-era5-single-levels", request).download(str(target))
PY

    mv "$part" "$target"
}

process_block() {
    local year="$1"
    local month_num="$2"
    local first_day="$3"
    local last_day="$4"

    local month first last tag expected_days
    month=$(printf "%02d" "$month_num")
    first=$(printf "%02d" "$first_day")
    last=$(printf "%02d" "$last_day")
    tag="${year}${month}_${first}-${last}"
    expected_days=$(( last_day - first_day + 1 ))

    local raw_grib="${RAW_DIR}/era5_u10_v10_${tag}.grib"
    local block_out="${BLOCK_DIR}/sfcWind_daily_global_${tag}.nc"

    # Reuse a valid completed block.
    if [[ -f "$block_out" ]]; then
        local ntime
        ntime=$(cdo -s ntime "$block_out" 2>/dev/null || echo 0)
        if [[ "$ntime" -eq "$expected_days" ]]; then
            echo "  SKIP block ${first}-${last}: valid output exists (${ntime} days)."
            return
        fi
        echo "  WARNING: invalid block output; rebuilding $block_out"
        rm -f "$block_out"
    fi

    echo "  Block ${year}-${month}-${first} to ${year}-${month}-${last}"

    if [[ ! -f "$raw_grib" ]]; then
        echo "    Downloading ${expected_days} days in one CDS request..."
        download_block "$year" "$month_num" "$first_day" "$last_day" "$raw_grib"
    else
        echo "    Using existing raw block."
    fi

    CURRENT_TMP=$(mktemp -d "${TMP_ROOT}/wind_${tag}_XXXXXX")

    local raw_nc="${CURRENT_TMP}/era5_uv.nc"
    local u_nc="${CURRENT_TMP}/u10.nc"
    local v_nc="${CURRENT_TMP}/v10.nc"
    local hourly_speed="${CURRENT_TMP}/sfcWind_hourly.nc"
    local daily_tmp="${CURRENT_TMP}/sfcWind_daily.nc"
    local final_tmp="${CURRENT_TMP}/sfcWind_final.nc"

    echo "    Converting GRIB to NetCDF..."
    cdo -O -f nc copy "$raw_grib" "$raw_nc"

    # Variable names may be 10u/10v or u10/v10 depending on ecCodes/CDO.
    local names u_name v_name
    names=$(cdo -s showname "$raw_nc" | tr ' ' '\n' | sed '/^$/d')
    u_name=$(printf "%s\n" "$names" | grep -E '^(10u|u10)$' | head -n 1 || true)
    v_name=$(printf "%s\n" "$names" | grep -E '^(10v|v10)$' | head -n 1 || true)

    if [[ -z "$u_name" || -z "$v_name" ]]; then
        echo "ERROR: unable to identify u10/v10 in $raw_nc" >&2
        echo "Variables found:" >&2
        cdo -s showname "$raw_nc" >&2
        exit 1
    fi

    cdo -O selname,"$u_name" "$raw_nc" "$u_nc"
    cdo -O selname,"$v_name" "$raw_nc" "$v_nc"

    echo "    Computing hourly scalar speed..."
    cdo -O -b F32 sqrt \
        -add \
        -sqr "$u_nc" \
        -sqr "$v_nc" \
        "$hourly_speed"

    echo "    Computing daily means..."
    cdo -O -b F32 daymean "$hourly_speed" "$daily_tmp"

    cdo -O -b F32 setname,sfcWind "$daily_tmp" "$final_tmp"

    # Optional metadata cleanup.
    if command -v ncatted >/dev/null 2>&1; then
        ncatted -O \
            -a standard_name,sfcWind,o,c,"wind_speed" \
            -a long_name,sfcWind,o,c,"Daily mean 10 m wind speed" \
            -a units,sfcWind,o,c,"m s-1" \
            "$final_tmp"
    fi

    local actual_days
    actual_days=$(cdo -s ntime "$final_tmp")
    if [[ "$actual_days" -ne "$expected_days" ]]; then
        echo "ERROR: block ${tag} has ${actual_days} daily records; expected ${expected_days}." >&2
        exit 1
    fi

    mv "$final_tmp" "$block_out"
    echo "    Saved processed block: $block_out"

    rm -rf "$CURRENT_TMP"
    CURRENT_TMP=""

    if [[ "$KEEP_RAW" != "true" ]]; then
        rm -f "$raw_grib"
    fi
}

build_month() {
    local year="$1"
    local month_num="$2"

    local month total_days monthly_out
    month=$(printf "%02d" "$month_num")
    total_days=$(days_in_month "$year" "$month_num")
    total_days=$((10#$total_days))
    monthly_out="${MONTHLY_DIR}/sfcWind_daily_global_${year}${month}.nc"

    if [[ -f "$monthly_out" ]]; then
        local ntime
        ntime=$(cdo -s ntime "$monthly_out" 2>/dev/null || echo 0)
        if [[ "$ntime" -eq "$total_days" ]]; then
            echo "SKIP ${year}-${month}: valid monthly file exists (${ntime} days)."
            return
        fi
        echo "WARNING: invalid monthly file; rebuilding $monthly_out"
        rm -f "$monthly_out"
    fi

    echo
    echo "=================================================================="
    echo "MONTH ${year}-${month}"
    echo "=================================================================="

    local block_files=()
    local first last first_fmt last_fmt

    first=1
    while (( first <= total_days )); do
        last=$(( first + BLOCK_DAYS - 1 ))
        if (( last > total_days )); then
            last="$total_days"
        fi

        process_block "$year" "$month_num" "$first" "$last"

        first_fmt=$(printf "%02d" "$first")
        last_fmt=$(printf "%02d" "$last")
        block_files+=(
            "${BLOCK_DIR}/sfcWind_daily_global_${year}${month}_${first_fmt}-${last_fmt}.nc"
        )

        first=$(( last + 1 ))
    done

    for f in "${block_files[@]}"; do
        if [[ ! -f "$f" ]]; then
            echo "ERROR: missing expected processed block: $f" >&2
            exit 1
        fi
    done

    local monthly_tmp="${monthly_out}.tmp.nc"
    cdo -O mergetime "${block_files[@]}" "$monthly_tmp"

    local actual_days
    actual_days=$(cdo -s ntime "$monthly_tmp")
    if [[ "$actual_days" -ne "$total_days" ]]; then
        echo "ERROR: ${year}-${month} has ${actual_days} daily records; expected ${total_days}." >&2
        rm -f "$monthly_tmp"
        exit 1
    fi

    mv "$monthly_tmp" "$monthly_out"
    echo "Saved monthly file: $monthly_out"

    if [[ "$KEEP_BLOCKS" != "true" ]]; then
        rm -f "${block_files[@]}"
    fi
}

build_year() {
    local year="$1"
    local yearly_out="${YEARLY_DIR}/sfcWind_daily_global_${year}.nc"

    local monthly_files=()
    local month
    for month in $(seq 1 12); do
        monthly_files+=(
            "${MONTHLY_DIR}/sfcWind_daily_global_${year}$(printf '%02d' "$month").nc"
        )
    done

    for f in "${monthly_files[@]}"; do
        if [[ ! -f "$f" ]]; then
            echo "Year ${year}: monthly files incomplete; yearly merge postponed."
            return
        fi
    done

    local expected_days
    if date -d "${year}-02-29" >/dev/null 2>&1; then
        expected_days=366
    else
        expected_days=365
    fi

    if [[ -f "$yearly_out" ]]; then
        local ntime
        ntime=$(cdo -s ntime "$yearly_out" 2>/dev/null || echo 0)
        if [[ "$ntime" -eq "$expected_days" ]]; then
            echo "SKIP year ${year}: valid yearly file exists (${ntime} days)."
            return
        fi
        rm -f "$yearly_out"
    fi

    local yearly_tmp="${yearly_out}.tmp.nc"
    cdo -O mergetime "${monthly_files[@]}" "$yearly_tmp"

    local actual_days
    actual_days=$(cdo -s ntime "$yearly_tmp")
    if [[ "$actual_days" -ne "$expected_days" ]]; then
        echo "ERROR: year ${year} has ${actual_days} days; expected ${expected_days}." >&2
        rm -f "$yearly_tmp"
        exit 1
    fi

    mv "$yearly_tmp" "$yearly_out"
    echo "Saved yearly file: $yearly_out"
}

for year in $(seq "$START_YEAR" "$END_YEAR"); do
    for month in $(seq "$START_MONTH" "$END_MONTH"); do
        build_month "$year" "$month"
    done

    if [[ "$BUILD_YEARLY" == "true" ]]; then
        build_year "$year"
    fi
done

echo
echo "DONE."
echo "Global monthly products: $MONTHLY_DIR"
echo "Global yearly products:  $YEARLY_DIR"
echo "Block length used:        $BLOCK_DAYS days"
echo
echo "No spatial subset and no remapping have been applied."
