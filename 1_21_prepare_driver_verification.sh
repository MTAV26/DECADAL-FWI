#!/usr/bin/env bash
set -euo pipefail
shopt -s nullglob

# ============================================================
# Prepare FY1-FY5 verification series for the meteorological
# drivers of FWI:
#
#   tas      : °C
#   hurs     : %
#   pr       : mm/day
#   sfcWind  : m/s
#
# Verification:
#   forecast window 1961-1965 -> verification year 1963
#   ...
#   forecast window 2020-2024 -> verification year 2022
#
# OBS/HIST:
#   annual mean -> 5-year running mean
#
# INIT:
#   annual means FY1...FY5 -> mean of the five annual values
# ============================================================


ROOT="/diskonfire/Decadal"

OUT="${ROOT}/verification_drivers_new"
TMP="${OUT}/.tmp"

MODEL="CMCC-CM2-SR5"

VARS=(tas hurs pr sfcWind)


# ============================================================
# MEMBERS
# ============================================================

INIT_MEMBERS=()
for n in $(seq 17 40); do
    INIT_MEMBERS+=("r${n}i1p1f1")
done

HIST_MEMBERS=()
for n in $(seq 2 11); do
    HIST_MEMBERS+=("r${n}i1p2f1")
done


# ============================================================
# ERA5 SOURCES
# ============================================================

declare -A OBS_FILE

OBS_FILE[tas]="/diskonfire/ERA5/daily_values/tasmean/tasmean_remapped.nc"
OBS_FILE[hurs]="/diskonfire/ERA5/daily_values/hurs_corrected/rhmean_remapped_corrected.nc"
OBS_FILE[pr]="/diskonfire/ERA5/daily_values/pr/precip_remapped.nc"
OBS_FILE[sfcWind]="/diskonfire/ERA5/daily_values/sfcWind_v2/sfcWind_remapped.nc"


# ============================================================
# HELPERS
# ============================================================

die() {
    echo
    echo "ERROR: $*" >&2
    exit 1
}


valid60() {

    local f="$1"

    [[ -f "$f" ]] || return 1

    local n
    n=$(cdo -s ntime "$f" 2>/dev/null || echo 0)

    [[ "$n" -eq 60 ]]
}


# ------------------------------------------------------------
# Process OBS/HIST complete 1961-2024 time series
#
# daily -> annual mean -> 5-year running mean
#
# This produces 60 values:
# 1961-65 ... 2020-24
# ------------------------------------------------------------

process_long_series() {

    local src="$1"
    local var="$2"
    local kind="$3"
    local outfile="$4"

    [[ -f "$src" ]] || die "Missing source: $src"

    case "${var}:${kind}" in

        tas:*)
            # K -> degC
            cdo -O -b F32 \
                settaxis,1963-07-01,00:00:00,1year \
                -setname,tas \
                -subc,273.15 \
                -runmean,5 \
                -yearmean \
                -selyear,1961/2024 \
                "$src" \
                "$outfile"
            ;;


        pr:OBS)
            # ERA5 daily accumulated precipitation:
            # m/day -> mm/day
            cdo -O -b F32 \
                settaxis,1963-07-01,00:00:00,1year \
                -setname,pr \
                -mulc,1000 \
                -runmean,5 \
                -yearmean \
                -selyear,1961/2024 \
                "$src" \
                "$outfile"
            ;;


        pr:HIST)
            # CMIP6 precipitation flux:
            # kg m-2 s-1 -> mm/day
            cdo -O -b F32 \
                settaxis,1963-07-01,00:00:00,1year \
                -setname,pr \
                -mulc,86400 \
                -runmean,5 \
                -yearmean \
                -selyear,1961/2024 \
                "$src" \
                "$outfile"
            ;;


        hurs:*)
            # Relative humidity already in %
            cdo -O -b F32 \
                settaxis,1963-07-01,00:00:00,1year \
                -setunit,% \
                -setname,hurs \
                -runmean,5 \
                -yearmean \
                -selyear,1961/2024 \
                "$src" \
                "$outfile"
            ;;

        sfcWind:*)
            # Wind speed already in m/s
            cdo -O -b F32 \
                settaxis,1963-07-01,00:00:00,1year \
                -setunit,"m s-1" \
                -setname,sfcWind \
                -runmean,5 \
                -yearmean \
                -selyear,1961/2024 \
                "$src" \
                "$outfile"
            ;;


        *)
            die "Unknown variable/kind: ${var}:${kind}"
            ;;
    esac
}


# ------------------------------------------------------------
# Process one initialized forecast:
#
# select FY1-FY5
# -> annual means
# -> equal-weight mean of the five annual means
# ------------------------------------------------------------

process_init_window() {

    local src="$1"
    local var="$2"
    local y1="$3"
    local y5="$4"
    local center="$5"
    local outfile="$6"

    [[ -f "$src" ]] || die "Missing INIT source: $src"

    case "$var" in

        tas)
            cdo -O -b F32 \
                settaxis,${center}-07-01,00:00:00,1year \
                -setname,tas \
                -subc,273.15 \
                -timmean \
                -yearmean \
                -selyear,${y1}/${y5} \
                "$src" \
                "$outfile"
            ;;


        hurs)
            cdo -O -b F32 \
                settaxis,${center}-07-01,00:00:00,1year \
                -setname,hurs \
                -timmean \
                -yearmean \
                -selyear,${y1}/${y5} \
                "$src" \
                "$outfile"
            ;;


        pr)
            cdo -O -b F32 \
                settaxis,${center}-07-01,00:00:00,1year \
                -setname,pr \
                -mulc,86400 \
                -timmean \
                -yearmean \
                -selyear,${y1}/${y5} \
                "$src" \
                "$outfile"
            ;;


        sfcWind)
            cdo -O -b F32 \
                settaxis,${center}-07-01,00:00:00,1year \
                -setname,sfcWind \
                -timmean \
                -yearmean \
                -selyear,${y1}/${y5} \
                "$src" \
                "$outfile"
            ;;


        *)
            die "Unknown INIT variable: $var"
            ;;
    esac
}


# ============================================================
# DIRECTORIES
# ============================================================

mkdir -p "$OUT" "$TMP"

for var in "${VARS[@]}"; do

    mkdir -p \
        "$OUT/OBS/$var" \
        "$OUT/INIT/$var" \
        "$OUT/HIST/$var"

done


# ============================================================
# 1. OBS / ERA5
# ============================================================

echo
echo "============================================================"
echo "1. ERA5 / OBS"
echo "============================================================"

for var in "${VARS[@]}"; do

    echo
    echo "=== OBS: $var ==="

    src="${OBS_FILE[$var]}"

    [[ -f "$src" ]] || die "Missing ERA5 source: $src"

    final="$OUT/OBS/$var/${var}_OBS_FY1-FY5_1963-2022.nc"

    if valid60 "$final"; then

        echo "SKIP: valid 60-step file already exists"

    else

        rm -f "$final"

        process_long_series \
            "$src" \
            "$var" \
            "OBS" \
            "$final"

        echo "Saved: $final"

    fi

done


# ============================================================
# 2. INITIALIZED DCPP-A
# ============================================================

echo
echo "============================================================"
echo "2. INITIALIZED / INIT"
echo "============================================================"

for var in "${VARS[@]}"; do

    echo
    echo "############################################################"
    echo "INIT VARIABLE: $var"
    echo "############################################################"

    for member in "${INIT_MEMBERS[@]}"; do

        echo
        echo "=== $var : $member ==="

        final="$OUT/INIT/$var/${var}_INIT_FY1-FY5_${member}_1963-2022.nc"

        if valid60 "$final"; then

            echo "SKIP: valid 60-step file already exists"
            continue

        fi

        rm -f "$final"

        tdir="$TMP/INIT_${var}_${member}"

        rm -rf "$tdir"
        mkdir -p "$tdir"


        for start in $(seq 1960 2019); do

            y1=$((start + 1))
            y5=$((start + 5))
            center=$((start + 3))


            # Exact pattern verified on server:
            #
            # tas_day_CMCC-CM2-SR5_dcppA-hindcast_
            # s1960-r17i1p1f1_gn_19601101-19701231_europe.nc

            matches=(
                "${ROOT}/${var}_day/${var}_day_${MODEL}_dcppA-hindcast_s${start}-${member}_gn_"*_europe.nc
            )

            if [[ ${#matches[@]} -ne 1 ]]; then

                echo
                echo "Expected exactly one INIT file."
                echo "var     = $var"
                echo "start   = $start"
                echo "member  = $member"
                echo "matches = ${#matches[@]}"

                if [[ ${#matches[@]} -gt 0 ]]; then
                    printf '  %s\n' "${matches[@]}"
                fi

                exit 1
            fi

            src="${matches[0]}"

            echo "  init ${start}: ${y1}-${y5} -> ${center}"

            tmpfile="$tdir/${var}_${member}_${center}.nc"

            process_init_window \
                "$src" \
                "$var" \
                "$y1" \
                "$y5" \
                "$center" \
                "$tmpfile"

        done


        # Merge 60 initialization windows
        cdo -O mergetime \
            "$tdir"/${var}_${member}_*.nc \
            "$final"


        n=$(cdo -s ntime "$final")

        if [[ "$n" -ne 60 ]]; then
            die "INIT $var $member has $n timesteps instead of 60"
        fi


        rm -rf "$tdir"

        echo "Saved: $final"

    done

done


# ============================================================
# 3. HISTORICAL + SSP2-4.5 / NO-INIT
# ============================================================

echo
echo "============================================================"
echo "3. HIST+245 / NO-INIT"
echo "============================================================"

for var in "${VARS[@]}"; do

    echo
    echo "############################################################"
    echo "HIST VARIABLE: $var"
    echo "############################################################"

    for member in "${HIST_MEMBERS[@]}"; do

        echo
        echo "=== $var : $member ==="


        # Exact pattern verified on server:
        #
        # tas_merged_CMCC-CM2-SR5_r2i1p2f1_
        # 19600101-20241231_regridded.nc

        src="${ROOT}/processed_NO-INIT/${var}_merged_${MODEL}_${member}_19600101-20241231_regridded.nc"

        [[ -f "$src" ]] || die "Missing HIST source: $src"


        final="$OUT/HIST/$var/${var}_HIST_FY1-FY5_${member}_1963-2022.nc"


        if valid60 "$final"; then

            echo "SKIP: valid 60-step file already exists"
            continue

        fi


        rm -f "$final"

        process_long_series \
            "$src" \
            "$var" \
            "HIST" \
            "$final"


        n=$(cdo -s ntime "$final")

        if [[ "$n" -ne 60 ]]; then
            die "HIST $var $member has $n timesteps instead of 60"
        fi


        echo "Saved: $final"

    done

done


# ============================================================
# 4. FINAL VALIDATION
# ============================================================

echo
echo "============================================================"
echo "4. FINAL VALIDATION"
echo "============================================================"

count=0
bad=0

while IFS= read -r f; do

    count=$((count + 1))

    n=$(cdo -s ntime "$f" 2>/dev/null || echo 0)

    if [[ "$n" -ne 60 ]]; then

        echo "BAD ntime=$n : $f"
        bad=$((bad + 1))

    fi

done < <(
    find "$OUT" \
        -type f \
        -name "*.nc" \
        ! -path "$TMP/*" \
        | sort
)


nobs=$(find "$OUT/OBS"  -type f -name "*.nc" | wc -l)
ninit=$(find "$OUT/INIT" -type f -name "*.nc" | wc -l)
nhist=$(find "$OUT/HIST" -type f -name "*.nc" | wc -l)


echo
echo "Total NetCDF files : $count"
echo "Expected           : 140"

echo
echo "OBS  : $nobs  (expected 4)"
echo "INIT : $ninit (expected 96)"
echo "HIST : $nhist (expected 40)"

echo
echo "Files with ntime != 60: $bad"


if [[ "$count" -ne 140 ]]; then
    die "Expected 140 NetCDF files, found $count"
fi

if [[ "$nobs" -ne 4 ]]; then
    die "Expected 4 OBS files, found $nobs"
fi

if [[ "$ninit" -ne 96 ]]; then
    die "Expected 96 INIT files, found $ninit"
fi

if [[ "$nhist" -ne 40 ]]; then
    die "Expected 40 HIST files, found $nhist"
fi

if [[ "$bad" -ne 0 ]]; then
    die "$bad files have incorrect time dimension"
fi


# ============================================================
# 5. QUICK METADATA CHECK
# ============================================================

echo
echo "============================================================"
echo "5. SAMPLE METADATA"
echo "============================================================"

for var in "${VARS[@]}"; do

    f="$OUT/OBS/$var/${var}_OBS_FY1-FY5_1963-2022.nc"

    echo
    echo "--- $var ---"

    echo -n "name  : "
    cdo -s showname "$f"

    echo -n "unit  : "
    cdo -s showunit "$f"

    echo -n "ntime : "
    cdo -s ntime "$f"

    echo -n "years : "
    cdo -s showyear "$f"

done


rm -rf "$TMP"


echo
echo "============================================================"
echo "DONE"
echo "============================================================"
echo
echo "Output:"
echo "$OUT"
echo

