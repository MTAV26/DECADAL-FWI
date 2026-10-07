#!/bin/bash

# Directorios de entrada y salida
indir="/diskonfire/Decadal/FWI"
outdir="/diskonfire/Decadal/processed_INIT"
mkdir -p "$outdir"

# Miembros (ajusta el rango según tus datos)
members=$(seq -f "r%gi1p1f1" 17 40)

# Forecast years que quieras procesar
forecast_years=(fy1 fy2 fy3 fy4 fy5)

# Años para calcular el percentil 95 (siempre 1991–2020)
years_pctl=$(seq 1991 2020)

# Parámetros del periodo móvil de 60 años
years_all_start=1961
window=60

for member in $members; do
  echo "=== Procesando miembro $member ==="

  for fy in "${forecast_years[@]}"; do
    echo "-- Forecast year: $fy --"

    # Crear dir temporal
    fy_outdir="${outdir}/${member}/${fy}"
    mkdir -p "$fy_outdir"

    # Extrae el número (1…5)
    fynum=${fy/fy/}

    #
    # 1) Extraer cada año 1991–2020 para el percentil
    #
    for year in $years_pctl; do
      start_year=$(( year - fynum ))
      end_year=$(( start_year + 10 ))
      infile="${indir}/FWI_CMCC-CM2-SR5_dcppA-hindcast_s${start_year}-${member}_gn_${start_year}1101-${end_year}1231_europe.nc"
      tmp_pctl="${fy_outdir}/FWI_${fy}_${year}.nc"

      if [ -f "$infile" ]; then
        echo "  [PCTL] Extrayendo $year → $(basename "$tmp_pctl")"
        cdo seldate,${year}-01-01,${year}-12-31 "$infile" "$tmp_pctl"
      else
        echo "  ⚠️ No existe $infile"
      fi
    done

    # 2) Concatenar 1991–2020 y calcular pctl95
    merged_pctl="${fy_outdir}/FWI_${member}_${fy}_1991-2020.nc"
    echo "  Concatenando 1991–2020 → $(basename "$merged_pctl")"
    cdo mergetime "${fy_outdir}/FWI_${fy}_"*.nc "$merged_pctl"

    file_pctl95="${outdir}/FWI_${member}_${fy}_pctl95.nc"
    echo "  Calculando 95º percentil → $(basename "$file_pctl95")"
    cdo timpctl,95 "$merged_pctl" -timmin "$merged_pctl" -timmax "$merged_pctl" "$file_pctl95"
    #
    # 3) Periodo móvil de 60 años para FWI95d
    #
    start_all=$(( years_all_start + fynum - 1 ))
    end_all=$(( start_all + window - 1 ))
    echo "  Periodo de estudio para FWI95d: ${start_all}–${end_all}"

    # Borra viejos temporales
    rm -f "${fy_outdir}/FWI_${fy}_"*"_FWI95d.nc"

    for year in $(seq "$start_all" "$end_all"); do
      start_year=$(( year - fynum ))
      end_year=$(( start_year + 10 ))
      infile="${indir}/FWI_CMCC-CM2-SR5_dcppA-hindcast_s${start_year}-${member}_gn_${start_year}1101-${end_year}1231_europe.nc"
      tmp_year="${fy_outdir}/FWI_${fy}_${year}_tmp.nc"
      out_year="${fy_outdir}/FWI_${fy}_${year}_FWI95d.nc"

      if [ -f "$infile" ]; then
        echo "    Año $year: extrayendo → $(basename "$tmp_year")"
        cdo seldate,${year}-01-01,${year}-12-31 "$infile" "$tmp_year"

        echo "    Año $year: contando días > pctl95 → $(basename "$out_year")"
        cdo -yearsum -gt "$tmp_year" "$file_pctl95" "$out_year"

        rm -f "$tmp_year"
      else
        echo "    ⚠️ No existe $infile"
      fi
    done

    # 4) Concatenar los 60 archivos FWI95d
    FWI95d_out="${outdir}/FWI_${member}_${fy}_FWI95d_${start_all}-${end_all}.nc"
    echo "  Concatenando FWI95d → $(basename "$FWI95d_out")"
    cdo mergetime "${fy_outdir}/FWI_${fy}_"*"_FWI95d.nc" "$FWI95d_out"

    echo "-- Fin de fy $fy --"
    echo
  done
done

