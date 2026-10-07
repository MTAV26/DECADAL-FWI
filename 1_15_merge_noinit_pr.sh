#!/bin/bash

# === CONFIGURACIÓN ===
model="CMCC-CM2-SR5"
var="pr"
grid="gn"
dir_hist="/diskonfire/Decadal/hist_pr_day"
dir_ssp_base="/diskonfire/Decadal/LEONE/daily_ssp245_leone/CMIP6/ScenarioMIP/CMCC/${model}/ssp245"
mask_file="/diskonfire/Decadal/landsea_mask.nc"
output_dir="/diskonfire/Decadal/processed_NO-INIT"
mkdir -p ${output_dir}

# === LISTA DE REALIZACIONES ===
members=(r2i1p2f1 r3i1p2f1 r4i1p2f1 r5i1p2f1 r6i1p2f1 r7i1p2f1 r8i1p2f1 r9i1p2f1 r10i1p2f1 r11i1p2f1)

# === BUCLE PRINCIPAL ===
for realization in "${members[@]}"; do
  echo "\n>>> Procesando miembro: ${realization}"

  hist_merged="${output_dir}/${var}_historical_${model}_${realization}_19600101-20141231.nc"
  ssp_merged="${output_dir}/${var}_ssp245_${model}_${realization}_20150101-20241231.nc"
  final_merged="${output_dir}/${var}_merged_${model}_${realization}_19600101-20241231.nc"
  final_regridded="${output_dir}/${var}_merged_${model}_${realization}_19600101-20241231_regridded.nc"

  if [ -f "${final_regridded}" ]; then
    echo "✅ Ya existe: ${final_regridded}. Saltando..."
    continue
  fi

  # === 1. UNIR HISTORICAL 1960-2014 ===
  echo "Merging historical files (1960–2014)..."
  cdo mergetime \
    ${dir_hist}/${var}_day_${model}_historical_${realization}_${grid}_19500101-19741231.nc \
    ${dir_hist}/${var}_day_${model}_historical_${realization}_${grid}_19750101-19991231.nc \
    ${dir_hist}/${var}_day_${model}_historical_${realization}_${grid}_20000101-20141231.nc \
    ${hist_merged}

  # === 2. RECORTAR SSP245 a 2015–2024 ===
  echo "Selecting 2015–2024 from SSP245..."
  ssp_file="${dir_ssp_base}/${realization}/day/${var}/${grid}/${var}_day_${model}_ssp245_${realization}_${grid}_20150101-20391231.nc"
  cdo selyear,2015/2024 ${ssp_file} ${ssp_merged}

  # === 3. CONCATENAR HISTORICAL + SSP245 ===
  echo "Concatenating historical + SSP245..."
  cdo mergetime ${hist_merged} ${ssp_merged} ${final_merged}

  # === 4. REGRILLEAR ===
  echo "Regridding to match landsea mask grid..."
  cdo remapbil,${mask_file} ${final_merged} ${final_regridded}

  # === 5. LIMPIEZA ===
  echo "Cleaning up temporary files..."
  rm -f ${hist_merged} ${ssp_merged} ${final_merged}

  echo "✅ Finalizado: ${final_regridded}"
done
