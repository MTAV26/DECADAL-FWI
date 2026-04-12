# DECADAL-FWI
#============================================================================================================================================================================
#============================================================================================================================================================================
                                        
                                                                    Multi-year predictions of European extreme fire weather conditions.

#============================================================================================================================================================================
#============================================================================================================================================================================

                                                 Miguel Ángel Torres-Vázquez^1*, Marco Turco^2, Carlos Delgado-Torres^3, Francesca Di Giuseppe^4,
                                                         Panos J Athanasiadis^5, Leone Cavicchia^5, Dario Nicolì^5, Enrico Scoccimarro^5

(1)	Universidad de Alcalá, Environmental Remote Sensing Research Group, Department of Geology, Geography and the Environment, Calle Colegios 2, Alcalá de Henares 28801, Spain.
(2)	Regional Atmospheric Modelling (MAR) Group, Department of Physics, Regional Campus of International Excellence Campus Mare Nostrum (CEIR), University of Murcia, 30100 Murcia, Spain.
(3)	Barcelona Supercomputing Center (BSC), Barcelona, Spain.
(4)	European Center for Medium-range Weather Forecast (ECMWF), Reading, UK
(5)	CMCC Foundation - Euro-Mediterranean Center on Climate Change, Italy.



                                                               Corresponding author: Miguel Ángel Torres-Vázquez (miguela.torres@uah.es)

#===================================================================================================================================================================================


OVERVIEW
This repository contains the R scripts used to compute Fire Weather Index (FWI)-based products and to evaluate deterministic and probabilistic skill over Europe using initialized (INIT) and uninitialized/historical (NO-INIT or HIST+245) simulations.

From the scripts available in this repository, the workflow can be reconstructed at two main levels:

1. Production of daily FWI fields from meteorological inputs.
2. Validation of smoothed FWI-derived predictors using correlation-based deterministic metrics and RPSS-based probabilistic metrics.

FILES INCLUDED
1. calculate_FWI.R
   Computes FWI for initialized decadal simulations by looping over daily tas, hurs, pr and sfcWind files and exporting one NetCDF per run.

2. calculate_FWI_ERA5.R
   Same logic and structure as calculate_FWI.R. It appears to be another version of the initialized FWI production script and may be redundant.

3. calculate_FWI_NO-INIT.R
   Computes FWI for the 10 historical/uninitialized CMCC-CM2-SR5 members (r2i1p2f1 to r11i1p2f1).

4. calculate_FWI_daily_ERA5.R
   Computes daily FWI from ERA5 proxy data remapped to the target grid. This script is useful as a template for the observational reference, but it is currently incomplete because the variable `wind_file` is not assigned.

5. CorrEno.R
   Custom function for ensemble-mean correlation with significance estimated using effective sample size (via `s2dv::Eno`).

6. B1_deterministic_validation_COR.R
   Deterministic validation script. Produces four correlation-based European maps and a CSV with domain diagnostics.

7. B2_probabilistic_validation_RPSS.R
   Probabilistic validation script. Produces two RPSS maps and prints diagnostic summaries for full domain and plotted bounding box.

SCIENTIFIC LOGIC OF THE WORKFLOW
The scientific workflow implied by the scripts is the following:

Step 1. Build daily FWI fields.
- ERA5-based observational or proxy FWI is computed from temperature, relative humidity, precipitation and wind.
- Initialized model-member FWI is computed for each available run.
- Historical/uninitialized model-member FWI is computed for 10 members.

Step 2. Build comparable derived FWI products for verification.
The validation scripts do not use the raw daily FWI fields directly. Instead, they expect already-prepared annual smoothed products in NetCDF format, labelled as `FWI95d` and `runmean5`.

Step 3. Deterministic skill assessment.
`B1_deterministic_validation_COR.R` computes:
- correlation INIT vs OBS
- correlation HIST+245 vs OBS
- residual correlation
- difference in correlation (INIT - HIST+245)

Step 4. Probabilistic skill assessment.
`B2_probabilistic_validation_RPSS.R` computes:
- RPSS of INIT against climatology
- RPSS of INIT against HIST+245

Step 5. Figures and diagnostics.
The validation scripts export publication-style PDFs in Robinson projection and, for the deterministic branch, an additional CSV summary with domain statistics.

SOFTWARE REQUIREMENTS
Recommended environment:
- R >= 4.2
- Linux recommended, although parts of the code were also edited under Windows

R packages required across the repository:
- loadeR
- transformeR
- visualizeR
- loadeR.2nc
- fireDanger
- ncdf4
- abind
- s2dv
- multiApply
- fields
- sf
- ggplot2
- reshape2
- scales
- RColorBrewer
- rnaturalearth
- grid
- maps

Example installation:

install.packages(c(
  "ncdf4", "abind", "fields", "sf", "ggplot2", "reshape2",
  "scales", "RColorBrewer", "rnaturalearth", "maps"
))

# Depending on your setup, the following may need specific repositories
# or manual installation:
# loadeR, transformeR, visualizeR, loadeR.2nc, fireDanger, s2dv, multiApply

EXPECTED DIRECTORY STRUCTURE
The scripts use hard-coded absolute paths. Before running them, edit the path variables so they match your system.

A reproducible directory layout would be:

project_root/
├── scripts/
│   ├── calculate_FWI.R
│   ├── calculate_FWI_ERA5.R
│   ├── calculate_FWI_NO-INIT.R
│   ├── calculate_FWI_daily_ERA5.R
│   ├── CorrEno.R
│   ├── B1_deterministic_validation_COR.R
│   └── B2_probabilistic_validation_RPSS.R
├── data/
│   ├── landsea_mask.nc
│   ├── ERA5/
│   │   ├── tasmean/
│   │   ├── hurs/
│   │   ├── pr/
│   │   └── wind/
│   ├── INIT/
│   │   ├── tas_day/
│   │   ├── hurs_day/
│   │   ├── pr_day/
│   │   └── sfcWind_day/
│   ├── NO_INIT/
│   │   ├── tas/
│   │   ├── hurs/
│   │   ├── pr/
│   │   └── sfcWind/
│   └── runmean5/
│       ├── FWI_obs_FWI95d_1961-2024_fixed_runmean5.nc
│       ├── FWI95d_reconstructed_<member>_1961-2024_runmean5.nc
│       └── FWI_CMCC-CM2-SR5_<member>_1960_2024_FWI95d_runmean5.nc
└── outputs/
    ├── FWI/
    └── plots/

MINIMUM INPUT DATA REQUIRED
1. Land-sea mask
- landsea_mask.nc

2. Initialized model daily meteorological fields
- tas_day_*.nc
- hurs_day_*.nc
- pr_day_*.nc
- sfcWind_day_*.nc

3. Historical/uninitialized daily meteorological fields
- tas_merged_CMCC-CM2-SR5_<member>_19600101-20241231_regridded.nc
- hurs_merged_CMCC-CM2-SR5_<member>_19600101-20241231_regridded.nc
- pr_merged_CMCC-CM2-SR5_<member>_19600101-20241231_regridded.nc
- sfcWind_merged_CMCC-CM2-SR5_<member>_19600101-20241231_regridded.nc

4. ERA5 reference inputs
- tasmean_remapped.nc
- rhmean_remapped.nc
- precip_remapped.nc
- remapped wind file corresponding to the variable used in `calculate_FWI_daily_ERA5.R`

5. Final validation-ready NetCDF files
These are required by B1 and B2 and are not generated directly by the uploaded validation scripts:
- observed FWI95d runmean5 file
- reconstructed INIT FWI95d runmean5 files
- historical/uninitialized FWI95d runmean5 files

HOW TO RUN THE ANALYSIS
A. Compute initialized-member FWI
Edit the paths in `calculate_FWI.R` and then run:

Rscript calculate_FWI.R

Expected behavior:
- loops through available initialized tas files
- checks whether matching hurs, pr and sfcWind files exist
- logs missing files
- computes FWI with `makeMultiGrid()` + `fwiGrid()`
- writes one NetCDF output per run

B. Compute NO-INIT / HIST+245 FWI
Edit the paths in `calculate_FWI_NO-INIT.R` and run:

Rscript calculate_FWI_NO-INIT.R

Expected behavior:
- loops over members r2i1p2f1 to r11i1p2f1
- computes daily FWI
- exports one NetCDF per member

C. Compute observational/reference FWI from ERA5
Preferred script template:
- `calculate_FWI_daily_ERA5.R`

Then run:

Rscript calculate_FWI_daily_ERA5.R

D. Build validation-ready products
This step is required but not fully documented in the uploaded scripts.
You must derive the `FWI95d` and `runmean5` NetCDF products used by B1 and B2. In practice, this means creating temporally aligned files for observations, initialized reconstructions and historical members on the same lon-lat-time grid.

E. Run deterministic validation
Edit `obs_dir`, `out_dir`, and ensure `source("CorrEno.R")` points correctly.
Then run:

Rscript B1_deterministic_validation_COR.R

Expected outputs:
- map_skill_init_robinson_steps_landonly.pdf
- map_skill_hist245_robinson_steps_landonly.pdf
- map_residual_skill_robinson_steps_landonly.pdf
- map_added_value_robinson_steps_landonly.pdf
- debug_domain_full_vs_bbox.csv

F. Run probabilistic validation
Edit `obs_dir`, `out_dir`, and verify data alignment.
Then run:

Rscript B2_probabilistic_validation_RPSS.R

Expected outputs:
- map_rpss_climatology_robinson_steps.pdf
- map_rpss_hist_robinson_steps.pdf
- diagnostic tables printed to console
- optional CSV export is prepared in the script but commented out

METHOD NOTES
1. CorrEno
`CorrEno.R` computes the ensemble mean first, then estimates Pearson correlation and significance using an effective sample size correction based on `s2dv::Eno`.

2. Deterministic validation
The deterministic script uses:
- `CorrEno()` for INIT vs OBS and HIST+245 vs OBS
- `s2dv::ResidualCorr()` for residual skill
- `s2dv::DiffCorr()` for added value (difference in correlation)

3. Probabilistic validation
The RPSS script uses tercile thresholds:
- `prob_thresholds = c(1/3, 2/3)`
and defines climatology over indices 31:60.

4. Plotting domain
The maps are cropped to a European window:
- longitude: -15 to 45
- latitude: 30 to 72


REPRODUCIBILITY STATEMENT
With the scripts currently available, a reader can reproduce the main validation logic, the figure generation, and much of the FWI production workflow, provided they already have access to the underlying meteorological NetCDF inputs and to the missing intermediate preprocessing that generates the `FWI95d/runmean5` validation files.

If strict full reproducibility is required for publication, the missing preprocessing code and final data-assembly step should be added to the repository.
