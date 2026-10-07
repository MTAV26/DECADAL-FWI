# DECADAL-FWI
#========================
                                        
Multi-year predictions of European extreme fire weather conditions.

#==========================

Miguel Ángel Torres-Vázquez^1*, Marco Turco^2, Carlos Delgado-Torres^3, Francesca Di Giuseppe^4,
Panos J Athanasiadis^5, Leone Cavicchia^5, Dario Nicolì^5, Enrico Scoccimarro^5

(1)	Universidad de Alcalá, Environmental Remote Sensing Research Group, Department of Geology, Geography and the Environment, Calle Colegios 2, Alcalá de Henares 28801, Spain.
(2)	Regional Atmospheric Modelling (MAR) Group, Department of Physics, Regional Campus of International Excellence Campus Mare Nostrum (CEIR), University of Murcia, 30100 Murcia, Spain.
(3)	Barcelona Supercomputing Center (BSC), Barcelona, Spain.
(4)	European Center for Medium-range Weather Forecast (ECMWF), Reading, UK
(5)	CMCC Foundation - Euro-Mediterranean Center on Climate Change, Italy.



                                                               Corresponding author: Miguel Ángel Torres-Vázquez (miguela.torres@uah.es)
                                                               Contact: Marco Turco, University of Murcia (marco.turco@um.es)

#################################################################################################

#################################################################################################
# A. General instructions
#################################################################################################

The workflow uses shell, Python and R scripts. Run scripts in numerical order.

Project structure:

  scripts/
    1_data_preparation/   data preparation and verification files
    2_analysis/           skill, detrending and trend analyses
    3_figures/            figures only; no scientific calculations
    functions/            shared R functions
    utilities/            auxiliary checks

The study period is 1961-2024. Verification is performed on 60 overlapping FY1-FY5
five-year windows, labelled by their central years 1963-2022.

Initialized predictions use the mean over forecast years 1-5 (FY1-FY5).
FWI95d is the annual number of days above the local 1991-2020 FWI 95th percentile.


For detrended analyses:
  - ERA5/OBS is linearly detrended independently at each grid point.
  - For INIT and HIST+245, the linear trend estimated from the corresponding ensemble
    mean is removed from every ensemble member before verification.
  - This preserves the ensemble spread while removing the common linear trend.

Raw ERA5 and model data are not distributed because of their size. Prepared verification
files are sufficient to run the analysis and figure scripts.

#################################################################################################
# B. Data preparation
#################################################################################################

Main preparation scripts:

  1_01-1_08   ERA5 inputs, daily FWI and observed FWI95d
  1_09-1_19   initialized and uninitialized model processing
  1_20        common FWI95d verification series
  1_21        meteorological-driver verification series

The prepared FWI95d verification data are stored under:

  data/verification_new/

The prepared meteorological-driver verification data are stored under:

  data/verification_drivers_new/

#################################################################################################
# C. Analysis workflow
#################################################################################################

Run from <PROJECT_ROOT>/scripts/:

  Rscript 2_analysis/2_01_calculate_deterministic_skill.R
  Rscript 2_analysis/2_02_calculate_probabilistic_skill.R
  Rscript 2_analysis/2_03_calculate_detrended_skill.R
  Rscript 2_analysis/2_04_calculate_driver_skill.R
  Rscript 2_analysis/2_05_calculate_fwi95d_trends.R
  Rscript 2_analysis/2_06_calculate_driver_skill_detrended.R

Analysis scripts:

  2_01   Raw deterministic skill
         - INIT vs ERA5 Pearson correlation
         - residual correlation relative to HIST+245
         - difference in correlation: INIT minus HIST+245

  2_02   Raw probabilistic skill
         - tercile RPSS vs climatology
         - tercile RPSS vs HIST+245

  2_03   Detrended verification
         - detrended INIT vs ERA5 correlation
         - detrended residual correlation
         - detrended difference in correlation
         - detrended RPSS vs climatology
         - detrended RPSS vs HIST+245

  2_04   Raw deterministic skill of the four FWI meteorological drivers:
         temperature, relative humidity, precipitation and wind speed

  2_05   Linear trends in five-year mean FWI95d for ERA5, INIT and HIST+245

  2_06   Detrended deterministic skill of the four FWI meteorological drivers


#################################################################################################
# D. Figure workflow
#################################################################################################

Figure scripts READ existing analysis outputs and DO NOT recalculate skill.

Run from <PROJECT_ROOT>/scripts/:

  Rscript 3_figures/3_01_plot_deterministic_skill_robinson.R
  Rscript 3_figures/3_02_plot_detrended_skill_robinson.R
  Rscript 3_figures/3_03_plot_rpss_robinson.R
  Rscript 3_figures/3_04_plot_driver_skill_robinson.R
  Rscript 3_figures/3_05_plot_fwi95d_trends_robinson.R
  Rscript 3_figures/3_06_plot_driver_skill_detrended_robinson.R
  Rscript 3_figures/3_07_plot_rpss_detrended_robinson.R

Current figure mapping:

Main manuscript
  3_01 -> Figure 2: raw deterministic FWI95d skill
  3_02 -> Figure 3: detrended deterministic FWI95d skill
  3_03 -> Figure 4: raw probabilistic FWI95d skill (RPSS)

Supplementary Material
  3_04 -> Figure S1: raw deterministic skill of the FWI meteorological drivers
  3_05 -> Figure S2: linear trends in FWI95d for ERA5, INIT and HIST+245
  3_06 -> Figure S3: detrended deterministic skill of the FWI meteorological drivers
  3_07 -> Figure S4: detrended RPSS vs climatology and HIST+245
