Anablepsoides cearensis — Ecological Niche Modeling and Climate Projections
This repository contains R scripts, occurrence records, spatial data, and supporting results associated with ecological niche modeling, distribution assessment, and future climate projections for Anablepsoides cearensis.

1. Repository structure

- 01_environmental_selection.R

Evaluates environmental predictor collinearity using Pearson correlation coefficients and variance inflation factors (VIF).

- 02_distribution_models.R

Calibrates and evaluates MaxEnt models using three successive occurrence datasets, corresponding to historical records and two stages of field surveys.

- 03_current_future_projections.R

Calibrates the species-specific MaxEnt model using A. cearensis occurrences, selects the optimal configuration, generates current and future habitat suitability predictions, and produces ensemble projections across global climate models.

Supplementary MaxEnt runs are used to inspect diagnostic outputs, including the 10th-percentile training presence threshold.

- 04_tss_kappa_evaluation.R

Calculates the True Skill Statistic (TSS) and Cohen's Kappa for the selected species-specific model using presence records and background points treated as pseudoabsences.

# The data/ directory contains the occurrence records and spatial files required by the analyses. 

# The Results/ directory contains supporting model outputs and evaluation results.

2. Environmental data sources

2.1. Current climate data

Current environmental variables were obtained from WorldClim version 2.1 (Fick & Hijmans, 2017), representing average climatic conditions for the period 1970–2000, at a spatial resolution of 30 arc-seconds (approximately 1 km at the equator).

The original dataset includes 19 bioclimatic variables (BIO1–BIO19). Elevation data provided by WorldClim 2.1, were also considered during the initial environmental variable assessment.

Environmental predictors were assessed for collinearity using Pearson correlation coefficients and variance inflation factors (VIF), followed by preliminary model evaluation.
Five bioclimatic variables were retained for the final ecological niche models:

- BIO4: Temperature Seasonality
- BIO7: Temperature Annual Range
- BIO15: Precipitation Seasonality
- BIO16: Precipitation of Wettest Quarter
- BIO18: Precipitation of Warmest Quarter

Official source: https://www.worldclim.org/data/worldclim21.html

Local data directory: data/worldclim_alt_present/

2.2. Future climate data

Future bioclimatic variables were obtained from WorldClim version 2.1, using downscaled climate projections from the Coupled Model Intercomparison Project Phase 6 (CMIP6), at a spatial resolution of 30 arc-seconds.
Three global climate models (GCMs) were used:

- BCC-CSM2-MR
- CNRM-ESM2-1
- MIROC6

Climate projections were evaluated under three Shared Socioeconomic Pathways (SSPs):

- SSP2-4.5
- SSP3-7.0
- SSP5-8.5

Two future periods were considered: 2041–2060 and 2081–2100.

This resulted in 18 climate datasets, representing all combinations of three GCMs, three SSPs, and two time periods.

For each dataset, the same five bioclimatic predictors used in the current models (BIO4, BIO7, BIO15, BIO16, and BIO18) were selected.

Official source: https://www.worldclim.org/data/cmip6/cmip6_clim30s.html

Local data directory: data/cmip6/

3. Data availability and reproducibility

The original environmental raster datasets are not included in this GitHub repository because of their large file sizes and redistribution considerations.

Users wishing to reproduce the analyses should download the corresponding current and future bioclimatic datasets directly from WorldClim using the links provided above.

Downloaded datasets should be placed in the indicated directories, preserving the file and folder names expected by the R scripts.

4. Running the analyses

After downloading and extracting the datasets:

- Place the current bioclimatic and elevation GeoTIFF files in data/worldclim_alt_present/.

- Place future climate GeoTIFF files in their respective subdirectories within data/cmip6/, organized by GCM, SSP, and time period.

- Preserve the original file and directory names expected by the R scripts.


The recommended workflow is:

1. Run 01_environmental_selection.R to examine predictor collinearity.

2. Run 02_distribution_models.R to calibrate the three distribution models.

3. Run 03_current_future_projections.R to calibrate the species-specific model and generate future climate projections.

4. Run 04_tss_kappa_evaluation.R after script 03, within the same R session, using the model objects created during calibration.

Scripts may generate or overwrite files in the Results/ directory.

4. Reference

Fick, S. E., & Hijmans, R. J. (2017). WorldClim 2: New 1-km spatial resolution climate surfaces for global land areas. International Journal of Climatology, 37(12), 4302–4315. https://doi.org/10.1002/joc.5086