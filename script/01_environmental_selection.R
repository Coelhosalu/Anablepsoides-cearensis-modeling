# -------------------------------------------------------------------------
# Pearson correlation analysis and environmental variable selection
# -------------------------------------------------------------------------

# Install and load required packages --------------------------------------

if(!require(raster)){install.packages("raster")}
library(raster)

# Load environmental raster layers ----------------------------------------


sta <- stack(list.files(path = "./data/worldclim_alt_present", 
                        pattern = ".tif$", full.names = TRUE))
names(sta)

# Load the shapefile defining the study area (freshwater ecoregions)
mapi.nemo <- raster::shapefile("./data/shapefiles/mapi_nemo.shp")
plot(mapi.nemo)

# Define the geographic extent of Brazil
extent_brasil <- extent(-75.05, -33.75, -34.80, 5.25)

# Crop environmental raster layers to the geographic extent of Brazil
bio_cropped <- crop(sta, extent_brasil)

# Mask and crop environmental layers to the study area
envs.ne <- raster::mask(bio_cropped, mapi.nemo)
envs.ne <- raster::crop(envs.ne, mapi.nemo)
plot(envs.ne[[19]])
names(envs.ne)


# Pearson correlation analysis --------------------------------------------
cor <- layerStats(envs.ne, 'pearson', na.rm=T)


# Export the Pearson correlation matrix
write.csv(cor$`pearson correlation coefficient`, file = 'Results/distribution/Dis-Pearson-Correlationship.csv', row.names = FALSE)

# Variance Inflation Factor (VIF) analysis ---------------------------------

# Select variables for VIF analysis following inspection of Pearson correlations
names(envs.ne)

# 14, 15, 16, 18, 19, 3, 4, 7, 9
#sta.selec <- envs.ne[[c(6,7,8,10,11,13,14,16,18)]]
sta.selec <- envs.ne[[c(1,4,6,7,8,10,11,13,14,16,19)]]

library(usdm)
library(terra)

# Convert RasterStack to SpatRaster for compatibility with usdm
sta.selec.terra <- terra::rast(sta.selec)

# Verify raster object class and layer names
class(sta.selec.terra)
names(sta.selec.terra)

# Set random seed and perform stepwise VIF selection (threshold = 3)
set.seed(1)
vif <- usdm::vifstep(sta.selec.terra, th =3)

# Extract the correlation matrix and inspect VIF results
matrix.cor <- vif@corMatrix
vif


# Final predictor positions used in ecological niche modeling: 7, 8, 10, 14, 16
# BIO19 and elevation were subsequently excluded based on negligible
# variable importance in preliminary MaxEnt calibrations.

# Pearson correlation coefficients among selected predictors:
# Variable pair       Pearson's r
# BIO15 x BIO18       -0.618
# BIO15 x BIO4        -0.446
# BIO16 x BIO4        -0.401
# BIO4  x BIO18        0.380
# BIO7  x BIO18        0.260


final_vars <- c(
  "wc2.1_30s_bio_4",
  "wc2.1_30s_bio_7",
  "wc2.1_30s_bio_15",
  "wc2.1_30s_bio_16",
  "wc2.1_30s_bio_18"
)

# Verify that final predictors passed the VIF filter
stopifnot(all(final_vars %in% vif@results$Variables))
envs.final <- envs.ne[[final_vars]]
names(envs.final)
plot(envs.final)

# Export VIF correlation matrix
write.csv(vif@corMatrix, file = 'Results/distribution/Dis-VIF-Correlationship.csv', row.names = FALSE)

# Export variables excluded due to multicollinearity
write.csv(vif@excluded, file = 'Results/distribution/Dis-VIF-collinearity-problem.csv', row.names = FALSE)

# Export VIF values for retained variables
write.csv(vif@results, file = 'Results/distribution/Dis-VIF-results.csv', row.names = FALSE)
