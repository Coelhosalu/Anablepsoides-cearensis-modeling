# Mudando o diretório de trabalho para o caminho desejado
setwd("~/UFRN/LISE/Projetos/Mestrado - Anablepsoides/Modelagem/Anablepsoides NE models")

load("~/UFRN/LISE/Projetos/Mestrado - Anablepsoides/Modelagem/Anablepsoides NE models/Results/future/current/models_predict.RData")

# ============================================================================
# Anablepsoides cearensis: MaxEnt calibration and climate projections
# ============================================================================

# Create input and output directories
dir.create(path = "data")
dir.create(path = "output")

# Allocate memory for the Java Virtual Machine (JVM)
options(java.parameters = "-Xmx4g")
memory.limit(1000000)

# -------------------------------------------------------------------------
# 1. Load required R packages
# -------------------------------------------------------------------------

install.packages("rlang")
library(rlang)

if (!require(dismo)) install.packages("dismo")
if (!require(maptools)) install.packages("maptools")
if (!require(rgdal)) install.packages("rgdal")
if (!require(sp)) install.packages("sp")
if (!require(ecospat)) install.packages("ecospat")
if (!require(ENMeval)) install.packages("ENMeval")
if (!require(raster)) install.packages("raster")
if (!require(dplyr)) install.packages("dplyr")

library(dismo)
library(maptools)
library(rgdal)
library(sp)
library(ecospat)
library(ENMeval)
library(raster)
library(dplyr)
library(rasterVis)
library(RColorBrewer)

# -------------------------------------------------------------------------
# 2. Import and filter A. cearensis occurrence records
# -------------------------------------------------------------------------

# Lendo os dados de ocorrências e removendo coordenadas duplicadas -------------
occ_total_raw <- read.csv("./data/occurrences/occ_ce_total.csv", header = T, sep = ";")
occ_total <- occ_total_raw[,2:3]

          # Remove occurrences falling within the same environmental raster cell

          occs.cells <- raster::extract(bio.climatic[[1]], occ_total, cellnumbers = TRUE)
          occs.cellDups <- duplicated(occs.cells[,1])
          occ_total <- occ_total[!occs.cellDups,]

          # Remove records with identical geographic coordinates
occ_total <- occ_total[!duplicated(occ_total),]

      # Create a spatial object from the original occurrence coordinates
      #coordinates used to detect collinearity variables
      occ_total_coord <- occ_total_raw[,2:3]
      occ_total_coord <- unique(occ_total_coord)
      coordinates(occ_total_coord) <- c("long","lat")
      proj4string(occ_total_coord) <- raster::projection(raster::raster())

# -------------------------------------------------------------------------
# 3. Load present-day WorldClim bioclimatic variables
# -------------------------------------------------------------------------
      
bio.climatic <- raster::stack(list.files(path = "./data/worldclim_alt_present", 
                                         pattern = ".tif$", full.names = TRUE))
names(bio.climatic)

# WorldClim v2.1; spatial resolution: 30 arc-seconds (~1 km at the equator)
# Source: https://www.worldclim.org/data/worldclim21.html

# -------------------------------------------------------------------------
# 4. Define the model calibration area and select predictors
# -------------------------------------------------------------------------

          # Load the spatial mask defining the original calibration area
mapi.nemo <- raster::shapefile("./data/shapefiles/mapi_nemo.shp")
plot(mapi.nemo)

# Predictor selection followed correlation/VIF analyses and preliminary model assessment

#BIO: 15,16,18,4,7,9 (9 excluida por baixa importancia nos modelos)
bio_vif <- bio.climatic[[c(7,8,10,14,16)]]

names(bio_vif)


# Select the five environmental predictors in their original layer order

# Inspect raster resolution and extent
res(bio_vif)
extent(bio_vif)


# Define a geographic extent covering Brazil
extent_brasil <- extent(-75.05, -33.75, -34.80, 5.25)

# Crop the selected environmental rasters to the geographic extent
bio_cropped <- crop(bio_vif, extent_brasil)

# Apply the original calibration-area mask and crop to its extent
envs.ne <- raster::mask(bio_cropped, mapi.nemo)
envs.ne <- raster::crop(envs.ne, mapi.nemo)
plot(envs.ne)

# Export selected present-day predictors for ancillary analyses and interpretation
names(envs.ne)
writeRaster(envs.ne[[1]], filename = "Results/discu/current/bio15", format = "raster", overwrite = TRUE)

writeRaster(envs.ne[[2]], filename = "Results/discu/current/bio16", format = "raster", overwrite = TRUE)

writeRaster(envs.ne[[3]], filename = "Results/discu/current/bio18", format = "raster", overwrite = TRUE)

writeRaster(envs.ne[[4]], filename = "Results/discu/current/bio04", format = "raster", overwrite = TRUE)

writeRaster(envs.ne[[5]], filename = "Results/discu/current/bio07", format = "raster", overwrite = TRUE)

#bio04, bio15, bio16


# -------------------------------------------------------------------------
# 5. Generate background points
# -------------------------------------------------------------------------


set.seed(1) 

# Sample 1,000 background points from the calibration area
bg.ne <- dismo::randomPoints(envs.ne, n = 1000) %>% as.data.frame()
colnames(bg.ne) <- colnames(occ_total)

# -------------------------------------------------------------------------
# 6. Calibrate candidate MaxEnt models with ENMeval
# -------------------------------------------------------------------------

e.mx.l <- ENMevaluate(occs = occ_total,
                      envs = envs.ne, 
                      bg = bg.ne, 
                      doClamp  = FALSE,
                      algorithm = 'maxent.jar', 
                      partitions = 'jackknife', #método de validação cruzada
                      #argumentos de ajuste do modelo
                      tune.args = list(fc = c("L", "LQ" ,"H", "LQH"),
                                       #feature class combinations
                                       rm = seq(0.5, 4, 0.5) #
                                       #regularization multiplayes 
                      ))


# -------------------------------------------------------------------------
# 7. Select the best-supported model
# -------------------------------------------------------------------------

# Create the directory for current-model outputs

ifelse(!dir.exists(file.path('Results/future/current')), 
       dir.create(file.path('Results/future/current')), FALSE)

# Save the complete ENMeval evaluation table

res <- eval.results(e.mx.l)
write.csv(res, 'Results/future/current/current_mapi.nemo.csv', row.names = FALSE)

res_filter <- res %>%
  filter(!is.na(AICc))

opt.aicc <- res_filter %>% filter(delta.AICc == 0)
opt.aicc

mod.seq <- eval.models(e.mx.l)[[opt.aicc$tune.args]]
plot(mod.seq)
mod.seq
eval.models(e.mx.l)

# select the model predictions for our optimal model

pred.aicc <- eval.predictions(e.mx.l)[[opt.aicc$tune.args]]
plot(pred.aicc)

  eval.models(e.mx.l)[[fc.L_rm.0.5]]
  
# -------------------------------------------------------------------------
# 8. Supplementary MaxEnt runs for threshold assessment
# -------------------------------------------------------------------------
  
  # Run MaxEnt separately to generate diagnostic output files and
  # examine the 10th percentile training presence threshold (10PTP).
  
  # The first run uses default feature and regularization settings.
  # The second run uses the feature class (L) and regularization
  # multiplier (0.5) selected by ENMeval.
  
  # These supplementary runs are not used to generate the main
  # current or future habitat suitability predictions.
  
  me.mx <- maxent(envs.ne, 
                  occ_total,
                         args = c("-J",
                                  "-P",
                                  "randomtestpoints=30", #50
                                  "replicates=20", #10
                                  "replicatetype=bootstrap",
                                  "randomseed",
                                  "removeduplicates",
                                  "doclamp=false", #none
                                  "extrapolate=true", #none
                                  "convergencethreshold=0.00001",
                                  "maximumiterations=10000",
                                  "maximumbackground=10000"),
                         path='Results/future/current/results_default')
  
  me.e.mx <- maxent(envs.ne, 
                    occ_total,
                    args = c("-J",
                             "-P",
                             "betamultiplier= 0.5",
                             "randomtestpoints=30", #50
                             "replicates=20", #10
                             "replicatetype=bootstrap",
                             "randomseed",
                             "removeduplicates",
                             "doclamp=false", #none
                             "extrapolate=true", #none
                             "convergencethreshold=0.00001",
                             "maximumiterations=10000",
                             "maximumbackground=10000"),
                    features=c("l"),
                    path='Results/future/current/results_L.0.5')
  
# Save best model

write.csv(opt.aicc, "Results/future/current/Opt_current.csv", row.names = FALSE)

# Export variable importance for the selected ENMeval model

Importance <- e.mx.l@variable.importance$fc.L_rm.0.5
write.csv(Importance, "Results/future/current/Current_Importance.csv")
Importance

# Save raster

writeRaster(e.mx.l@predictions$fc.L_rm.0.5,
            'Results/future/current/current.tif', overwrite = TRUE)
plot(e.mx.l@predictions$fc.L_rm.0.5)


# -------------------------------------------------------------------------
# 9. Project the ENMeval-selected model to future climates
# -------------------------------------------------------------------------

# Three CMIP6 GCMs: BCC-CSM2-MR, CNRM-ESM2-1, and MIROC6.
# Emission pathways: SSP2-4.5, SSP3-7.0, and SSP5-8.5.
# Periods: 2041-2060 and 2081-2100.


# Função para fazer as projeções com o modelo selecionado
do_projection <- function(model, raster_data) {
  pred <- dismo::predict(model, raster_data)
  return(pred)
}

# Use the ENMeval-selected model for all climate projections (L; RM = 0.5)
selected_model <- e.mx.l@models$fc.L_rm.0.5


# -------------------------------------------------------------------------
# 9.1. SSP2-4.5 climate projections
# -------------------------------------------------------------------------

#BIO: 15,16,18,4,7,9 (9 excluida por baixa importancia nos modelos)

ssp245.bcc.2041.raw <- raster::stack(list.files(path = "./data/cmip6/wc2.1_30s_bioc_BCC-CSM2-MR_ssp245_2041-2060", 
                                                pattern = ".tif$", full.names = TRUE))
names(ssp245.bcc.2041.raw)
ssp245.bcc.2041.raw <- ssp245.bcc.2041.raw[[c(15,16,18,4,7)]]
names(ssp245.bcc.2041.raw)
names(ssp245.bcc.2041.raw) <- names(envs.ne)
names(ssp245.bcc.2041.raw)

ssp245.bcc.2041.raw <- raster::crop(ssp245.bcc.2041.raw, extent_brasil)
ssp245.bcc.2041.raw <- raster::mask(ssp245.bcc.2041.raw, mapi.nemo)
ssp245.bcc.2041.raw <- raster::crop(ssp245.bcc.2041.raw, mapi.nemo)

  # Project SSP245 / BCC-CSM2-MR / 2041 
    proj_ssp245.bcc.2041 <- do_projection(selected_model, ssp245.bcc.2041.raw)
    plot(proj_ssp245.bcc.2041)
    writeRaster(proj_ssp245.bcc.2041, 'Results/future/ssp245/ssp245.bcc.2041.tif', overwrite = TRUE)

# Export BIO4, BIO15, and BIO16 for supplementary interpretation
    
    writeRaster(ssp245.bcc.2041.raw[[4]], 
                filename = "Results/discu/future/ssp245.bcc.2041_bio04", 
                format = "raster", overwrite = TRUE)
    writeRaster(ssp245.bcc.2041.raw[[1]], 
                filename = "Results/discu/future/ssp245.bcc.2041_bio15", 
                format = "raster", overwrite = TRUE)
    writeRaster(ssp245.bcc.2041.raw[[2]], 
                filename = "Results/discu/future/ssp245.bcc.2041_bio16", 
                format = "raster", overwrite = TRUE)
    
   
    # Project SSP245 / BCC-CSM2-MR / 2081
ssp245.bcc.2081.raw <- raster::stack(list.files(path = "./data/cmip6/wc2.1_30s_bioc_BCC-CSM2-MR_ssp245_2081-2100", 
                                                pattern = ".tif$", full.names = TRUE))
names(ssp245.bcc.2081.raw)
ssp245.bcc.2081.raw <- ssp245.bcc.2081.raw[[c(15,16,18,4,7)]]
names(ssp245.bcc.2081.raw)
names(ssp245.bcc.2081.raw) <- names(envs.ne)
names(ssp245.bcc.2081.raw)

ssp245.bcc.2081.raw <- raster::crop(ssp245.bcc.2081.raw, extent_brasil)
ssp245.bcc.2081.raw <- raster::mask(ssp245.bcc.2081.raw, mapi.nemo)
ssp245.bcc.2081.raw <- raster::crop(ssp245.bcc.2081.raw, mapi.nemo)

      # Projetando os rasters para o cenário futuro "ssp245.bcc.2081.raw"
      proj_ssp245.bcc.2081 <- do_projection(selected_model, ssp245.bcc.2081.raw)
      writeRaster(proj_ssp245.bcc.2081, 'Results/future/ssp245/ssp245.bcc.2081.tif', overwrite = TRUE)
      
      # Export BIO4, BIO15, and BIO16 for supplementary interpretation
      
      writeRaster(ssp245.bcc.2081.raw[[4]], 
                  filename = "Results/discu/future/ssp245.bcc.2081_bio04", 
                  format = "raster", overwrite = TRUE)
      writeRaster(ssp245.bcc.2081.raw[[1]], 
                  filename = "Results/discu/future/ssp245.bcc.2081_bio15", 
                  format = "raster", overwrite = TRUE)
      writeRaster(ssp245.bcc.2081.raw[[2]], 
                  filename = "Results/discu/future/ssp245.bcc.2081_bio16", 
                  format = "raster", overwrite = TRUE)
      
      


ssp245.cnrm.2041.raw <- raster::stack(list.files(path = "./data/cmip6/wc2.1_30s_bioc_CNRM-ESM2-1_ssp245_2041-2060", 
                                                 pattern = ".tif$", full.names = TRUE))
names(ssp245.cnrm.2041.raw)
ssp245.cnrm.2041.raw <- ssp245.cnrm.2041.raw[[c(15,16,18,4,7)]]
names(ssp245.cnrm.2041.raw)
names(ssp245.cnrm.2041.raw) <- names(envs.ne)
names(ssp245.cnrm.2041.raw)

ssp245.cnrm.2041.raw <- raster::crop(ssp245.cnrm.2041.raw, extent_brasil)
ssp245.cnrm.2041.raw <- raster::mask(ssp245.cnrm.2041.raw, mapi.nemo)
ssp245.cnrm.2041.raw <- raster::crop(ssp245.cnrm.2041.raw, mapi.nemo)

      # Projetando os rasters para o cenário futuro "ssp245.cnrm.2041.raw"
      proj_ssp245.cnrm.2041 <- do_projection(selected_model, ssp245.cnrm.2041.raw)
      writeRaster(proj_ssp245.cnrm.2041, 'Results/future/ssp245/ssp245.cnrm.2041.tif', overwrite = TRUE)
      
      #salvando para discu  bio04, bio15, bio16
      
      writeRaster(ssp245.cnrm.2041.raw[[4]], 
                  filename = "Results/discu/future/ssp245.cnrm.2041_bio04", 
                  format = "raster", overwrite = TRUE)
      writeRaster(ssp245.cnrm.2041.raw[[1]], 
                  filename = "Results/discu/future/ssp245.cnrm.2041_bio15", 
                  format = "raster", overwrite = TRUE)
      writeRaster(ssp245.cnrm.2041.raw[[2]], 
                  filename = "Results/discu/future/ssp245.cnrm.2041_bio16", 
                  format = "raster", overwrite = TRUE)
      
      

ssp245.cnrm.2081.raw <- raster::stack(list.files(path = "./data/cmip6/wc2.1_30s_bioc_CNRM-ESM2-1_ssp245_2081-2100", 
                                                 pattern = ".tif$", full.names = TRUE))
names(ssp245.cnrm.2081.raw )
ssp245.cnrm.2081.raw  <- ssp245.cnrm.2081.raw [[c(15,16,18,4,7)]]
names(ssp245.cnrm.2081.raw )
names(ssp245.cnrm.2081.raw ) <- names(envs.ne)
names(ssp245.cnrm.2081.raw )

ssp245.cnrm.2081.raw <- raster::crop(ssp245.cnrm.2081.raw, extent_brasil)
ssp245.cnrm.2081.raw <- raster::mask(ssp245.cnrm.2081.raw, mapi.nemo)
ssp245.cnrm.2081.raw <- raster::crop(ssp245.cnrm.2081.raw, mapi.nemo)

      # Projetando os rasters para o cenário futuro "ssp245.cnrm.2081.raw"
      proj_ssp245.cnrm.2081 <- do_projection(selected_model, ssp245.cnrm.2081.raw)
      writeRaster(proj_ssp245.cnrm.2081, 'Results/future/ssp245/ssp245.cnrm.2081.tif', overwrite = TRUE)
      
      
      #salvando para discu  bio04, bio15, bio16
      
      writeRaster(ssp245.cnrm.2081.raw[[4]], 
                  filename = "Results/discu/future/ssp245.cnrm.2081_bio04", 
                  format = "raster", overwrite = TRUE)
      writeRaster(ssp245.cnrm.2081.raw[[1]], 
                  filename = "Results/discu/future/ssp245.cnrm.2081_bio15", 
                  format = "raster", overwrite = TRUE)
      writeRaster(ssp245.cnrm.2081.raw[[2]], 
                  filename = "Results/discu/future/ssp245.cnrm.2081_bio16", 
                  format = "raster", overwrite = TRUE)
      
      
      
ssp245.miroc6.2041.raw <- raster::stack(list.files(path = "./data/cmip6/wc2.1_30s_bioc_MIROC6_ssp245_2041-2060", 
                                                   pattern = ".tif$", full.names = TRUE))

names(ssp245.miroc6.2041.raw)
ssp245.miroc6.2041.raw <- ssp245.miroc6.2041.raw[[c(15,16,18,4,7)]]
names(ssp245.miroc6.2041.raw)
names(ssp245.miroc6.2041.raw) <- names(envs.ne)
names(ssp245.miroc6.2041.raw)

ssp245.miroc6.2041.raw <- raster::crop(ssp245.miroc6.2041.raw, extent_brasil)
ssp245.miroc6.2041.raw <- raster::mask(ssp245.miroc6.2041.raw, mapi.nemo)
ssp245.miroc6.2041.raw <- raster::crop(ssp245.miroc6.2041.raw, mapi.nemo)


      # Projetando os rasters para o cenário futuro "ssp245.miroc6.2041.raw"
      proj_ssp245.miroc6.2041 <- do_projection(selected_model, ssp245.miroc6.2041.raw)
      writeRaster(proj_ssp245.miroc6.2041, 'Results/future/ssp245/ssp245.miroc6.2041.tif', overwrite = TRUE)
      
      
      
      #salvando para discu  bio04, bio15, bio16
      
      writeRaster(ssp245.miroc6.2041.raw[[4]], 
                  filename = "Results/discu/future/ssp245.miroc6.2041_bio04", 
                  format = "raster", overwrite = TRUE)
      writeRaster(ssp245.miroc6.2041.raw[[1]], 
                  filename = "Results/discu/future/ssp245.miroc6.2041_bio15", 
                  format = "raster", overwrite = TRUE)
      writeRaster(ssp245.miroc6.2041.raw[[2]], 
                  filename = "Results/discu/future/ssp245.miroc6.2041_bio16", 
                  format = "raster", overwrite = TRUE)
      
      
      
      

ssp245.miroc6.2081.raw <- raster::stack(list.files(path = "./data/cmip6/wc2.1_30s_bioc_MIROC6_ssp245_2081-2100", 
                                                   pattern = ".tif$", full.names = TRUE))

names(ssp245.miroc6.2081.raw)
ssp245.miroc6.2081.raw <- ssp245.miroc6.2081.raw[[c(15,16,18,4,7)]]
names(ssp245.miroc6.2081.raw)
names(ssp245.miroc6.2081.raw) <- names(envs.ne)
names(ssp245.miroc6.2081.raw)


ssp245.miroc6.2081.raw <- raster::crop(ssp245.miroc6.2081.raw, extent_brasil)
ssp245.miroc6.2081.raw <- raster::mask(ssp245.miroc6.2081.raw, mapi.nemo)
ssp245.miroc6.2081.raw <- raster::crop(ssp245.miroc6.2081.raw, mapi.nemo)


      # Projetando os rasters para o cenário futuro "ssp245.miroc6.2081.raw"
      proj_ssp245.miroc6.2081 <- do_projection(selected_model, ssp245.miroc6.2081.raw)
      writeRaster(proj_ssp245.miroc6.2081, 'Results/future/ssp245/ssp245.miroc6.2081.tif', overwrite = TRUE)
      

      
      
      
      #salvando para discu  bio04, bio15, bio16
      
      writeRaster(ssp245.miroc6.2081.raw[[4]], 
                  filename = "Results/discu/future/ssp245.miroc6.2081_bio04", 
                  format = "raster", overwrite = TRUE)
      writeRaster(ssp245.miroc6.2081.raw[[1]], 
                  filename = "Results/discu/future/ssp245.miroc6.2081_bio15", 
                  format = "raster", overwrite = TRUE)
      writeRaster(ssp245.miroc6.2081.raw[[2]], 
                  filename = "Results/discu/future/ssp245.miroc6.2081_bio16", 
                  format = "raster", overwrite = TRUE)
      
      
      
      
      
      
# -------------------------------------------------------------------------
# 9.2. SSP3-7.0 climate projections
# -------------------------------------------------------------------------
      
      
      
ssp370.bcc.2041.raw <- raster::stack(list.files (path = "./data/cmip6/wc2.1_30s_bioc_BCC-CSM2-MR_ssp370_2041-2060", 
                                                       pattern = ".tif$", 
                                                       full.names = TRUE))
          names(ssp370.bcc.2041.raw)
          ssp370.bcc.2041.raw <- ssp370.bcc.2041.raw[[c(15,16,18,4,7)]]
          names(ssp370.bcc.2041.raw)
          names(ssp370.bcc.2041.raw) <- names(envs.ne)
          names(ssp370.bcc.2041.raw)
          
          
          ssp370.bcc.2041.raw <- raster::crop(ssp370.bcc.2041.raw, extent_brasil)
          ssp370.bcc.2041.raw <- raster::mask(ssp370.bcc.2041.raw, mapi.nemo)
          ssp370.bcc.2041.raw <- raster::crop(ssp370.bcc.2041.raw, mapi.nemo)
          
          
          
          # Projetando os rasters para o cenário futuro "ssp245.miroc6.2081.raw"
          proj_ssp370.bcc.2041.raw <- do_projection(selected_model, ssp370.bcc.2041.raw)
          writeRaster(proj_ssp370.bcc.2041.raw, 'Results/future/ssp370/ssp370.bcc.2041.raw.tif', overwrite = TRUE)

          
          
          
          #salvando para discu  bio04, bio15, bio16
          
          writeRaster(ssp370.bcc.2041.raw[[4]], 
                      filename = "Results/discu/future/ssp370.bcc.2041_bio04", 
                      format = "raster", overwrite = TRUE)
          writeRaster(ssp370.bcc.2041.raw[[1]], 
                      filename = "Results/discu/future/ssp370.bcc.2041_bio15", 
                      format = "raster", overwrite = TRUE)
          writeRaster(ssp370.bcc.2041.raw[[2]], 
                      filename = "Results/discu/future/ssp370.bcc.2041_bio16", 
                      format = "raster", overwrite = TRUE)
          

ssp370.bcc.2081.raw <- raster::stack(list.files (path = "./data/cmip6/wc2.1_30s_bioc_BCC-CSM2-MR_ssp370_2081-2100", 
                                                       pattern = ".tif$", 
                                                       full.names = TRUE)) 
          names(ssp370.bcc.2081.raw)
          ssp370.bcc.2081.raw <- ssp370.bcc.2081.raw[[c(15,16,18,4,7)]]
          names(ssp370.bcc.2081.raw)
          names(ssp370.bcc.2081.raw) <- names(envs.ne)
          names(ssp370.bcc.2081.raw)
          
          
          ssp370.bcc.2081.raw <- raster::crop(ssp370.bcc.2081.raw, extent_brasil)
          ssp370.bcc.2081.raw <- raster::mask(ssp370.bcc.2081.raw, mapi.nemo)
          ssp370.bcc.2081.raw <- raster::crop(ssp370.bcc.2081.raw, mapi.nemo)
          
          
          
          # Projetando os rasters para o cenário futuro "ssp245.miroc6.2081.raw"
          proj_ssp370.bcc.2081.raw <- do_projection(selected_model, ssp370.bcc.2081.raw)
          writeRaster(proj_ssp370.bcc.2081.raw, 'Results/future/ssp370/ssp370.bcc.2081.raw.tif', overwrite = TRUE)



ssp370.cnrm.2041.raw <- raster::stack(list.files (path = "./data/cmip6/wc2.1_30s_bioc_CNRM-ESM2-1_ssp370_2041-2060", 
                                                        pattern = ".tif$", 
                                                        full.names = TRUE))
          names(ssp370.cnrm.2041.raw)
          ssp370.cnrm.2041.raw <- ssp370.cnrm.2041.raw[[c(15,16,18,4,7)]]
          names(ssp370.cnrm.2041.raw)
          names(ssp370.cnrm.2041.raw) <- names(envs.ne)
          names(ssp370.cnrm.2041.raw)
          
          
          ssp370.cnrm.2041.raw <- raster::crop(ssp370.cnrm.2041.raw, extent_brasil)
          ssp370.cnrm.2041.raw <- raster::mask(ssp370.cnrm.2041.raw, mapi.nemo)
          ssp370.cnrm.2041.raw <- raster::crop(ssp370.cnrm.2041.raw, mapi.nemo)
          
          
          
          # Projetando os rasters para o cenário futuro "ssp245.miroc6.2081.raw"
          proj_ssp370.cnrm.2041.raw <- do_projection(selected_model, ssp370.cnrm.2041.raw)
          writeRaster(proj_ssp370.cnrm.2041.raw, 'Results/future/ssp370/ssp370.cnrm.2041.raw.tif', overwrite = TRUE)



ssp370.cnrm.2081.raw <- raster::stack(list.files (path = "./data/cmip6/wc2.1_30s_bioc_CNRM-ESM2-1_ssp370_2081-2100", 
                                                        pattern = ".tif$", 
                                                        full.names = TRUE))
          names(ssp370.cnrm.2081.raw)
          ssp370.cnrm.2081.raw <- ssp370.cnrm.2081.raw[[c(15,16,18,4,7)]]
          names(ssp370.cnrm.2081.raw)
          names(ssp370.cnrm.2081.raw) <- names(envs.ne)
          names(ssp370.cnrm.2081.raw)
          
          
          ssp370.cnrm.2081.raw <- raster::crop(ssp370.cnrm.2081.raw, extent_brasil)
          ssp370.cnrm.2081.raw <- raster::mask(ssp370.cnrm.2081.raw, mapi.nemo)
          ssp370.cnrm.2081.raw <- raster::crop(ssp370.cnrm.2081.raw, mapi.nemo)
          
          
          
          # Projetando os rasters para o cenário futuro "ssp245.miroc6.2081.raw"
          proj_ssp370.cnrm.2081.raw <- do_projection(selected_model, ssp370.cnrm.2081.raw)
          writeRaster(proj_ssp370.cnrm.2081.raw, 'Results/future/ssp370/ssp370.cnrm.2081.raw.raw.tif', overwrite = TRUE)



ssp370.miroc6.2041.raw <- raster::stack(list.files (path = "./data/cmip6/wc2.1_30s_bioc_MIROC6_ssp370_2041-2060", 
                                                           pattern = ".tif$", 
                                                           full.names = TRUE))
          names(ssp370.miroc6.2041.raw)
          ssp370.miroc6.2041.raw <- ssp370.miroc6.2041.raw[[c(15,16,18,4,7)]]
          names(ssp370.miroc6.2041.raw)
          names(ssp370.miroc6.2041.raw) <- names(envs.ne)
          names(ssp370.miroc6.2041.raw)
          
          
          ssp370.miroc6.2041.raw <- raster::crop(ssp370.miroc6.2041.raw, extent_brasil)
          ssp370.miroc6.2041.raw <- raster::mask(ssp370.miroc6.2041.raw, mapi.nemo)
          ssp370.miroc6.2041.raw <- raster::crop(ssp370.miroc6.2041.raw, mapi.nemo)
          
          
          
          # Projetando os rasters para o cenário futuro "ssp245.miroc6.2081.raw"
          proj_ssp370.miroc6.2041.raw <- do_projection(selected_model, ssp370.miroc6.2041.raw)
          writeRaster(proj_ssp370.miroc6.2041.raw, 'Results/future/ssp370/ssp370.miroc6.2041.raw.tif', overwrite = TRUE)
          
          

ssp370.miroc6.2081.raw <- raster::stack(list.files (path = "./data/cmip6/wc2.1_30s_bioc_MIROC6_ssp370_2081-2100", 
                                                           pattern = ".tif$", 
                                                           full.names = TRUE)) 
        names(ssp370.miroc6.2081.raw)
        ssp370.miroc6.2081.raw <- ssp370.miroc6.2081.raw[[c(15,16,18,4,7)]]
        names(ssp370.miroc6.2081.raw)
        names(ssp370.miroc6.2081.raw) <- names(envs.ne)
        names(ssp370.miroc6.2081.raw)
        
        
        ssp370.miroc6.2081.raw <- raster::crop(ssp370.miroc6.2081.raw, extent_brasil)
        ssp370.miroc6.2081.raw <- raster::mask(ssp370.miroc6.2081.raw, mapi.nemo)
        ssp370.miroc6.2081.raw <- raster::crop(ssp370.miroc6.2081.raw, mapi.nemo)
        
        
        
        # Projetando os rasters para o cenário futuro "ssp245.miroc6.2081.raw"
        proj_ssp370.miroc6.2081.raw <- do_projection(selected_model, ssp370.miroc6.2081.raw)
        writeRaster(proj_ssp370.miroc6.2081.raw, 'Results/future/ssp370/ssp370.miroc6.2081.raw.tif', overwrite = TRUE)
        
        
        
        
        
# -------------------------------------------------------------------------
# 9.3. SSP5-8.5 climate projections
# -------------------------------------------------------------------------
        
        
        
ssp585.bcc.2041.raw <- raster::stack(list.files (path = "./data/cmip6/wc2.1_30s_bioc_BCC-CSM2-MR_ssp585_2041-2060", 
                                                         pattern = ".tif$", 
                                                         full.names = TRUE)) 
        names(ssp585.bcc.2041.raw)
        ssp585.bcc.2041.raw <- ssp585.bcc.2041.raw[[c(15,16,18,4,7)]]
        names(ssp585.bcc.2041.raw)
        names(ssp585.bcc.2041.raw) <- names(envs.ne)
        names(ssp585.bcc.2041.raw)
        
        
        ssp585.bcc.2041.raw <- raster::crop(ssp585.bcc.2041.raw, extent_brasil)
        ssp585.bcc.2041.raw <- raster::mask(ssp585.bcc.2041.raw, mapi.nemo)
        ssp585.bcc.2041.raw <- raster::crop(ssp585.bcc.2041.raw, mapi.nemo)
        
        
        
        # Projetando os rasters para o cenário futuro 
        proj_ssp585.bcc.2041.raw <- do_projection(selected_model, ssp585.bcc.2041.raw)
        writeRaster(proj_ssp585.bcc.2041.raw, 'Results/future/ssp585/ssp585.bcc.2041.raw.tif', overwrite = TRUE)


ssp585.bcc.2081.raw <- raster::stack(list.files (path = "./data/cmip6/wc2.1_30s_bioc_BCC-CSM2-MR_ssp585_2081-2100", 
                                                         pattern = ".tif$", 
                                                         full.names = TRUE)) 
        names(ssp585.bcc.2081.raw)
        ssp585.bcc.2081.raw <- ssp585.bcc.2081.raw[[c(15,16,18,4,7)]]
        names(ssp585.bcc.2081.raw)
        names(ssp585.bcc.2081.raw) <- names(envs.ne)
        names(ssp585.bcc.2081.raw)
        
        
        ssp585.bcc.2081.raw <- raster::crop(ssp585.bcc.2081.raw, extent_brasil)
        ssp585.bcc.2081.raw <- raster::mask(ssp585.bcc.2081.raw, mapi.nemo)
        ssp585.bcc.2081.raw <- raster::crop(ssp585.bcc.2081.raw, mapi.nemo)
        
        
        
        # Projetando os rasters para o cenário futuro 
        proj_ssp585.bcc.2081.raw <- do_projection(selected_model, ssp585.bcc.2081.raw)
        writeRaster(proj_ssp585.bcc.2081.raw, 'Results/future/ssp585/ssp585.bcc.2081.raw.tif', overwrite = TRUE)
        
        
        

ssp585.cnrm.2041.raw <- raster::stack(list.files (path = "./data/cmip6/wc2.1_30s_bioc_CNRM-ESM2-1_ssp585_2041-2060", 
                                                          pattern = ".tif$", 
                                                          full.names = TRUE))
        names(ssp585.cnrm.2041.raw)
        ssp585.cnrm.2041.raw  <- ssp585.cnrm.2041.raw[[c(15,16,18,4,7)]]
        names(ssp585.cnrm.2041.raw)
        names(ssp585.cnrm.2041.raw) <- names(envs.ne)
        names(ssp585.cnrm.2041.raw)
        
        
        ssp585.cnrm.2041.raw <- raster::crop(ssp585.cnrm.2041.raw, extent_brasil)
        ssp585.cnrm.2041.raw <- raster::mask(ssp585.cnrm.2041.raw, mapi.nemo)
        ssp585.cnrm.2041.raw <- raster::crop(ssp585.cnrm.2041.raw, mapi.nemo)
        
        
        
        # Projetando os rasters para o cenário futuro 
        proj_ssp585.cnrm.2041.raw <- do_projection(selected_model, ssp585.cnrm.2041.raw )
        writeRaster(proj_ssp585.cnrm.2041.raw , 'Results/future/ssp585/ssp585.cnrm.2041.raw.tif', overwrite = TRUE)




ssp585.cnrm.2081.raw <- raster::stack(list.files (path = "./data/cmip6/wc2.1_30s_bioc_CNRM-ESM2-1_ssp585_2081-2100", 
                                                          pattern = ".tif$", 
                                                          full.names = TRUE))
        names(ssp585.cnrm.2081.raw)
        ssp585.cnrm.2081.raw <- ssp585.cnrm.2081.raw[[c(15,16,18,4,7)]]
        names(ssp585.cnrm.2081.raw)
        names(ssp585.cnrm.2081.raw) <- names(envs.ne)
        names(ssp585.cnrm.2081.raw)
        
        
        ssp585.cnrm.2081.raw <- raster::crop(ssp585.cnrm.2081.raw, extent_brasil)
        ssp585.cnrm.2081.raw <- raster::mask(ssp585.cnrm.2081.raw, mapi.nemo)
        ssp585.cnrm.2081.raw <- raster::crop(ssp585.cnrm.2081.raw, mapi.nemo)
        
        
        
# Projetando os rasters para o cenário futuro 
proj_ssp585.cnrm.2081.raw  <- do_projection(selected_model, ssp585.cnrm.2081.raw)
writeRaster(proj_ssp585.cnrm.2081.raw, 'Results/future/ssp585/ssp585.cnrm.2081.raw.tif', overwrite = TRUE)




ssp585.miroc6.2041.raw <- raster::stack(list.files (path = "./data/cmip6/wc2.1_30s_bioc_MIROC6_ssp585_2041-2060", 
                                                            pattern = ".tif$", 
                                                            full.names = TRUE))
        names(ssp585.miroc6.2041.raw)
        ssp585.miroc6.2041.raw <- ssp585.miroc6.2041.raw[[c(15,16,18,4,7)]]
        names(ssp585.miroc6.2041.raw)
        names(ssp585.miroc6.2041.raw) <- names(envs.ne)
        names(ssp585.miroc6.2041.raw)
        
        
        ssp585.miroc6.2041.raw <- raster::crop(ssp585.miroc6.2041.raw, extent_brasil)
        ssp585.miroc6.2041.raw <- raster::mask(ssp585.miroc6.2041.raw, mapi.nemo)
        ssp585.miroc6.2041.raw <- raster::crop(ssp585.miroc6.2041.raw, mapi.nemo)
        
        
        
        # Projetando os rasters para o cenário futuro 
        proj_ssp585.miroc6.2041.raw  <- do_projection(selected_model, ssp585.miroc6.2041.raw)
        writeRaster(proj_ssp585.miroc6.2041.raw, 'Results/future/ssp585/ssp585.miroc6.2041.raw.tif', overwrite = TRUE)
        
        
        
        
ssp585.miroc6.2081.raw <- raster::stack(list.files (path = "./data/cmip6/wc2.1_30s_bioc_MIROC6_ssp585_2081-2100", 
                                                            pattern = ".tif$", 
                                                            full.names = TRUE))
        names(ssp585.miroc6.2081.raw)
        ssp585.miroc6.2081.raw <- ssp585.miroc6.2081.raw[[c(15,16,18,4,7)]]
        names(ssp585.miroc6.2081.raw)
        names(ssp585.miroc6.2081.raw) <- names(envs.ne)
        names(ssp585.miroc6.2081.raw)
        
        
        ssp585.miroc6.2081.raw <- raster::crop(ssp585.miroc6.2081.raw, extent_brasil)
        ssp585.miroc6.2081.raw <- raster::mask(ssp585.miroc6.2081.raw, mapi.nemo)
        ssp585.miroc6.2081.raw <- raster::crop(ssp585.miroc6.2081.raw, mapi.nemo)
        
        
        
        # Projetando os rasters para o cenário futuro 
        proj_ssp585.miroc6.2081.raw  <- do_projection(selected_model, ssp585.miroc6.2081.raw)
        writeRaster(proj_ssp585.miroc6.2081.raw, 'Results/future/ssp585/ssp585.miroc6.2081.raw.tif', overwrite = TRUE)
        
# -------------------------------------------------------------------------
# 10. Calculate multi-GCM ensemble means for each SSP and period
# -------------------------------------------------------------------------
        
        # Ensemble dos modelos SSP 245
        ensemble_245_2041 <- list(
          proj_ssp245.bcc.2041,
          proj_ssp245.cnrm.2041,
          proj_ssp245.miroc6.2041
        )
        ensemble_245_2081 <- list(
          proj_ssp245.bcc.2081,
          proj_ssp245.cnrm.2081,
          proj_ssp245.miroc6.2081
        )
        
        # Ensemble dos modelos SSP 370
        ensemble_370_2041 <- list(
          proj_ssp370.bcc.2041.raw,
          proj_ssp370.cnrm.2041.raw,
          proj_ssp370.miroc6.2041.raw
        )
        ensemble_370_2081 <- list(
          proj_ssp370.bcc.2081.raw,
          proj_ssp370.cnrm.2081.raw,
          proj_ssp370.miroc6.2081.raw
        )
        
        # Ensemble dos modelos SSP 585
        ensemble_585_2041 <- list(
          proj_ssp585.bcc.2041.raw,
          proj_ssp585.cnrm.2041.raw,
          proj_ssp585.miroc6.2041.raw
        )
        ensemble_585_2081 <- list(
          proj_ssp585.bcc.2081.raw,
          proj_ssp585.cnrm.2081.raw,
          proj_ssp585.miroc6.2081.raw
        )
        
        # Função para calcular a média das células dos rasters
        mean_fun <- function(...) {
          mean(c(...), na.rm = TRUE)
        }
        
        # Ensemble SSP245
        ensemble_mean_245_2041 <- do.call(overlay, c(ensemble_245_2041, fun = mean_fun))
        plot(ensemble_mean_245_2041)
        writeRaster(ensemble_mean_245_2041, filename = "Results/future/ssp245/2041", format = "raster", overwrite = TRUE)
        ensemble_mean_245_2081 <- do.call(overlay, c(ensemble_245_2081, fun = mean_fun))
        plot(ensemble_mean_245_2081)
        writeRaster(ensemble_mean_245_2081, filename = "Results/future/ssp245/2081", format = "raster", overwrite = TRUE)
        
        # Ensemble SSP370
        ensemble_mean_370_2041 <- do.call(overlay, c(ensemble_370_2041, fun = mean_fun))
        plot(ensemble_mean_370_2041)
        writeRaster(ensemble_mean_370_2041, filename = "Results/future/ssp370/2041", format = "raster", overwrite = TRUE)
        ensemble_mean_370_2081 <- do.call(overlay, c(ensemble_370_2081, fun = mean_fun))
        plot(ensemble_mean_370_2081)
        writeRaster(ensemble_mean_370_2081, filename = "Results/future/ssp370/2081", format = "raster", overwrite = TRUE)
        
        # Ensemble SSP585
        ensemble_mean_585_2041 <- do.call(overlay, c(ensemble_585_2041, fun = mean_fun))
        plot(ensemble_mean_585_2041)
        writeRaster(ensemble_mean_585_2041, filename = "Results/future/ssp585/2041", format = "raster", overwrite = TRUE)
        ensemble_mean_585_2081 <- do.call(overlay, c(ensemble_585_2081, fun = mean_fun))
        plot(ensemble_mean_585_2081)
        writeRaster(ensemble_mean_585_2081, filename = "Results/future/ssp585/2081", format = "raster", overwrite = TRUE)
        
