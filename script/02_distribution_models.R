
# Carregando os dados do modelo Anablepsoides NE 
#load("~/UFRN/LISE/Projetos/Mestrado - Anablepsoides/Modelagem/Anablepsoides NE models/Results/distribution/dsitribution-models.RData")

# Mudando o diretório de trabalho para o caminho desejado
#setwd("~/UFRN/LISE/Projetos/Mestrado - Anablepsoides/Modelagem/Anablepsoides NE models")

# -----------------------------------------------------------------------------
# Anablepsoides: sequential distribution model calibration and evaluation
# Original analysis: three occurrence datasets, ENMeval and MaxEnt (maxent.jar)
# -----------------------------------------------------------------------------

# 1. Load required packages --------------------------------------------------

# Set Java memory before initializing rJava
options(java.parameters = "-Xmx4g")


install.packages("rlang")
library(rlang)

if (!require(dismo)) install.packages("dismo")
if (!require(biomod2)) install.packages("biomod2")
if (!require(terra)) install.packages("terra")
if (!require(rgdal)) install.packages("rgdal")
if (!require(sp)) install.packages("sp")
if (!require(ecospat)) install.packages("ecospat")
if (!require(ENMeval)) install.packages("ENMeval")
if (!require(raster)) install.packages("raster")
if (!require(dplyr)) install.packages("dplyr")
if (!require(dplyr)) install.packages("rasterVis")
if (!require(rJava)) install.packages("rJava")

library(rJava)
library(dismo)
library(biomod2)
library(terra)
library(ecospat)
library(rgdal)
library(sp)
library(ENMeval)
library(raster)
library(dplyr)
library(rasterVis)
library(RColorBrewer)

# Additional packages should be loaded in the sections where needed.

# 2. Import occurrence data --------------------------------------------------


# Historical occurrence records
      occ_prev_raw <- read.csv("./data/occurrences/occ_cebavi_previous.csv", 
                               header = T, 
                               sep = ";")
      occ_prev <- occ_prev_raw[,2:3]
      
      # Removing occurrences that have the same coordinates
      occ_prev <- occ_prev[!duplicated(occ_prev),]
      
      # Criando um objeto "SpatialPointsDataFrame" com as coordenadas das ocorrências
      occ_prev_coord <- occ_prev_raw[,2:3]
      occ_prev_coord <- unique(occ_prev_coord)
      coordinates(occ_prev_coord) <- c("long","lat")
      proj4string(occ_prev_coord) <- raster::projection(raster::raster())
      
      
      
# Records available after the first field survey.
      occ_1.col_raw <- read.csv("./data/occurrences/occ_cebavi_previous_2.csv", 
                               header = T, 
                               sep = ";")
      occ_1.col <- occ_1.col_raw[,2:3]
      
      # Removing occurrences that have the same coordinates
      occ_1.col <- occ_1.col[!duplicated(occ_1.col),]
      
      # Criando um objeto "SpatialPointsDataFrame" com as coordenadas das ocorrências
      occ_1.col_coord <- occ_1.col_raw[,2:3]
      occ_1.col_coord <- unique(occ_1.col_coord)
      coordinates(occ_1.col_coord) <- c("long","lat")
      proj4string(occ_1.col_coord) <- raster::projection(raster::raster())
      
      
      
#Records available after the second field survey
      occ_2.col_raw <- read.csv("./data/occurrences/occ_anablepsoides_total.csv", 
                                header = T, 
                                sep = ";")
      occ_2.col <- occ_2.col_raw[,2:3]
      
      # Removing occurrences that have the same coordinates
      occ_2.col <- occ_2.col[!duplicated(occ_2.col),]
      
      # Criando um objeto "SpatialPointsDataFrame" com as coordenadas das ocorrências
      occ_2.col_coord <- occ_2.col_raw[,2:3]
      occ_2.col_coord <- unique(occ_2.col_coord)
      coordinates(occ_2.col_coord) <- c("long","lat")
      proj4string(occ_2.col_coord) <- raster::projection(raster::raster())
      
      
# 3. Import environmental predictors ---------------------------------------
      

      bio.climatic.elev <- raster::stack(
        list.files(
          path = "./data/worldclim_alt_present",
          pattern = "\\.tif$",
          full.names = TRUE
        )
      )
      names(bio.climatic.elev)
      
      
# 4. Select environmental predictors ---------------------------------------
      
      # Variables retained following Pearson correlations, VIF assessment,
      # and preliminary model evaluation (BIO19 and elevation were excluded).
      # Preserve the predictor order from the original index selection (7,8,10,14,16).
      selected_vars <- c(
        "wc2.1_30s_bio_15",
        "wc2.1_30s_bio_16",
        "wc2.1_30s_bio_18",
        "wc2.1_30s_bio_4",
        "wc2.1_30s_bio_7"
      )
      stopifnot(all(selected_vars %in% names(bio.climatic.elev)))
      env.clim.elev.prev <- bio.climatic.elev[[selected_vars]]
      names(env.clim.elev.prev)
      
      
# 5. Define the original calibration region --------------------------------

      # Polygon delimiting the calibration area along the northeastern coast.
      nordeste <- raster::shapefile("./data/shapefiles/nordeste.shp")
      plot(nordeste)
      
      # Inspect geographic metadata; the raster and polygon CRS must agree.
      res(bio.climatic.elev)
      extent(bio.climatic.elev)
      crs(bio.climatic.elev)
      crs(nordeste)
      
      extent_brasil <- extent(-75.05, -33.75, -34.80, 5.25)
      bio_cropped <- crop(env.clim.elev.prev, extent_brasil)
      envs.ne <- raster::mask(bio_cropped, nordeste)
      envs.ne <- raster::crop(envs.ne, nordeste)
      plot(envs.ne)
      
# 6. Generate background points --------------------------------------------

      set.seed(1)    
      bg <- dismo::randomPoints(envs.ne, n = 1000) %>% as.data.frame()
      colnames(bg) <- colnames(occ_prev)

# 7. Calibrate sequential MaxEnt models ------------------------------------

      # feature classes: L, LQ, H, LQH.
      # regularization multipliers: 0.5 to 5.0, in increments of 0.5.
      # jackknife partitioning.
      
      # Historical occurrences
      e.mx.l <- ENMevaluate(occs = occ_prev, 
                                      envs = envs.ne, 
                                      bg = bg, 
                                      doClamp = FALSE,
                                      algorithm = 'maxent.jar', 
                                      partitions = 'jackknife',
                                      tune.args = list(fc = c("L", "LQ" ,"H", "LQH"),
                                                       rm = seq(0.5, 5, 0.5)))
      # After the first field survey
      e.mx.2 <- ENMevaluate(occs = occ_1.col, 
                            envs = envs.ne, 
                            bg = bg, 
                            doClamp = FALSE,
                            algorithm = 'maxent.jar', 
                            partitions = 'jackknife',
                            tune.args = list(fc = c("L", "LQ" ,"H", "LQH"),
                                             rm = seq(0.5, 5, 0.5)))
      # After the second field survey
      e.mx.3 <- ENMevaluate(occs = occ_2.col, 
                            envs = envs.ne, 
                            bg = bg, 
                            doClamp = FALSE,
                            algorithm = 'maxent.jar', 
                            partitions = 'jackknife',
                            tune.args = list(fc = c("L", "LQ" ,"H", "LQH"),
                                             rm = seq(0.5, 5, 0.5)))

# 8. Select the models -----------------------------------------------------   
      #1
        res1 <- eval.results(e.mx.l)
        res_filter1 <- res1 %>% filter(!is.na(AICc))
        delta_AICc_filter1 <- res_filter1 %>% filter(delta.AICc == 0)
        
        # auc.val.avg and or.10p.avg
        opt_seq1 <- delta_AICc_filter1 %>%
          filter(auc.val.avg == max(auc.val.avg)) %>%
          filter(or.10p.avg == min(or.10p.avg))
        melhor_modelo1 <- eval.models(e.mx.l)[[opt_seq1$tune.args]]

      #2
        res2 <- eval.results(e.mx.2)
        res_filter2 <- res2 %>% filter(!is.na(AICc))
        delta_AICc_filter2 <- res_filter2 %>% filter(delta.AICc == 0)
        
        # auc.val.avg and or.10p.avg
        opt_seq2 <- delta_AICc_filter2 %>%
          filter(auc.val.avg == max(auc.val.avg)) %>%
          filter(or.10p.avg == min(or.10p.avg))
        melhor_modelo2 <- eval.models(e.mx.2)[[opt_seq2$tune.args]]

      #3
        res3 <- eval.results(e.mx.3)
        res_filter3 <- res3 %>% filter(!is.na(AICc))
        delta_AICc_filter3 <- res_filter3 %>% filter(delta.AICc == 0)
        
        # auc.val.avg and or.10p.avg
        opt_seq3 <- delta_AICc_filter3 %>%
          filter(auc.val.avg == max(auc.val.avg)) %>%
          filter(or.10p.avg == min(or.10p.avg))
        melhor_modelo3 <- eval.models(e.mx.3)[[opt_seq3$tune.args]]
        
        
# 9. Export calibration results ---------------------------------------------

        # Create directories for model outputs
        dir.create("Results/distribution/1", recursive = TRUE, showWarnings = FALSE)
        dir.create("Results/distribution/2", recursive = TRUE, showWarnings = FALSE)
        dir.create("Results/distribution/3", recursive = TRUE, showWarnings = FALSE)

        write.csv(res1, 'Results/distribution/1/all_results_prev1.csv', row.names = FALSE)
        write.csv(res2, 'Results/distribution/2/all_results_prev2.csv', row.names = FALSE)
        write.csv(res3, 'Results/distribution/3/all_results_prev3.csv', row.names = FALSE)

# select the model predictions for our optimal model
pred.seq1 <- e.mx.l@predictions$fc.LQ_rm.0.5
writeRaster(pred.seq1 , 'Results/distribution/1/pred.01.tif', overwrite = TRUE)

pred.seq2 <- e.mx.2@predictions$fc.LQH_rm.3
writeRaster(pred.seq2 , 'Results/distribution/2/pred.02.tif', overwrite = TRUE)

pred.seq3 <- e.mx.3@predictions$fc.LQ_rm.0.5
writeRaster(pred.seq3 , 'Results/distribution/3/pred.03.tif', overwrite = TRUE)

plot(pred.seq3)

# Export variable importance for the selected models
Importance1 <- e.mx.l@variable.importance$fc.LQ_rm.0.5
write.csv(Importance1, "Results/distribution/1/Importance1.csv")
Importance1

Importance2 <- e.mx.2@variable.importance$fc.LQH_rm.3
write.csv(Importance2, "Results/distribution/2/Importance2.csv")
Importance2

Importance3 <- e.mx.3@variable.importance$fc.LQ_rm.0.5
write.csv(Importance3, "Results/distribution/3/Importance3.csv")
Importance3

