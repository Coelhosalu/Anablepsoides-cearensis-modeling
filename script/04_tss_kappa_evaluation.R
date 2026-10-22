# -------------------------------------------------------------------------
# Model evaluation: True Skill Statistic (TSS) and Cohen's Kappa
# -------------------------------------------------------------------------


# This script must be executed after the current model calibration script.
# Required objects: occ_total, bg.ne, and pred.aicc.




# Load required packages
library(raster)
library(dplyr)
library(ecospat)

# Combine presence and background coordinates
occ_tss <- rbind(occ_total, bg.ne)

# Extract predicted suitability values at sampling points
pred.tss <- raster::extract(pred.aicc, occ_tss)
pred.tss

# Assign observed values: 1 = presence, 0 = background
occ_total.tss <- occ_total %>% mutate(nova_coluna = 1)
bg.ne.tss <- bg.ne %>% mutate(nova_coluna = 0)

# Combine presence and background observations
occ_tss <- rbind(occ_total.tss,bg.ne.tss)

occ_tss.presence <- occ_tss$nova_coluna

# -------------------------------------------------------------------------
# True Skill Statistic (TSS)
# -------------------------------------------------------------------------

# Calculate TSS across classification thresholds
TSS.values <- ecospat.max.tss(pred.tss, occ_tss.presence)
TSS.values

# Export evaluation table, maximum TSS, and optimal threshold
write.csv(TSS.values$table,
          "Results/future/current/max_tss_table.csv",
          row.names = FALSE)

write.csv(TSS.values$max.TSS,
          "Results/future/current/max_tss.csv",
          row.names = FALSE)

write.csv(TSS.values$max.threshold,
          "Results/future/current/max_tss_threshold.csv",
          row.names = FALSE)

# -------------------------------------------------------------------------
# Cohen's Kappa
# -------------------------------------------------------------------------

# Calculate Kappa across classification thresholds
kappa.values <- ecospat.max.kappa(pred.tss, occ_tss.presence)

# Inspect results
kappa.values

# Export evaluation table, maximum Kappa, and optimal threshold
write.csv(kappa.values$table,
          "Results/future/current/max_kappa_table.csv",
          row.names = FALSE)

write.csv(kappa.values$max.Kappa,
          "Results/future/current/max_kappa.csv",
          row.names = FALSE)

write.csv(kappa.values$max.threshold,
          "Results/future/current/max_kappa_threshold.csv",
          row.names = FALSE)
