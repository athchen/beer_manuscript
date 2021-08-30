#' Run edgeR and beer on CoronaScan data. 

## Load packages
source(here("R", "load_packages.R"))

## Read in data
cs <- readRDS(here("data_raw", "coronascan.rds"))

## Run edgeR with beadsRR
cs_out <- edgeR(cs, assay.names = c("edgeR_logfc", "edgeR_prob"), 
                parallel = "multisession", beadsRR = TRUE)

## Run beer with beadsRR
beer_assays <- c(phi = NULL, phi_Z = "beer_logfc", Z = "beer_prob", 
                 c = "sampleInfo", pi = "sampleInfo")
cs_out <- brew(cs_out, assay.names = beer_assays, beadsRR = TRUE, 
                parallel = "multisession")

saveRDS(cs_out, here("data_processed", "coronascan_results.rds"))