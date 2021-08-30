#' Run edgeR and beer on HIV EC data. 

## Load packages
source(here("R", "load_packages.R"))

## Read in data
hiv <- readRDS(here("data_raw", "hiv.rds"))

## Run edgeR with beadsRR
hiv_out <- edgeR(hiv, assay.names = c("edgeR_logfc", "edgeR_prob"), 
                 parallel = "multisession", beadsRR = TRUE)

## Run beer with beadsRR
beer_assays <- c(phi = NULL, phi_Z = "beer_logfc", Z = "beer_prob", 
                 c = "sampleInfo", pi = "sampleInfo")
hiv_out <- brew(hiv_out, assay.names = beer_assays, beadsRR = TRUE, 
                parallel = "multisession")

saveRDS(hiv_out, here("data_processed", "hiv_results.rds"))