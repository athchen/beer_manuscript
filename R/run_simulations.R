#' Run edgeR and BEER for simulated data sets. 
#' This script was run on a cluster, which passed a number corresponding to 
#' which condition to run for the job. 
#' 
#' For example, to run the condition with edgeR starting parameters and only 2
#' beads-only samples, the command `Rscript run_simulations.R 1` would indicate
#' that setting. See the conditions data frame for all conditions. 

## Setup ----------
source(here("R", "load_packages.R"))
source(here("data_raw", "simulations.R"))

## Define conditions ----------
conditions <- data.frame(
    num_beads = rep(2^seq(3), each = 4),
    ab_method = rep(c("edgeR", "mle", "mom", "truth"), times = 3)
)

cmd_args <- commandArgs(trailingOnly = TRUE)
condition_num <- as.numeric(cmd_args)[1]

## Define simulations ----------
num_sims <- 10
set.seed(20210223)
sim_seeds <- c(20210223, sample(-1e6:1e6, num_sims - 1))

## Define simulation conditions
sim_names <- str_pad(seq(num_sims), 3, pad = "0")
sim_conditions <- data.frame(
    set = seq(num_sims),
    seed = sim_seeds,
    directory = paste0("sim_", sim_names)
)

## Set-up directories ----------
parent_directory <- file.path(
    here("data_processed"),
    paste0(
        "simulation_", 
        conditions$num_beads[condition_num], "beads_", 
        conditions$ab_method[condition_num]
    )
)
if(!dir.exists(parent_directory)){ dir.create(parent_directory, recursive = TRUE) }

## Run BEER and edgeR -----------
tmp <- lapply(seq(10), function(sim_num){
    print(sim_num)
    sim_data <- simulate_data(sim_conditions$seed[sim_num])
    
    ## Subset to the proper # of beads-only samples
    beads_index <- which(sim_data$group == "beads")
    sample_index <- which(sim_data$group != "beads")
    subset_index <- c(beads_index[seq(conditions$num_beads[condition_num])], 
                      sample_index)
    sim_data <- sim_data[, subset_index]
    
    ## UNCOMMENT THESE LINES 58, 59, AND 77 TO SAVE MCMC SAMPLES
    # tmpdir <- file.path(parent_directory, "tmp")
    # if(!dir.exists(tmpdir)){ dir.create(tmpdir, recursive = TRUE) }
    
    ## Run edgeR
    sim_out <- edgeR(sim_data, assay.names = c("edgeR_logfc", "edgeR_prob"))
    
    ## Run BEER
    beer_assays <- c(phi = NULL, phi_Z = "beer_logfc", Z = "beer_prob", 
                     c = "sampleInfo", pi = "sampleInfo")
    prior_params <- if(conditions$ab_method[condition_num] != "truth"){
        list(method = conditions$ab_method[condition_num])
    } else {
        list(method = "custom",
             a_0 = peptideInfo(sim_data)$a_0, b_0 = peptideInfo(sim_data)$b_0)
    }
    
    sim_out <- brew(
        sim_out, 
        # sample.dir = file.path(tmpdir, sim_conditions$directory[sim_num]),
        jags.params = list(seed = sim_conditions$seed[sim_num]), 
        prior.params = c(prior_params, 
                         list(a_pi = 2, b_pi = 300, a_phi = 1.25, b_phi = 0.1,
                              a_c = 80, b_c = 20, fc = 1)), 
        assay.names = beer_assays, 
        parallel = "multisession", 
    )
    
    ## Save output
    saveRDS(sim_out, file.path(
        parent_directory, 
        paste0(sim_conditions$directory[sim_num], ".rds")
    ))
})