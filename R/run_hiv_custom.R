#' 

# Load packages ------------
if(!"here" %in% installed.packages()){
    install.packages(here)
}
source(here::here("R", "load_packages.R"))

## Read in data
hiv <- readRDS(here("data_raw", "hiv.rds"))

# Run custom brew and save RDS ------------
# Define sample of interest
sample_num <- which(hiv$group != "beads")[1]

# Define prior parameters
beads_ab0 <- beer:::getAB(subsetBeads(hiv))

a_pi <- 2
b_pi <- 300

a_phi <- 1.25
b_phi <- 0.1

a_c <- 80
b_c <- 20

# Adjust for super-enriched peptides
se_peps <- beer:::guessEnriched(hiv, method = "mle", beads.prior = beads_ab0)[, sample_num]
beads_ab0 <- getAB(subsetBeads(hiv)[!se_peps, ])

# Run JAGS
# Define data
data_list <- list(N = 1,
                  P = nrow(hiv) - sum(se_peps), 
                  B = 0, 
                  n = sum(counts(hiv)[!se_peps, sample_num]),
                  Y = counts(hiv)[!se_peps, sample_num, drop = FALSE], 
                  a_0 = beads_ab0$a_0, 
                  b_0 = beads_ab0$b_0, 
                  a_c = a_c, 
                  b_c = b_c, 
                  a_pi = a_pi, 
                  b_pi = b_pi, 
                  a_phi = a_phi, 
                  b_phi = b_phi, 
                  fc = 1)

# Define initial values
inits_list <- beer:::guessInits(hiv[, sample_num], beads.prior = beads_ab0)
inits_list$`.RNG.name` <- "base::Wichmann-Hill"
inits_list$`.RNG.seed` <- 20200416

model_file <- system.file("extdata/phipseq_model.bugs", package = "beer")

# Compile and run model
jags_model <- jags.model(file = model_file,
                         data = data_list,
                         inits = inits_list,
                         n.chains = 1)

jags_samples <- coda.samples(jags_model,
                             variable.names = c("c", "pi", "Z", "phi","theta"),
                             n.iter = 2e4)

saveRDS(jags_samples, here("data_processed", "hiv_samples.rds"))
