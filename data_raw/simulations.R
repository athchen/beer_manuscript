#' Function to generate simulated data sets based on HIV EC beads-only samples.
#' 
#' @param seed RNG seed
#' @return PhIPData object with the true parameters saved in assays/metadata 
simulate_data <- function(seed){
    set.seed(seed)
    
    ## Read-in HIV EC data
    hiv <- readRDS(here("data_raw", "hiv.rds"))
    
    ## Constants
    N <- 20
    P <- 1000
    B <- c(rep(1, 9), rep(0, 11))
    n <- round(abs(rnorm(N, 1e6, 1e5)))
    
    ## Beads-only parameters, sort peptides by proportion
    beads_params <- getAB(subsetBeads(hiv[1:P, ]), lower = 1) %>%
        as_tibble() %>%
        mutate(mean = a_0/(a_0 + b_0)) %>%
        arrange(mean)
    
    a_0 <- beads_params$a_0
    b_0 <- beads_params$b_0
    
    ## Proportion of peptides enriched
    Z <- cbind(matrix(0, nrow = P, ncol = 9), 
               sapply(1:10, function(x) sample(c(rep(1, 50), rep(0, P - 50)))))
    pi <- colMeans(Z)
    
    ## Fold-changes
    ## - 10 peptides with fc between 1-2
    ## - 10 peptides with fc between 2-4
    ## - 10 peptides with fc between 4-8
    ## - 10 peptides with fc between 8-10
    ## - 10 peptides with fc between 16-32
    ## Arrange the enriched peptides so that they are spaced across the range of
    ## proportion of reads pulled
    phi <- apply(Z, 2, function(x) {
        phi_j <- rep(1, P)
        phi_cat <- as.vector(matrix(c(runif(10, 1, 2), runif(10, 2, 4), 
                                      runif(10, 4, 8), runif(10, 8, 16), 
                                      runif(10, 16, 32)), 
                                    nrow = 5, byrow = TRUE))
        phi_j[which(x!= 0)] <- phi_cat
        phi_j
    })
    
    ## Duplicate last column for sample 20, add small noise for fc. 
    Z <- cbind(Z, Z[, N - 1])
    pi <- c(pi, pi[N - 1])
    phi <- cbind(phi, phi[, 19] + rnorm(P, 0, 0.01)*Z[, 19])
    
    ## Generate remaining parameters
    a <- matrix(NA, nrow = P, ncol = N)
    b <- matrix(NA, nrow = P, ncol = N)
    theta <- matrix(NA, nrow = P, ncol = N)
    Y <- matrix(NA, nrow = P, ncol = N)
    
    for(j in 1:N){
        for(i in 1:P) {
            if (Z[i, j] == 1 & B[j] == 0) {
                mean_e <- min(phi[i, j]*a_0[i]/(a_0[i] + b_0[i]), 1)
                var_e <- a_0[i]*b_0[i]/((a_0[i] + b_0[i])^2*
                                            (a_0[i] + b_0[i] + 1))
                
                a[i, j] <- max((1-mean_e)*mean_e^2/var_e - mean_e, 1)
                b[i, j] <- a[i, j]*(1/mean_e - 1)
                
            } else if (Z[i, j] == 0 & B[j] == 1) {
                
                a[i, j] <- a_0[i]
                b[i, j] <- b_0[i]
                
            } else {
                mean_ne <- a_0[i]/(a_0[i] + b_0[i])
                var_ne <- (a_0[i]*b_0[i])/((a_0[i]+b_0[i])^2*
                                               (a_0[i] + b_0[i] + 1))
                
                a[i, j] <- max((1-mean_ne)*mean_ne^2/var_ne - var_ne, 1)
                b[i, j] <- a[i, j]*(1/mean_ne - 1)
                
            }
            
            theta[i, j] <- rbeta(1, a[i, j], b[i, j])
            if(is.na(theta[i, j])){
                print(paste0("Sample ", j, " and peptide ", i, 
                             " failed to create theta."))
            }
            
            Y[i, j] <- rbinom(1, n[j], theta[i, j])
        }
    }
    
    ## Define new n
    n_init <- n
    n <- colSums(Y)
    
    ## Define estimate of c
    c <- sapply(1:N, function(x){
        if(x %in% 10:N) {
            ne_index <- which(Z[, x] == 0)
            data_sub <- data.frame(Y = Y[ne_index, x],
                                   expected_rc = a_0[ne_index]/
                                       (a_0[ne_index] + b_0[ne_index])*n[x])
            lm_fit <- lm(Y ~ expected_rc - 1, data = data_sub)
            coef(lm_fit)[[1]]
        } else 1
    })
    
    ## Create PhIPData object
    sample_params <- data.frame(group = c(rep("beads", sum(B) - 1), 
                                          rep("sample", N - sum(B) + 1)), 
                                n_init = n_init, n = n, 
                                true_c = c, true_pi = pi)
    pep_params <- data.frame(a_0 = a_0, b_0 = b_0)
    prior_params <- list(seed = seed,
                         a_pi = 1000, b_pi = 2.4e4, 
                         a_phi = 1.25, b_phi = 0.1, 
                         a_c = 33.02, b_c = 17.92, 
                         fc = 1)
    sim_data <- PhIPData(counts = Y, 
                         peptideInfo = pep_params, 
                         sampleInfo = sample_params, 
                         metadata = prior_params)
    assays(sim_data)[c("true_Z", "true_phi", "true_a", 
                       "true_b", "true_theta")] <-
        list(true_Z = Z, true_phi = phi, true_a = a, 
             true_b = b, true_theta = theta)
    
    sim_data
}
