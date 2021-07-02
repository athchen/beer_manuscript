source(file.path("R", "load_packages.R"))
sim_dirs <- list.files("data_processed", "simulation_[0-9]", full.names = TRUE)
sim_list <- sapply(sim_dirs, list.files, full.names = TRUE) %>% as.vector()

# Read and tidy data ---------
sim_data <- lapply(sim_list, function(sim_num){
    
    # Extract setting (data set, # beads used, method of estimating a_0, b_0)
    data_num <- str_match(sim_num, "sim_([0-9]*)\\.rds")[, 2] %>% as.numeric()
    num_beads <- str_match(sim_num, "([0-9])beads")[, 2] %>% as.numeric()
    ab_method <- str_match(sim_num, "[0-9]beads_(.*)/")[, 2]

    # Read in data
    all_results <- readRDS(sim_num)

    N <- ncol(all_results)
    P <- nrow(all_results)
    B <- sum(all_results$group == getBeadsName())
    n <- librarySize(all_results)

    a_pi <- metadata(all_results)$a_pi
    b_pi <- metadata(all_results)$b_pi
    a_c <- metadata(all_results)$a_c
    b_c <- metadata(all_results)$b_c
    a_phi <- metadata(all_results)$a_phi
    b_phi <- metadata(all_results)$b_phi

    a_0 <- peptideInfo(all_results)$a_0
    b_0 <- peptideInfo(all_results)$b_0

    pi <- all_results$pi
    c <- all_results$c

    Z <- assay(all_results, "true_Z")
    phi <- assay(all_results, "true_phi")
    a <- assay(all_results, "true_a")
    b <- assay(all_results, "true_b")
    theta <- assay(all_results, "true_theta")
    Y <- counts(all_results)

    # Arrange in a df
    data.frame(sim_num = data_num,
               num_beads = num_beads, 
               ab_method = ab_method, 
               sample = rep(1:N, each = P),
               peptide = rep(1:P, times = N),
               n = rep(n, each = P),
               c = rep(c, each = P),
               a_0 = a_0,
               b_0 = b_0,
               phi = as.vector(phi),
               Z = as.vector(Z),
               Y = as.vector(Y),
               post_prob = as.vector(assay(all_results, "beer_prob")),
               est_phi = as.vector(assay(all_results, "beer_fc")),
               est_phiZ = as.vector(assay(all_results, "beer_fcZ")),
               edgeR_fc = as.vector(assay(all_results, "edgeR_logfc")),
               edgeR_pval = as.vector(assay(all_results, "edgeR_prob"))) %>%
        mutate(is_se = ifelse(is.na(post_prob & sample %in% (B+1):N), 1, 0),
               post_prob = ifelse(is_se, 1, post_prob),
               pred_enriched = ifelse(post_prob > 0.5, 1, 0)) %>%
        group_by(peptide) %>%
        mutate(expected_prop = mean(Y[sample %in% 1:B]/n[sample %in% 1:B]),
               expected_rc = expected_prop*n) %>%
        ungroup()
})

names(sim_data) <- sim_list

# Clean-up variable space
rm(list = ls()[ls() != "sim_data"])
