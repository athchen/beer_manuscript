#' load_curves.R
#' 
#' Code to process and load sens/spec/ppv data.

# Set-up ----------
source(file.path("R", "load_simulations.R"))
source(file.path("R", "helper_functions.R"))

# Create ROC/PRC data for each data set and simulation condition ------------
# Note: the following code takes a while to run. 
fc_thresholds <- c(1, 2, 4, 8, 16, 32, NA)

beer_roc_by_fc <- lapply(sim_data, function(sim){
    
    roc_by_fc <- map_dfr(2:length(fc_thresholds), function(pos){
        
        fc_min <- ifelse(is.na(fc_thresholds[pos]), 1, fc_thresholds[pos - 1])
        fc_max <- ifelse(is.na(fc_thresholds[pos]), 32, fc_thresholds[pos])
        
        n_beads <- unique(sim$num_beads)
        n_samples <- max(sim$sample)
        
        sim %>%
            filter(((phi <= fc_max & phi > fc_min) | Z == 0) &
                       (sample %in% seq(n_beads + 1, n_samples))) %>%
            select(post_prob, Z) %>%
            dplyr::rename(prop_enriched = post_prob) %>%
            get_roc(., min_cutoff = 0, max_cutoff = 1 + 1e-6) %>%
            mutate(approach = "BEER",
                   sim_num = unique(sim$sim_num), 
                   num_beads = n_beads,
                   ab_method = unique(sim$ab_method), 
                   group = ifelse(is.na(fc_thresholds[pos]), "full data",
                                  paste0(fc_min, "<phi<=", fc_max)))
    }) %>%
        mutate(group_lab = factor(group, 
                                  levels = c("full data", "1<phi<=2", "2<phi<=4",
                                             "4<phi<=8", "8<phi<=16", "16<phi<=32"), 
                                  labels = c(TeX("full data"), 
                                             TeX("$1<\\phi_{ij}\\leq 2$"), 
                                             TeX("$2<\\phi_{ij}\\leq 4$"),
                                             TeX("$4<\\phi_{ij}\\leq 8$"),
                                             TeX("$8<\\phi_{ij}\\leq 16$"),
                                             TeX("$16<\\phi_{ij}\\leq 32$"))))
}) %>% plyr::ldply(.id = NULL)

# Simulation output for edgeR data is run for each method of estimating a_0, b_0,
# so to avoid duplication, we just use the data from the "truth" runs with
# 4 and 8 beads-only samples. 
edgeR_roc_by_fc <- lapply(sim_data[c(1:10, 41:50)], function(sim){
    
    roc_by_fc <- map_dfr(2:length(fc_thresholds), function(pos){
        
        fc_min <- ifelse(is.na(fc_thresholds[pos]), 1, fc_thresholds[pos - 1])
        fc_max <- ifelse(is.na(fc_thresholds[pos]), 32, fc_thresholds[pos])
        
        n_beads <- unique(sim$num_beads)
        n_samples <- max(sim$sample)
        
        sim %>%
            filter(((phi <= fc_max & phi > fc_min) | Z == 0) &
                       (sample %in% seq(n_beads + 1, n_samples))) %>%
            select(edgeR_pval, Z) %>%
            dplyr::rename(prop_enriched = edgeR_pval) %>%
            get_roc(., min_cutoff = 0, max_cutoff = max(.$prop_enriched)) %>%
            mutate(approach = "edgeR",
                   sim_num = unique(sim$sim_num), 
                   num_beads = n_beads,
                   ab_method = unique(sim$ab_method), 
                   group = ifelse(is.na(fc_thresholds[pos]), "full data",
                                  paste0(fc_min, "<phi<=", fc_max)))
    }) %>%
        mutate(group_lab = factor(group, 
                                  levels = c("full data", "1<phi<=2", "2<phi<=4",
                                             "4<phi<=8", "8<phi<=16", "16<phi<=32"), 
                                  labels = c(TeX("full data"), 
                                             TeX("$1<\\phi_{ij}\\leq 2$"), 
                                             TeX("$2<\\phi_{ij}\\leq 4$"),
                                             TeX("$4<\\phi_{ij}\\leq 8$"),
                                             TeX("$8<\\phi_{ij}\\leq 16$"),
                                             TeX("$16<\\phi_{ij}\\leq 32$"))))
}) %>% plyr::ldply(.id = NULL)

roc_by_fc <- bind_rows(beer_roc_by_fc, edgeR_roc_by_fc) %>%
    mutate(method = ifelse(approach == "edgeR", "edgeR", paste0(approach, "_", ab_method)),
           method = factor(method,
                           levels = c("BEER_truth", "BEER_mom", "BEER_mle",
                                      "BEER_edgeR", "edgeR")))

# Average individual curves -----------
# Remote fc > 16 since the curves are essentially perfect which causes problems 
# for interpolation.
interpolate <- roc_by_fc %>%
    filter(group!= "16<phi<=32") %>%
    group_by(sim_num, group, group_lab, method, num_beads) %>%
    group_split() %>%
    map_dfr(function(sub_df){
        
        specm_range <- seq(0, 1, by = 0.01)
        sens_approx <-  approx(c(0, 1 - sub_df$spec, 1), 
                               c(0, sub_df$sens, 1), specm_range, 
                               yleft = 0, yright = 1, rule = 1)$y
        ppv_approx <-  approx(c(0, sub_df$sens), 
                              c(1, sub_df$ppv), specm_range, 
                              yleft = 1, yright = 0, rule = 1)$y
        
        data.frame(sim_num = unique(sub_df$sim_num),
                   group = unique(sub_df$group), 
                   group_lab = unique(sub_df$group_lab), 
                   method = unique(sub_df$method),
                   num_beads = unique(sub_df$num_beads),
                   x = specm_range, 
                   sens = sens_approx, 
                   ppv = ppv_approx)
    })
