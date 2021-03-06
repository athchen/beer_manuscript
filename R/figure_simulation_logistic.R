#' figure_simulation_logistic.R
#' 
#' Code to generate figures:
#' - Figure 2: simulation_logistic_fc.png

# Set-up --------------
if(!"here" %in% installed.packages()){
    install.packages(here)
}

if(file.exists(here::here("data_processed", "simulation_curves.rda"))){
    load(here::here("data_processed", "simulation_curves.rda"))
    source(here::here("R", "load_packages.R"))
} else {
    # Takes a while to run - only need to run once. 
    source(here::here("R", "load_curves.R"))
}

# Get posterior probability cutoff with FDR rate of around 5\%
# Comment out lines 21 & 22 to see variation in cutoffs by simulation
beer_fdr <- beer_roc_by_fc %>%
    filter(approach == "BEER" & ab_method == "edgeR" &
               group == "full data" & num_beads == 8) %>%
    group_by(approach, ab_method, sim_num, num_beads) %>%
    # mutate(dist_fdr = abs(ppv - 0.95)) %>%
    # filter(dist_fdr == min(dist_fdr, na.rm = TRUE)) %>%
    summarize(mean_cutoff = mean(cutoff),
              median_cutoff = median(cutoff), 
              range_cutoff = max(cutoff) - min(cutoff), 
              .groups = "drop") 

rm(list = setdiff(ls(), c("beer_fdr", "sim_data", lsf.str())))

# Tidy simulation data, add edgeR and bh hits
sim_tidy <- lapply(sim_data[grepl("8beads_edgeR", names(sim_data))], function(sim){
    sim %>% 
        left_join(beer_fdr[, c("sim_num", "mean_cutoff")], by = "sim_num") %>%
        mutate(beer_hits = ifelse(post_prob > mean_cutoff, 1, 0), 
               rpm = Y/n*1e6) %>%
        group_by(sample) %>%
        mutate(edgeR_bh = p.adjust(10^(-edgeR_pval), method = "BH")) %>%
        ungroup()
})

# Figure 2: simulation_logistic_fc.png ----------
# BEER
beer_pred <- lapply(sim_tidy, function(sim){
    fit <- glm(beer_hits ~ phi, family = binomial(link = "logit"),
               data = sim %>% filter(Z == 1))
    
    penriched_fit(fit, data.frame(phi = seq(1, 100, length = 1000))) %>%
        mutate(sim_num = unique(sim$sim_num), 
               approach = "BEER", 
               threshold = "FDR - 0.05")
}) %>% plyr::ldply(.id = NULL)

# edgeR - BH 0.05
edgeR_05_pred <- lapply(sim_tidy, function(sim){
    fit <- glm(edgeR_hits ~ phi, family = binomial(link = "logit"),
               data = sim %>% 
                   mutate(edgeR_hits = ifelse(edgeR_bh < 0.05, 1, 0)) %>%
                   filter(Z == 1))
    
    penriched_fit(fit, data.frame(phi = seq(1, 100, length = 1000))) %>%
        mutate(sim_num = unique(sim$sim_num), 
               approach = "edgeR", 
               threshold = "BH - 0.05")
}) %>% plyr::ldply(.id = NULL)

logistic_fc <- bind_rows(beer_pred, edgeR_05_pred) %>%
    # group_by(approach, threshold, phi) %>%
    # summarize(mean_prob = mean(predict_p), 
    #           max_prob = max(predict_p), 
    #           min_prob = min(predict_p), 
    #           .groups = "drop") %>%
    ggplot(aes(x = log2(phi), y = predict_p)) +
    geom_line(aes(color = approach, 
                  group = paste0(approach, threshold, sim_num)), 
              size = 0.1, alpha = 0.5) +
    geom_line(aes(x = log2(phi), y = mean_prob, color = approach), 
              size = 1,
              data = bind_rows(beer_pred, edgeR_05_pred) %>%
                  group_by(approach, threshold, phi) %>%
                  summarize(mean_prob = mean(predict_p), .groups = "drop")) +
    labs(x = "true fold change",
         y = "probability of being classified as enriched",
         color = "method") +
    scale_color_manual(breaks = c("BEER", "edgeR"),
                       values = c("red", "black"),
                       labels = c("BEER", "edgeR")) +
    scale_y_continuous(breaks = seq(0, 1, 0.1),
                       minor_breaks = seq(0, 1, 0.05)) +
    scale_x_continuous(breaks = 0:4,
                       labels = 2^(0:4), 
                       limits = c(0, 4)) +
    theme_bw() +
    theme(aspect.ratio = 1)

ggsave("figures/simulation_logistic_fc.png", logistic_fc, dpi = 600, 
       units = "in", width = 6, height = 4.5)
