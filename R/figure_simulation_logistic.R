#' figure_simulation_logistic.R
#' 
#' Code to generate figures:
#' - simulation_logistic_fc.png
#' - simulation_logistic_rc.png

# Set-up --------------
# Get posterior probability cutoff with FDR rate of around 5\%
# Comment out lines 21 & 22 to see variation in cutoffs by simulation
source(file.path("R", "load_curves.R"))
beer_fdr <- beer_roc_by_fc %>%
    filter(approach == "BEER" & ab_method == "edgeR" &
               group == "full data" & num_beads == 8) %>%
    group_by(approach, ab_method, sim_num, num_beads) %>%
    mutate(dist_fdr = abs(ppv - 0.95)) %>%
    filter(dist_fdr == min(dist_fdr, na.rm = TRUE)) %>%
    summarize(mean_cutoff = mean(cutoff),
              median_cutoff = median(cutoff), 
              range_cutoff = max(cutoff) - min(cutoff), 
              .groups = "drop") %>%
    filter(sim_num == 1) %>%
    pull(mean_cutoff)

rm(list = setdiff(ls(), c("beer_fdr", lsf.str())))

# Read in simulation data, make sure sim # corresponds to the beer_fdr
sim_out <- readRDS(file.path("data_processed", "simulation_8beads_edgeR",
                             "sim_001.rds"))

sim_tidy <- as_df(sim_out, metadata = TRUE) %>%
    group_by(peptide) %>%
    mutate(is_se = ifelse(beads != "beads" & is.na(beer_prob), TRUE, FALSE), 
           expected_prop = mean(counts[beads == "beads"]/n[beads == "beads"]), 
           expected_rc = expected_prop*n, 
           beer_hits = ifelse(beer_prob > beer_fdr, 1, 0)) %>%
    ungroup() %>%
    group_by(sample) %>%
    mutate(edgeR_bh = p.adjust(10^(-edgeR_prob), method = "BH")) %>%
    ungroup()

# Figure: simulation_logistic_fc.png ----------
# BEER
beer_fit <- glm(beer_hits ~ true_phi, family = binomial(link = "logit"),
                data = sim_tidy %>% filter(true_Z == 1))

beer_pred <- penriched_fit(beer_fit, 
                           data.frame(true_phi = seq(1, 100, length = 1000))) %>%
    mutate(approach = "BEER", 
           threshold = "FDR - 0.05")

# edgeR - Bonferroni
edgeR_bon_fit <- glm(edgeR_hits ~ true_phi, family = binomial(link = "logit"),
                    data = sim_tidy %>% 
                        mutate(edgeR_hits = ifelse(edgeR_prob > -log10(5e-5), 1, 0)) %>%
                        filter(true_Z == 1))

edgeR_bon_pred <- penriched_fit(edgeR_bon_fit, 
                                data.frame(true_phi = seq(1, 100, length = 1000))) %>%
    mutate(approach = "edgeR", 
           threshold = "Bonferroni")

# edgeR - BH 0.05
edgeR_05_fit <- glm(edgeR_hits ~ true_phi, family = binomial(link = "logit"),
                    data = sim_tidy %>% 
                        mutate(edgeR_hits = ifelse(edgeR_bh < 0.05, 1, 0)) %>%
                        filter(true_Z == 1))

edgeR_05_pred <- penriched_fit(edgeR_05_fit, 
                               data.frame(true_phi = seq(1, 100, length = 1000))) %>%
    mutate(approach = "edgeR", 
           threshold = "BH - 0.05")

logistic_fc <- bind_rows(beer_pred, edgeR_bon_pred, edgeR_05_pred) %>%
    ggplot(aes(x = log2(true_phi), y = predict_p, group = paste0(approach, threshold))) +
    geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci, fill = approach), alpha = 0.2) +
    geom_line(aes(color = approach, linetype = threshold)) +
    labs(x = "true fold change",
         y = "probability of being classified as enriched",
         color = "approach") +
    scale_color_manual(breaks = c("BEER", "edgeR"),
                       values = c("red", "black"),
                       labels = c("BEER", "edgeR")) +
    scale_fill_manual(breaks = c("BEER", "edgeR"),
                      values = c("red", "black"),
                      labels = c("BEER", "edgeR")) +
    scale_y_continuous(breaks = seq(0, 1, 0.1),
                       minor_breaks = seq(0, 1, 0.05)) +
    scale_x_continuous(breaks = 0:6,
                       labels = 2^(0:6)) +
    scale_linetype_manual(values = c("dashed", "dotted", "solid"), 
                          breaks = c("FDR - 0.05", "BH - 0.05", "Bonferroni")) +
    theme_bw() +
    theme(aspect.ratio = 1)

ggsave("figures/simulation_logistic_fc.png", logistic_fc,
       units = "in", width = 6, height = 4.5)

# Figure: simulation_logistic_rc.png ----------
# BEER
beer_fit <- glm(beer_hits ~ rpm, family = binomial(link = "logit"),
                data = sim_tidy %>% filter(true_Z == 1) %>%
                    mutate(rpm = counts/n*1e6))

beer_pred <- penriched_fit(beer_fit, 
                           data.frame(rpm = seq(1, 2.5*1e5, length = 1000))) %>%
    mutate(approach = "BEER", 
           threshold = "FDR - 0.05")

# edgeR - Bonferroni
edgeR_bon_fit <- glm(edgeR_hits ~ rpm, family = binomial(link = "logit"),
                     data = sim_tidy %>% 
                         mutate(edgeR_hits = ifelse(edgeR_prob > -log10(5e-5), 1, 0), 
                                rpm = counts/n*1e6) %>%
                         filter(true_Z == 1))

edgeR_bon_pred <- penriched_fit(edgeR_bon_fit, 
                                data.frame(rpm = seq(1, 2.5*1e5, length = 1000))) %>%
    mutate(approach = "edgeR", 
           threshold = "Bonferroni")

# edgeR - BH 0.05
edgeR_05_fit <- glm(edgeR_hits ~ rpm, family = binomial(link = "logit"),
                    data = sim_tidy %>% 
                        mutate(edgeR_hits = ifelse(edgeR_bh < 0.05, 1, 0), 
                               rpm = counts/n*1e6) %>%
                        filter(true_Z == 1))

edgeR_05_pred <- penriched_fit(edgeR_05_fit, 
                               data.frame(rpm = seq(1, 2.5*1e5, length = 1000))) %>%
    mutate(approach = "edgeR", 
           threshold = "BH - 0.05")

logistic_rc <- bind_rows(beer_pred, edgeR_bon_pred, edgeR_05_pred) %>%
    ggplot(aes(x = log10(rpm), y = predict_p, group = paste0(approach, threshold))) +
    geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci, fill = approach), alpha = 0.2) +
    geom_line(aes(color = approach, linetype = threshold)) +
    labs(x = "reads per million",
         y = "probability of being classified as enriched",
         color = "approach") +
    scale_color_manual(breaks = c("BEER", "edgeR"),
                       values = c("red", "black"),
                       labels = c("BEER", "edgeR")) +
    scale_fill_manual(breaks = c("BEER", "edgeR"),
                      values = c("red", "black"),
                      labels = c("BEER", "edgeR")) +
    scale_y_continuous(breaks = seq(0, 1, 0.1),
                       minor_breaks = seq(0, 1, 0.05)) +
    scale_x_continuous(breaks = 0:5,
                       minor_breaks = seq(0, 5, 0.5),
                       labels = c(1, 10, 100, sapply(paste0("$10^{", 3:5, "}$"), TeX))) +
    scale_linetype_manual(values = c("dashed", "dotted", "solid"), 
                          breaks = c("FDR - 0.05", "BH - 0.05", "Bonferroni")) +
    theme_bw() +
    theme(aspect.ratio = 1)

ggsave("figures/simulation_logistic_rc.png", logistic_rc,
       units = "in", width = 6, height = 4.5)
