#' figure_simulation_postprob.R
#' 
#' Code to generate figures:
#' - simulation_postprob.png

# Set-up --------------
source(file.path("R", "load_packages.R"))
source(file.path("R", "helper_functions.R"))

sim_out <- readRDS(file.path("data_processed", "simulation_8beads_edgeR",
                             "sim_001.rds"))
sim_tidy <- as_df(sim_out, metadata = TRUE) %>%
    group_by(peptide) %>%
    mutate(sample = factor(sample, 1:ncol(sim_out)), 
           is_se = ifelse(sample != "beads" & is.na(beer_prob), TRUE, FALSE), 
           expected_prop = mean(counts[beads == "beads"]/n[beads == "beads"]), 
           expected_rc = expected_prop*n) %>%
    ungroup()

# Figure: simulation_postprob.png ------------
# Peptides from beads-only samples are excluded below. Super-enriched peptides 
# are colored in grey. 
sample_names <- c("beads_only", paste0("sample ", 1:10), "replicate of 10")
names(sample_names) <- 9:20

sim_tidy %>% 
    filter(beads != "beads") %>%
    arrange(beer_prob) %>%
    ggplot(aes(x = log10(expected_rc + 1), y = log10(counts + 1),
               color = beer_prob, group = sample)) +
    geom_point(size = 1) +
    geom_abline(aes(slope = 1, intercept = log10(true_c)), 
                size = 0.25, color = "black") +
    facet_wrap(sample ~., ncol = 4, 
               labeller = labeller(sample = sample_names)) +
    coord_fixed(xlim = c(0, 6), ylim = c(0, 6)) +
    labs(title = "Posterior probabilities of enrichment for simulated data",
         x = "log10(expected read counts + 1)",
         y = "log10(observed read counts + 1)",
         color = "posterior probability of enrichment") +
    scale_color_gradientn(colors = hot_cold_cols,
                          values = c(0, 0.15, 0.35, 0.45, 0.49, 0.51, 0.75, 1)) +
    theme_bw() +
    theme(legend.position = "bottom")

ggsave("figures/simulation_postprob.png", units = "in", width = 7.5)
