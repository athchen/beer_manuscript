#' figure_simulation_enrichments.R
#' 
#' Code to generate figures:
#' - simulation_enrichments_e.png
#' - simulation_enrichments_ne.png

# Set-up --------------
source(here("R", "load_packages.R"))
source(here("R", "helper_functions.R"))

# Define colors and facet names
sample_names <- c("beads_only", paste0("sample ", 1:10), "replicate of 10")
names(sample_names) <- 9:20

color_vector <- c("red4", "red", "blue", "grey50", "grey80")
names(color_vector) <- c("both", "both at 0.1", "BEER only", 
                         "not enriched in both", "not highlighted")

# Read in data
sim_out <- readRDS(file.path("data_processed", "simulation_8beads_edgeR",
                             "sim_001.rds"))

# BH correction and label predictions
sim_tidy <- as_df(sim_out, metadata = TRUE) %>%
    group_by(peptide) %>%
    mutate(sample = factor(sample, levels = 1:ncol(sim_out)), 
           is_se = ifelse(sample != "beads" & is.na(beer_prob), TRUE, FALSE), 
           expected_prop = mean(counts[beads == "beads"]/n[beads == "beads"]), 
           expected_rc = expected_prop*n) %>%
    ungroup() %>%
    group_by(sample) %>%
    mutate(bh_pvalues = p.adjust(10^(-edgeR_prob), method = "BH"), 
           edgeR_vs_Bayes = case_when((beer_prob > 0.5 | is_se) & bh_pvalues < 0.05 ~ "both", 
                                      (beer_prob > 0.5 | is_se) & bh_pvalues < 0.1 ~ "both at 0.1", 
                                      (beer_prob > 0.5 | is_se) ~ "BEER only", 
                                      beer_prob <= 0.5 & bh_pvalues < 0.05 ~ "edgeR only", 
                                      beer_prob <= 0.5 & bh_pvalues < 0.1 ~ "edgeR only, 0.1",
                                      TRUE ~ "not enriched in both"))
# Note: only certain BEER/edgeR enrichment combinations are present in this data
table(sim_tidy$edgeR_vs_Bayes)

# Figure: simulation_enrichments_e.png ------------
sim_tidy %>% 
    filter(beads != "beads") %>%
    mutate(point_color = ifelse(true_Z == 1, edgeR_vs_Bayes, 
                                "not highlighted")) %>%
    arrange(desc(point_color)) %>%
    ggplot(aes(x = log10(expected_rc + 1), y = log10(counts + 1), 
               color = point_color, group = sample)) +
    geom_point(size = 1, alpha = 0.75) +
    geom_abline(aes(slope = 1, intercept = log10(true_c)), 
                size = 0.25, color = "black") +
    coord_fixed(xlim = c(0, 6), ylim = c(0, 6)) +
    facet_wrap(sample ~., ncol = 4, 
               labeller = labeller(sample = sample_names)) +
    labs(title = "Enriched peptides",
         x = "log10(expected read count + 1)",
         y = "log10(observed read count + 1)",
         color = "predicted enrichment") +
    scale_color_manual(values = color_vector,
                       breaks = c("both", "both at 0.1", "BEER only",
                                  "not enriched in both"), 
                       labels = c("enriched in BEER and edgeR (FDR 0.05)", 
                                  "enriched in BEER and edgeR (FDR 0.10)", 
                                  "enriched in BEER only", 
                                  "not enriched in BEER and not enriched in edgeR")) +
    theme_bw() +
    guides(color = guide_legend(ncol = 1)) +
    theme(legend.position = "bottom")

ggsave("figures/simulation_enrichments_e.png", units = "in", 
       height = 6.75, width = 6)

# Figure: simulation_enrichments_ne.png ------------
sim_tidy %>% 
    filter(beads != "beads") %>%
    mutate(point_color = ifelse(true_Z == 0, edgeR_vs_Bayes,
                                "not highlighted")) %>%
    arrange(desc(point_color)) %>%
    ggplot(aes(x = log10(expected_rc + 1), y = log10(counts + 1), 
               color = point_color, group = sample)) +
    geom_point(size = 1) +
    geom_abline(aes(slope = 1, intercept = log10(true_c)), 
                size = 0.25, color = "black") +
    coord_fixed(xlim = c(0, 6), ylim = c(0, 6)) +
    facet_wrap(sample ~., ncol = 4, 
               labeller = labeller(sample = sample_names)) +
    labs(title = "Non-enriched peptides",
         x = "log10(expected read count + 1)",
         y = "log10(observed read count + 1)",
         color = "predicted enrichment") +
    scale_color_manual(values = color_vector,
                       breaks = c("both", "both at 0.1", "BEER only",
                                  "not enriched in both"), 
                       labels = c("enriched in BEER and edgeR (FDR 0.05)", 
                                  "enriched in BEER and edgeR (FDR 0.10)", 
                                  "enriched in BEER only", 
                                  "not enriched in BEER and not enriched in edgeR")) +
    theme_bw() +
    guides(color = guide_legend(ncol = 1)) +
    theme(legend.position = "bottom")

ggsave("figures/simulation_enrichments_ne.png", units = "in",
       height = 6.75, width = 6)
