#' figure_simulation_fc.R
#' 
#' Code to generate figures:
#' - simulation_fc.png

# Set-up --------------
source(here("R", "load_packages.R"))
source(here("R", "helper_functions.R"))

sim_out <- readRDS(here("data_processed", "simulation_8beads_edgeR",
                        "sim_001.rds"))

sim_tidy <- as_df(sim_out, metadata = TRUE) %>%
    group_by(peptide) %>%
    mutate(is_se = ifelse(beads != "beads" & is.na(beer_prob), TRUE, FALSE), 
           expected_prop = mean(counts[beads == "beads"]/n[beads == "beads"]), 
           expected_rc = expected_prop*n) %>%
    ungroup()

# Calculate true_c
true_c <- sim_tidy %>%
    group_by(sample) %>%
    group_split() %>%
    map_dfr(function(df){
        c_coef <- if(unique(df$beads) != "beads" &
                     as.numeric(unique(df$sample)) != 9){
            lm_fit <- lm(counts ~ expected_rc, data = df %>% filter(true_Z == 0))
            coef(lm_fit)[["expected_rc"]]
        } else 1
        data.frame(sample = unique(df$sample), 
                   true_c = c_coef)
    }) %>%
    arrange(as.numeric(sample))

# Figure: simulation_fc.png ------------
# Peptides from beads-only samples and super-enriched peptides are excluded
# from the following plots
beer_phi <- sim_tidy %>%
    filter(beads == "sample" & !is.na(beer_fc)) %>%
    mutate(cond_fc = ifelse(beer_prob < 0.5, 1, beer_fcZ)) %>%
    arrange(beer_prob) %>%
    ggplot(aes(x = log2(true_phi), y = log2(cond_fc), color = beer_prob)) +
    geom_point(alpha = 0.75) +
    labs(title = "BEER with edgeR", 
         x = TeX("$\\log_2$ actual fold change"),
         y = TeX("$\\log_2$ predicted fold change | $\\hat{Z}_{ij}$ "),
         color = "posterior probability") +
    geom_abline(aes(intercept = 0, slope = 1), color = "black") +
    geom_hline(aes(yintercept = log2(1)), linetype = "dashed") +
    geom_vline(aes(xintercept = log2(1)), linetype = "dashed") +
    scale_y_continuous(limits = c(0, 4.5), oob = scales::squish) +
    scale_x_continuous(limits = c(0, 4.5), oob = scales::squish) +
    scale_color_gradientn(colors = hot_cold_cols,
                          values = c(0, 0.15, 0.35, 0.45, 0.49, 0.51, 0.75, 1)) +
    coord_fixed(ratio = 1) +
    guides(color = guide_colorbar(barwidth = unit(1.75, "in"))) +
    theme_bw() +
    theme(legend.position = "none")

beer_cphi <- sim_tidy %>%
    left_join(true_c, by = c("sample")) %>%
    filter(beads == "sample" & !is.na(beer_fc)) %>%
    mutate(cond_fc = ifelse(beer_prob < 0.5, 1, beer_fcZ)) %>%
    arrange(beer_prob) %>%
    ggplot(aes(x = log2(true_phi*true_c), y = log2(cond_fc*est_c), 
               color = beer_prob)) +
    geom_point(alpha = 0.75) +
    labs(title = "BEER with edgeR", 
         x = TeX("$c_{j}\\cdot\\log_2$ actual fold change"),
         y = TeX("$\\hat{c}_{j}\\cdot\\log_2$ predicted fold change | $\\hat{Z}_{ij}$ "),
         color = "posterior probability") +
    geom_abline(aes(intercept = 0, slope = 1), color = "black") +
    geom_hline(aes(yintercept = log2(1)), linetype = "dashed") +
    geom_vline(aes(xintercept = log2(1)), linetype = "dashed") +
    scale_y_continuous(limits = c(0, 4.5), oob = scales::squish) +
    scale_x_continuous(limits = c(0, 4.5), oob = scales::squish) +
    scale_color_gradientn(colors = hot_cold_cols,
                          values = c(0, 0.15, 0.35, 0.45, 0.49, 0.51, 0.75, 1)) +
    coord_fixed(ratio = 1) +
    guides(color = guide_colorbar(barwidth = unit(1.75, "in"))) +
    theme_bw() +
    theme(legend.position = "bottom", 
          legend.text = element_text(size = 8, angle = 45, hjust = 1), 
          legend.title = element_text(vjust = 0.85))

edgeR_phi <- sim_tidy %>%
    filter(beads == "sample" & !is.na(beer_fc)) %>%
    group_by(sample) %>%
    mutate(edgeR_prob_bh = p.adjust(10^(-edgeR_prob), method = "BH"), 
           edgeR_logprob_bh = -log10(edgeR_prob_bh)) %>%
    arrange(edgeR_logprob_bh) %>%
    ggplot(aes(x = log2(true_phi), y = edgeR_logfc, color = edgeR_logprob_bh)) +
    geom_point(alpha = 0.75) +
    labs(title = "edgeR", 
         x = TeX("$\\log_2$ actual fold change"),
         y = TeX("$\\log_2$ predicted fold change"),
         color = "BH adjusted p-values") +
    geom_abline(aes(intercept = 0, slope = 1), color = "black") +
    geom_hline(aes(yintercept = log2(1)), linetype = "dashed") +
    geom_vline(aes(xintercept = log2(1)), linetype = "dashed") +
    scale_y_continuous(limits = c(0, 4.5), oob = scales::squish) +
    scale_x_continuous(limits = c(0, 4.5), oob = scales::squish) +
    scale_color_gradientn(colors = hot_cold_cols,
                          limits = c(0, 3), 
                          breaks = c(0, 1, -log10(0.05), 2, 3),
                          labels = c(1, 0.1, 0.05, 0.01, 0.001), 
                          values = scales::rescale(c(0, 0.5, 0.75, 1, 1.25, 1.5, 2, 3)), 
                          oob = scales::squish) +
    coord_fixed(ratio = 1) +
    guides(color = guide_colorbar(barwidth = unit(1.75, "in"))) +
    theme_bw() +
    theme(legend.position = "none")

edgeR_cphi <- sim_tidy %>%
    left_join(true_c, by = c("sample")) %>%
    filter(beads == "sample" & !is.na(beer_fc)) %>%
    group_by(sample) %>%
    mutate(edgeR_prob_bh = p.adjust(10^(-edgeR_prob), method = "BH"), 
           edgeR_logprob_bh = -log10(edgeR_prob_bh)) %>%
    arrange(edgeR_logprob_bh) %>%
    ggplot(aes(x = log2(true_phi*true_c), y = edgeR_logfc, color = edgeR_logprob_bh)) +
    geom_point(alpha = 0.75) +
    labs(title = "edgeR", 
         x = TeX("$c_{j}\\cdot\\log_2$ actual fold change"),
         y = TeX("$\\log_2$ predicted fold change"),
         color = "BH adjusted p-values") +
    geom_abline(aes(intercept = 0, slope = 1), color = "black") +
    geom_hline(aes(yintercept = log2(1)), linetype = "dashed") +
    geom_vline(aes(xintercept = log2(1)), linetype = "dashed") +
    scale_y_continuous(limits = c(0, 4.5), oob = scales::squish) +
    scale_x_continuous(limits = c(0, 4.5), oob = scales::squish) +
    scale_color_gradientn(colors = hot_cold_cols,
                          limits = c(0, 3), 
                          breaks = c(0, 1, -log10(0.05), 2, 3),
                          labels = c(1, 0.1, 0.05, 0.01, 0.001), 
                          values = scales::rescale(c(0, 0.5, 0.75, 1, 1.25, 1.5, 2, 3)), 
                          oob = scales::squish) +
    coord_fixed(ratio = 1) +
    guides(color = guide_colorbar(barwidth = unit(1.75, "in"))) +
    theme_bw() +
    theme(legend.position = "bottom", 
          legend.text = element_text(size = 8, angle = 45, hjust = 1), 
          legend.title = element_text(vjust = 0.85))

simulation_fc <- ggarrange(beer_phi, edgeR_phi, beer_cphi, edgeR_cphi, 
                           nrow = 2, ncol = 2, align = "v", 
                           heights = c(4, 4.9), labels = "AUTO")

ggsave("figures/simulation_fc.png", simulation_fc, 
       units = "in", width = 8, height = 9)
