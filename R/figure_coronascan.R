#' figure_coronascan.R
#' 
#' Code to generate figures:
#' - coronascan_protein.png
#' - coronascan_replicates.png
#' - coronascan_ranked_prob.png

# Set-up --------------
source(here("R", "load_packages.R"))
source(here("R", "helper_functions.R"))

# Read in data
cs <- readRDS(here("data_processed", "coronascan_results.rds"))

# Convert to tidy format and add hits for each approach
cs_tidy <- as(cs, "DataFrame") %>%
    as_tibble() %>%
    group_by(sample) %>%
    mutate(is_se = ifelse(sample != "beads" & is.na(beer_prob), TRUE, FALSE), 
           sample = factor(sample, levels = colnames(cs)),
           beer_hits = ifelse(beer_prob > 0.5 | is_se, 1, 0), 
           edgeR_bh = p.adjust(10^(-edgeR_prob), method = "BH"), 
           edgeR_hits = ifelse(edgeR_bh < 0.05, 1, 0)) %>%
    ungroup()

# Figure: coronascan_protein.png ------------
# Color palette
cs_species <- unique(cs_tidy$organism)
grey_palette <- palette(gray(seq(0.1, 0.8, len = (length(cs_species) - 1))))
num_sarscov2 <- grep("SARS-CoV2", cs_species)

# Facet labels
facet_labels <- paste0(colnames(cs), " (", 
                       ifelse(sampleInfo(cs)$group == "vrc", 
                              "VRC", sampleInfo(cs)$group), 
                       ")")
names(facet_labels) <- colnames(cs)

cs_tidy %>%
    filter(group != "beads") %>%
    select(sample, peptide, organism, protein_name, beer_hits, edgeR_hits) %>%
    group_by(sample, organism, protein_name) %>%
    summarize(prot_prop_Bayes = mean(beer_hits),
              prot_prop_edgeR = mean(edgeR_hits),
              num_peptides = n(), .groups = "drop") %>%
    arrange(desc(num_peptides)) %>%
    ggplot(aes(x = 0.5*(prot_prop_Bayes + prot_prop_edgeR),
               y = prot_prop_Bayes - prot_prop_edgeR,
               color = organism,
               size = num_peptides/2)) +
    facet_wrap(sample ~., ncol = 5, 
               labeller = labeller(sample = facet_labels)) +
    geom_hline(aes(yintercept = 0), size = 0.5, color = "grey50") +
    geom_vline(aes(xintercept = 0), size = 0.5, color = "grey50") +
    geom_point(alpha = 0.8) +
    labs(x = TeX("$\\frac{1}{2}$ (BEER + edgeR)"),
         y = "BEER - edgeR",
         color = "Organism",
         size = "# peptides") +
    scale_x_continuous(trans = "mysqrt", 
                       limits = c(0, 1), 
                       breaks = seq(0, 1, by = 0.2)) +
    scale_color_manual(values = c(grey_palette[1:(num_sarscov2 - 1)], 
                                  "firebrick2", 
                                  grey_palette[-(1:(num_sarscov2 - 1))]), 
                       breaks = cs_species) +
    theme_bw() +
    theme(aspect.ratio = 1, 
          legend.title = element_text(size = 10), 
          legend.text = element_text(size = 8), 
          legend.key.size = unit(0.75, "lines"), 
          legend.position = "bottom") +
    guides(color = guide_legend(ncol = 1, order = 1), 
           size = guide_legend(ncol = 1, order = 2))

ggsave("figures/coronascan_protein.png", units = "in", width = 10.5, height = 7)

# Figure: coronascan_replicates.png -----------
# Define sample to look at
cs_sample <- "CS 1"

# Plot for proportion of reads
prop_reads <- cs_tidy %>% 
    filter(sample == cs_sample) %>%
    mutate(prop_reads = counts/n) %>%
    select(pair_id, pair_num, prop_reads) %>%
    pivot_wider(names_from = "pair_num", 
                values_from = "prop_reads") %>%
    dplyr::rename(pair_1 = `1`, pair_2 = `2`) %>%
    ggplot(aes(x = log10(pair_1), 
               y = log10(pair_2))) + 
    geom_point(size = 1) +
    geom_abline(aes(intercept = 0, slope = 1)) + 
    labs(title = "Proportion of reads pulled", 
         x = "peptide pair #1, log10(proportion)", 
         y = "peptide pair #2, log10(proportion)") +
    scale_x_continuous(breaks = seq(-5, -2, by = 1),
                       minor_breaks = seq(-5, -2, by = 0.5), 
                       labels = TeX(paste0("$10^{", -5:-2, "}$"))) +
    scale_y_continuous(breaks = seq(-5, -2, by = 1),
                       minor_breaks = seq(-5, -2, by = 0.5), 
                       labels = TeX(paste0("$10^{", -5:-2, "}$"))) +
    theme_bw() +
    theme(aspect.ratio = 1, 
          title = element_text(size = 10))

# CAT curve for BEER posterior probabilities
#   For a fair comparison, set super-enriched peptides in the replicate sample
#   to have a rank of 1 as all of these peptides are identified as enriched
#   by BEER and edgeR. 
se_peps <- cs_tidy %>% 
    filter(sample == cs_sample & is_se) %>%
    group_by(pair_id) %>%
    filter(row_number() == 1) %>%
    ungroup() %>%
    mutate(rank = 1:n())

#   Ties are broken by the posterior marginal estimates of phi
beer_ranks <- cs_tidy %>%
    filter(sample == cs_sample) %>%
    group_by(pair_num) %>%
    arrange(desc(beer_prob), desc(beer_logfc), pair_id, .by_group = TRUE) %>%
    mutate(rank_se = nrow(se_peps) + cumsum(!pair_id %in% se_peps$pair_id)) %>%
    left_join(se_peps %>% select(pair_id, rank), by = c("pair_id")) %>%
    mutate(rank = ifelse(pair_id %in% se_peps$pair_id, rank, rank_se)) %>%
    select(pair_id, pair_num, rank) %>%
    pivot_wider(names_from = "pair_num", 
                values_from = "pair_id") %>%
    dplyr::rename(pair_1 = `1`, pair_2 = `2`) %>%
    arrange(rank)

beer_conc <- sapply(1:nrow(beer_ranks), function(rank){
    length(intersect(beer_ranks$pair_1[beer_ranks$rank <= rank], 
                     beer_ranks$pair_2[beer_ranks$rank <= rank]))/rank
})

# CAT curve for edgeR p-values
edgeR_ranks <- cs_tidy %>%
    filter(sample == cs_sample) %>%
    group_by(pair_num) %>%
    arrange(edgeR_bh, desc(edgeR_logfc), pair_id, .by_group = TRUE) %>%
    mutate(rank_se = nrow(se_peps) + cumsum(!pair_id %in% se_peps$pair_id)) %>%
    left_join(se_peps %>% select(pair_id, rank), by = c("pair_id")) %>%
    mutate(rank = ifelse(pair_id %in% se_peps$pair_id, rank, rank_se)) %>%
    select(pair_id, pair_num, rank) %>%
    pivot_wider(names_from = "pair_num", 
                values_from = "pair_id") %>%
    dplyr::rename(pair_1 = `1`, pair_2 = `2`) %>%
    arrange(rank)

edgeR_conc <- sapply(1:nrow(edgeR_ranks), function(rank){
    length(intersect(edgeR_ranks$pair_1[edgeR_ranks$rank <= rank], 
                     edgeR_ranks$pair_2[edgeR_ranks$rank <= rank]))/rank
})

cat_plot <- data.frame(rank = 1:length(beer_conc), 
                       BEER = beer_conc, 
                       edgeR = edgeR_conc) %>%
    pivot_longer(cols = c("BEER", "edgeR"), 
                 names_to = "approach", 
                 values_to = "concordance") %>%
    ggplot(aes(x = rank, y = concordance, color = approach)) +
    geom_line() +
    labs(title = "Concordance at the top", 
         y = "Concordance", 
         x = "Rank") +
    coord_cartesian(ylim = c(0, 1), xlim = c(0, 100)) +
    scale_x_continuous(breaks = seq(0, 100, by = 10)) +
    scale_y_continuous(breaks = seq(0, 1, by = 0.2)) +
    scale_color_manual(values = c("red", "black")) +
    theme_bw() + 
    theme(aspect.ratio = 1, 
          title = element_text(size = 10),
          legend.background = element_rect(color = "black", size = 0.3), 
          legend.position = c(0.8, 0.2))

ggarrange(prop_reads, cat_plot,
          nrow = 1, ncol = 2)

ggsave("figures/coronascan_replicates.png", units = "in", width = 8, height = 4)

# Figure: coronascan_ranked_prob.png -----------
post_prob <- cs_tidy %>%
    filter(sample == cs_sample) %>%
    select(pair_num, beer_prob) %>%
    group_by(pair_num) %>%
    arrange(desc(beer_prob), .by_group = TRUE) %>%
    mutate(rank = 1:n()) %>%
    ggplot(aes(x = rank, y = beer_prob, group = pair_num, color = pair_num)) + 
    geom_point(size = 1) +
    geom_line() +
    coord_cartesian(xlim = c(0, 100), ylim = c(0, 1)) +
    labs(title = "BEER", 
         x = "rank", 
         y = "posterior probabilities", 
         color = "peptide") +
    scale_color_manual(values = c("red", "black")) +
    theme_bw() + 
    theme(aspect.ratio = 1, 
          title = element_text(size = 10), 
          legend.background = element_rect(color = "black", size = 0.3), 
          legend.position = "none")

p_value <- cs_tidy %>%
    filter(sample == cs_sample) %>%
    select(pair_num, edgeR_prob) %>%
    group_by(pair_num) %>%
    arrange(desc(edgeR_prob), .by_group = TRUE) %>%
    mutate(rank = 1:n()) %>%
    ggplot(aes(x = rank, y = edgeR_prob, group = pair_num, color = pair_num)) + 
    geom_point(size = 1) +
    geom_line() +
    labs(title = "edgeR", 
         x = "rank", 
         y = "-log10(p-values)", 
         color = "peptide") +
    coord_cartesian(xlim = c(0, 100), ylim = c(0, 6)) +
    scale_y_continuous(oob = scales::oob_censor_any) +
    scale_color_manual(values = c("red", "black")) +
    theme_bw() +
    theme(aspect.ratio = 1, 
          title = element_text(size = 10), 
          legend.background = element_rect(color = "black", size = 0.3), 
          legend.position = c(0.825, 0.825))

ggarrange(post_prob, p_value, nrow = 1)

ggsave("figures/coronascan_ranked_prob.png", units = "in", width = 8, height = 4)

# Table: coronascan_replicates
cs_tidy %>%
    select(sample, pair_id, pair_num, beer_hits, edgeR_hits) %>%
    pivot_longer(cols = contains("hits"), 
                 names_to = "method", 
                 values_to = "hits", 
                 names_pattern = "([A-Za-z]*)_hits") %>%
    pivot_wider(names_from = pair_num, values_from = hits) %>%
    mutate(concordant = ifelse(`1` == `2`, 1, 0))%>%
    group_by(sample, method) %>%
    summarize(num = sum(concordant), 
              total_pep = n(), 
              prop = num/total_pep, .groups = "drop") %>%
    pivot_longer(cols = c("num", "prop"), 
                 names_to = "param", 
                 values_to = "value") %>%
    unite("method_param", method, param) %>%
    pivot_wider(names_from = "method_param", values_from = "value") %>%
    select(-total_pep) %>%
    mutate(across(contains(c("beer", "edgeR")), round, digits = 3))
