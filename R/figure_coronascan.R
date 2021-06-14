#' figure_coronascan.R
#' 
#' Code to generate figures:
#' - coronascan_protein.png
#' - coronascan_replicates.png

# Set-up --------------
source(file.path("R", "load_packages.R"))
source(file.path("R", "helper_functions.R"))

# Read in data
cs <- readRDS(file.path("data_processed", "coronascan_132.rds"))

sample_names <- paste0("CS ", 1:ncol(cs))
names(sample_names) <- colnames(cs)

# Add peptide information, add hits for each approach
cs_tidy <- as_df(cs) %>%
    left_join(as_tibble(mcols(peptideInfo(cs))) %>%
                  mutate(organism = gsub("\\(Severe .*)","", organism)), 
              by = c("peptide" = "pep_num")) %>%
    group_by(sample) %>%
    mutate(is_se = ifelse(sample != "beads" & is.na(beer_prob), TRUE, FALSE), 
           sample_id = sample_names[sample], 
           sample_id = factor(sample_id, levels = sample_names),
           beer_hits = ifelse(beer_prob > 0.5 | is_se, 1, 0), 
           edgeR_bh = p.adjust(10^(-edgeR_prob), method = "BH"), 
           edgeR_hits = ifelse(edgeR_bh < 0.05, 1, 0)) %>%
    ungroup()

# Figure: coronascan_protein.png ------------
cs_tidy %>%
    filter(beads != "beads") %>%
    select(sample, sample_id, peptide, organism, protein_name, beer_hits, edgeR_hits) %>%
    group_by(sample, sample_id, organism, protein_name) %>%
    summarize(prot_prop_Bayes = mean(beer_hits),
              prot_prop_edgeR = mean(edgeR_hits),
              num_peptides = n(), .groups = "drop") %>%
    arrange(desc(num_peptides)) %>%
    ggplot(aes(x = 0.5*(prot_prop_Bayes + prot_prop_edgeR),
               y = prot_prop_Bayes - prot_prop_edgeR,
               color = organism,
               size = num_peptides/2)) +
    facet_wrap(sample_id ~., ncol = 4) +
    geom_hline(aes(yintercept = 0), size = 0.5, color = "grey50") +
    geom_vline(aes(xintercept = 0), size = 0.5, color = "grey50") +
    geom_point() +
    labs(title = "Proportion of enriched peptides by protein",
         x = TeX("$\\frac{1}{2}$ (BEER + edgeR)"),
         y = "BEER - edgeR",
         color = "Organism",
         size = "# peptides") +
    theme_bw() +
    theme(legend.title = element_text(size = 10), 
          legend.text = element_text(size = 8), 
          legend.key.size = unit(0.75, "lines"), 
          legend.position = "bottom") +
    guides(color = guide_legend(ncol = 1, order = 1), 
           size = guide_legend(ncol = 1, order = 2))

ggsave("figures/coronascan_protein.png", units = "in", width = 7.5, height = 8)

# Figure: coronascan_replicates.png -----------
# Define sample to look at
cs_sample <- "CS 9"

# Plot for proportion of reads
prop_reads <- cs_tidy %>% 
    filter(sample_id == cs_sample) %>%
    mutate(prop_reads = counts/n) %>%
    select(pair_id, pair_num, prop_reads) %>%
    pivot_wider(names_from = "pair_num", 
                values_from = "prop_reads") %>%
    rename(pair_1 = `1`, pair_2 = `2`) %>%
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
          title = element_text(size = 8))

# Plot for BEER posterior probabilities
concordance_beer <- cs_tidy %>%
    filter(sample_id == cs_sample) %>%
    select(pair_id, pair_num, beer_hits) %>%
    pivot_wider(names_from = "pair_num",
                values_from = "beer_hits") %>%
    rename(pair_1 = `1`, pair_2 = `2`) %>%
    group_by(pair_1, pair_2) %>%
    summarize(num_pairs = n(), .groups = "drop") %>%
    bind_cols(xpos = c(-0.03, -0.03, 1.03, 1.03),
              ypos =  c(-0.05, 1.05,-0.05, 1.05),
              hjustvar = c(0, 0, 1, 1),
              vjustvar = c(0, 1, 0, 1)) %>%
    mutate(annotation = paste0(c("Both are not enriched (",
                                 "One is enriched (",
                                 "One is enriched (",
                                 "Both are enriched ("),
                                 num_pairs, rep(")", 4)))
    
post_prob <- cs_tidy %>% 
    filter(sample_id == cs_sample) %>%
    select(sample_id, pair_id, pair_num, beer_prob) %>%
    pivot_wider(names_from = "pair_num", 
                values_from = "beer_prob") %>%
    rename(pair_1 = `1`, pair_2 = `2`) %>%
    ggplot(aes(x = pair_1, y = pair_2)) + 
    geom_point() +
    labs(title = "BEER posterior probabilities", 
         x = "peptide pair #1", 
         y = "peptide pair #2") +
    geom_hline(aes(yintercept = 0.5), linetype = "dashed") +
    geom_vline(aes(xintercept = 0.5), linetype = "dashed") +
    geom_text(aes(x = xpos, y = ypos,
                  hjust = hjustvar, vjust = vjustvar,
                  label = annotation), color = "red", size = 3,
              data = concordance_beer) +
    coord_fixed(xlim = c(0, 1), ylim = c(0, 1)) +
    scale_x_continuous(breaks = seq(0, 1, by = 0.2), 
                       minor_breaks = seq(0, 1, by = 0.1), 
                       expand = c(0.035, 0.035)) +
    scale_y_continuous(breaks = seq(0, 1, by = 0.2), 
                       minor_breaks = seq(0, 1, by = 0.1), 
                       expand = c(0.035, 0.035)) +
    theme_bw() + 
    theme(title = element_text(size = 10), 
          aspect.ratio = 1)

# Plot for edgeR p-values
concordance_edgeR <- cs_tidy %>%
    filter(sample_id == cs_sample) %>%
    select(pair_id, pair_num, edgeR_hits) %>%
    pivot_wider(names_from = "pair_num",
                values_from = "edgeR_hits") %>%
    rename(pair_1 = `1`, pair_2 = `2`) %>%
    group_by(pair_1, pair_2) %>%
    summarize(num_pairs = n(), .groups = "drop") %>%
    bind_cols(xpos = c(-0.25, -0.25, 6.25, 6.25),
              ypos =  c(-0.25, 6.25, -0.25, 6.25),
              hjustvar = c(0, 0, 1, 1),
              vjustvar = c(0, 1, 0, 1)) %>%
    mutate(annotation = paste0(c("Both are not enriched (",
                                 "One is enriched (",
                                 "One is enriched (",
                                 "Both are enriched ("),
                               num_pairs, rep(")", 4)))

bh_cutoffs <- cs_tidy %>%
    filter(sample_id == cs_sample) %>%
    group_by(sample_id) %>%
    filter(edgeR_bh == max(edgeR_bh[edgeR_bh < 0.05])) %>%
    select(sample_id, edgeR_prob)
    
edgeR_pval <- cs_tidy %>% 
    filter(sample_id == cs_sample) %>%
    select(pair_id, pair_num, edgeR_prob) %>%
    pivot_wider(names_from = "pair_num", 
                values_from = "edgeR_prob") %>%
    rename(pair_1 = `1`, pair_2 = `2`) %>%
    ggplot(aes(x = pair_1, y = pair_2)) + 
    geom_point() +
    labs(title = "edgeR -log10(p-values)", 
         x = "peptide pair #1", 
         y = "peptide pair #2") +
    geom_hline(aes(yintercept = edgeR_prob),
               linetype = "dashed", data = bh_cutoffs) +
    geom_vline(aes(xintercept = edgeR_prob),
               linetype = "dashed", data = bh_cutoffs) +
    geom_text(aes(x = xpos, y = ypos,
                  hjust = hjustvar, vjust = vjustvar,
                  label = annotation), color = "red", size = 3,
              data = concordance_edgeR) +
    coord_fixed(xlim = c(-0.3, 6.3), ylim = c(-0.3, 6.3)) +
    scale_x_continuous(breaks = seq(0, 6, by = 1), 
                       minor_breaks = seq(0, 6, by = 0.5), 
                       expand = c(0.01, 0.01)) +
    scale_y_continuous(breaks = seq(0, 6, by = 1), 
                       minor_breaks = seq(0, 6, by = 0.5), 
                       expand = c(0.01, 0.01)) +
    theme_bw() +
    theme(title = element_text(size = 10), 
          aspect.ratio = 1)

# CAT curve for BEER posterior probabilities
#   For a fair comparison, set super-enriched peptides in the replicate sample
#   to have a rank of 1 as all of these peptides are identified as enriched
#   by BEER and edgeR. 
se_peps <- cs_tidy %>% 
    filter(sample_id == cs_sample & is_se & pair_num == 1) %>%
    mutate(rank = 1:n())

#   Ties are broken by the posterior marginal estimates of phi
beer_ranks <- cs_tidy %>%
    filter(sample_id == cs_sample) %>%
    group_by(pair_num) %>%
    arrange(desc(beer_prob), desc(beer_logfc), pair_id, .by_group = TRUE) %>%
    mutate(rank_se = nrow(se_peps) + cumsum(!pair_id %in% se_peps$pair_id)) %>%
    left_join(se_peps %>% select(pair_id, rank), by = c("pair_id")) %>%
    mutate(rank = ifelse(pair_id %in% se_peps$pair_id, rank, rank_se)) %>%
    select(pair_id, pair_num, rank) %>%
    pivot_wider(names_from = "pair_num", 
                values_from = "pair_id") %>%
    rename(pair_1 = `1`, pair_2 = `2`) %>%
    arrange(rank)

beer_conc <- sapply(1:nrow(beer_ranks), function(rank){
    length(intersect(beer_ranks$pair_1[beer_ranks$rank <= rank], 
                     beer_ranks$pair_2[beer_ranks$rank <= rank]))/rank
})

# CAT curve for edgeR p-values
edgeR_ranks <- cs_tidy %>%
    filter(sample_id == cs_sample) %>%
    group_by(pair_num) %>%
    arrange(edgeR_bh, desc(edgeR_logfc), pair_id, .by_group = TRUE) %>%
    mutate(rank_se = nrow(se_peps) + cumsum(!pair_id %in% se_peps$pair_id)) %>%
    left_join(se_peps %>% select(pair_id, rank), by = c("pair_id")) %>%
    mutate(rank = ifelse(pair_id %in% se_peps$pair_id, rank, rank_se)) %>%
    select(pair_id, pair_num, rank) %>%
    pivot_wider(names_from = "pair_num", 
                values_from = "pair_id") %>%
    rename(pair_1 = `1`, pair_2 = `2`) %>%
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
          legend.background = element_rect(color = "black", size = 0.3), 
          legend.position = c(0.8, 0.2))

ggarrange(prop_reads, cat_plot, post_prob, edgeR_pval, 
          nrow = 2, ncol = 2)

ggsave("figures/coronascan_replicates.png", units = "in", width = 8, height = 8)
