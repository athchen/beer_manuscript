#' figure_hiv.R
#' 
#' Code to generate figures:
#' - hiv_protein.png
#' - hiv_protein_noB.png
#' - hiv_protein_A.png
#' - hiv_replicates.png
#' - hiv_ranked_prob.png

# Set-up --------------
source(file.path("R", "load_packages.R"))
source(here("R", "helper_functions.R"))

# Read in data
hiv <- readRDS(here::here("data_processed", "hiv_results.rds"))

# Convert to tidy format and add hits for each approach
hiv_tidy <- as(hiv, "DataFrame") %>%
    as_tibble() %>%
    group_by(sample) %>%
    mutate(is_se = ifelse(sample != "beads" & is.na(beer_prob), TRUE, FALSE), 
           sample = factor(sample, levels = colnames(hiv)),
           beer_hits = ifelse(beer_prob > 0.5 | is_se, 1, 0), 
           edgeR_bh = p.adjust(10^(-edgeR_prob), method = "BH"), 
           edgeR_hits = ifelse(edgeR_bh < 0.05, 1, 0)) %>%
    ungroup()

# Figure: hiv_protein.png ------------
# Color palette
hiv_subtype <- unique(hiv_tidy$taxon_species)
grey_palette <- palette(gray(seq(0.1, 0.8, len = (length(hiv_subtype) - 1))))
num_B <- grep("HIV type 1 group M subtype B", hiv_subtype)
subtype_order <- c(hiv_subtype[num_B], hiv_subtype[-num_B])

hiv_tidy %>%
    filter(group != "beads" & !grepl("9", sample)) %>%
    select(sample, peptide, UniProt_acc, taxon_species,
           beer_hits, edgeR_hits) %>%
    group_by(sample, UniProt_acc, taxon_species) %>%
    dplyr::summarize(prot_prop_Bayes = mean(beer_hits),
              prot_prop_edgeR = mean(edgeR_hits),
              num_peptides = n(), .groups = "drop") %>%
    mutate(taxon_species = factor(taxon_species, subtype_order)) %>%
    arrange(desc(taxon_species), desc(num_peptides)) %>%
    ggplot(aes(x = 0.5*(prot_prop_Bayes + prot_prop_edgeR), 
               y = prot_prop_Bayes - prot_prop_edgeR,
               color = taxon_species,
               size = num_peptides)) +
    facet_wrap(sample ~., ncol = 4) +
    geom_hline(aes(yintercept = 0), size = 0.5, color = "grey50") +
    geom_vline(aes(xintercept = 0), size = 0.5, color = "grey50") +
    geom_point(alpha = 0.8) +
    labs(x = TeX("$\\frac{1}{2}$ (BEER + edgeR)"),
         y = "BEER - edgeR",
         color = "HIV strain",
         size = "# peptides") +
    scale_x_continuous(trans = "mysqrt", 
                       limits = c(0, 1), 
                       breaks = seq(0, 1, by = 0.2)) +
    scale_y_continuous(breaks = seq(-0.25, 1, by = 0.25), 
                       limits = c(-0.25, 1)) +
    scale_color_manual(values = c( "firebrick2", grey_palette),
                       breaks = subtype_order) +
    theme_bw() +
    theme(legend.title = element_text(size = 10), 
          legend.text = element_text(size = 8), 
          legend.key.size = unit(0.75, "lines"), 
          legend.position = "bottom") +
    guides(size = guide_legend(ncol = 1), 
           color = guide_legend(ncol = 2))

ggsave("figures/hiv_protein.png", units = "in", width = 8.5, height = 6.5)

# Figure: hiv_protein_noB.png ----------
hiv_tidy %>%
    filter(group != "beads" & !grepl("9", sample) &
               !grepl("HIV type 1 group M subtype B", taxon_species)) %>%
    select(sample, peptide, UniProt_acc, taxon_species,
           beer_hits, edgeR_hits) %>%
    group_by(sample, UniProt_acc, taxon_species) %>%
    summarize(prot_prop_Bayes = mean(beer_hits),
              prot_prop_edgeR = mean(edgeR_hits),
              num_peptides = n(), .groups = "drop") %>%
    mutate(taxon_species = factor(taxon_species, subtype_order)) %>%
    arrange(desc(taxon_species), desc(num_peptides)) %>%
    ggplot(aes(x = 0.5*(prot_prop_Bayes + prot_prop_edgeR), 
               y = prot_prop_Bayes - prot_prop_edgeR,
               color = taxon_species,
               size = num_peptides)) +
    facet_wrap(sample ~., ncol = 4) +
    geom_hline(aes(yintercept = 0), size = 0.5, color = "grey50") +
    geom_vline(aes(xintercept = 0), size = 0.5, color = "grey50") +
    geom_point(alpha = 0.8) +
    labs(x = TeX("$\\frac{1}{2}$ (BEER + edgeR)"),
         y = "BEER - edgeR",
         color = "HIV strain",
         size = "# peptides") +
    scale_x_continuous(trans = "mysqrt", 
                       limits = c(0, 1), 
                       breaks = seq(0, 1, by = 0.2)) +
    scale_y_continuous(breaks = seq(-0.25, 1, by = 0.25), 
                       limits = c(-0.25, 1)) +
    scale_color_manual(values = c( "firebrick2", grey_palette),
                       breaks = subtype_order) +
    theme_bw() +
    theme(legend.title = element_text(size = 10), 
          legend.text = element_text(size = 8), 
          legend.key.size = unit(0.75, "lines"), 
          legend.position = "bottom") +
    guides(size = guide_legend(ncol = 1), 
           color = guide_legend(ncol = 2))

ggsave("figures/hiv_protein_noB.png", units = "in", width = 8.5, height = 6.5)

# Figure: hiv_protein_A.png ----------
# Color palette
num_A <- grep("HIV type 1 group M subtype A", hiv_subtype)
subtype_order <- c(hiv_subtype[c(num_B, num_A)], hiv_subtype[-c(num_B, num_A)])

hiv_tidy %>%
    filter(grepl("9", sample)) %>%
    select(sample, peptide, UniProt_acc, taxon_species,
           beer_hits, edgeR_hits) %>%
    group_by(sample, UniProt_acc, taxon_species) %>%
    summarize(prot_prop_Bayes = mean(beer_hits),
              prot_prop_edgeR = mean(edgeR_hits),
              num_peptides = n(), .groups = "drop") %>%
    mutate(taxon_species = factor(taxon_species, subtype_order)) %>%
    arrange(desc(taxon_species), desc(num_peptides)) %>%
    ggplot(aes(x = 0.5*(prot_prop_Bayes + prot_prop_edgeR), 
               y = prot_prop_Bayes - prot_prop_edgeR,
               color = taxon_species,
               size = num_peptides)) +
    facet_wrap(sample ~., ncol = 4) +
    geom_hline(aes(yintercept = 0), size = 0.5, color = "grey50") +
    geom_vline(aes(xintercept = 0), size = 0.5, color = "grey50") +
    geom_point(alpha = 0.8) +
    labs(x = TeX("$\\frac{1}{2}$ (BEER + edgeR)"),
         y = "BEER - edgeR",
         color = "HIV strain",
         size = "# peptides") +
    scale_x_continuous(trans = "mysqrt", 
                       limits = c(0, 1), 
                       breaks = seq(0, 1, by = 0.2)) +
    scale_y_continuous(breaks = seq(-0.25, 1, by = 0.25), 
                       limits = c(-0.25, 1)) +
    scale_color_manual(values = c("firebrick2", "blue", grey_palette[-num_A]),
                       breaks = subtype_order) +
    theme_bw() +
    theme(legend.title = element_text(size = 10), 
          legend.text = element_text(size = 8), 
          legend.key.size = unit(0.75, "lines"), 
          legend.position = "bottom") +
    guides(size = guide_legend(ncol = 1), 
           color = guide_legend(ncol = 2))

ggsave("figures/hiv_protein_A.png", units = "in", width = 8.5, height = 6.5)

# Figure: hiv_replicates.png -----------
# Plot for proportion of reads
prop_reads <- hiv_tidy %>%
    mutate(prop_reads = counts/n) %>%
    filter(sample %in% c("HIV EC 9", "replicate of HIV EC 9")) %>%
    select(sample, peptide, prop_reads) %>%
    pivot_wider(names_from = "sample", 
                values_from = "prop_reads") %>%
    dplyr::rename(sample_1 = `HIV EC 9`, sample_2 = `replicate of HIV EC 9`) %>%
    ggplot(aes(x = log10(sample_1), y = log10(sample_2))) + 
    geom_abline(aes(intercept = 0, slope = 1)) +
    geom_point() +
    labs(title = "Proportion of reads pulled", 
         x = "sample, log10(proportion)", 
         y = "replicate, log10(proportion)") +
    scale_x_continuous(breaks = seq(-5, -1, by = 1),
                       minor_breaks = seq(-5, -1, by = 0.5), 
                       labels = TeX(paste0("$10^{", -5:-1, "}$"))) +
    scale_y_continuous(breaks = seq(-5, -1, by = 1),
                       minor_breaks = seq(-5, -1, by = 0.5), 
                       labels = TeX(paste0("$10^{", -5:-1, "}$"))) +
    theme_bw() +
    theme(aspect.ratio = 1, 
          title = element_text(size = 10))

# CAT curve for BEER posterior probabilities
#   For a fair comparison, set super-enriched peptides in the replicate sample
#   to have a rank of 1 as all of these peptides are identified as enriched
#   by BEER and edgeR. 
se_peps <- hiv_tidy %>% 
    filter(is_se & grepl("9", sample)) %>%
    group_by(peptide) %>%
    filter(row_number() == 1) %>%
    ungroup() %>%
    mutate(rank = 1:n())

#   Ties are broken by the posterior marginal estimates of phi
beer_ranks <- hiv_tidy %>%
    filter(sample %in% c("HIV EC 9", "replicate of HIV EC 9")) %>%
    group_by(sample) %>%
    arrange(desc(beer_prob), desc(beer_logfc), peptide, .by_group = TRUE) %>%
    mutate(rank_se = nrow(se_peps) + cumsum(!peptide %in% se_peps$peptide)) %>%
    left_join(se_peps %>% select(peptide, rank), by = c("peptide")) %>%
    mutate(rank = ifelse(peptide %in% se_peps$peptide, rank, rank_se)) %>%
    ungroup() %>%
    select(sample, peptide, rank) %>%
    pivot_wider(names_from = "sample", 
                values_from = "peptide") %>%
    dplyr::rename(sample_1 = `HIV EC 9`, sample_2 = `replicate of HIV EC 9`) %>%
    arrange(rank)

beer_conc <- sapply(1:nrow(beer_ranks), function(rank){
    length(intersect(beer_ranks$sample_1[beer_ranks$rank <= rank], 
                     beer_ranks$sample_2[beer_ranks$rank <= rank]))/rank
})

# CAT curve for edgeR p-values
edgeR_ranks <- hiv_tidy %>%
    filter(sample %in% c("HIV EC 9", "replicate of HIV EC 9")) %>%
    group_by(sample) %>%
    arrange(edgeR_bh, desc(edgeR_logfc), .by_group = TRUE) %>%
    mutate(rank_se = 8 + cumsum(!peptide %in% se_peps$peptide)) %>%
    left_join(se_peps %>% select(peptide, rank), by = c("peptide")) %>%
    mutate(rank = ifelse(peptide %in% se_peps$peptide, rank, rank_se)) %>%
    ungroup() %>%
    select(sample, peptide, rank) %>%
    pivot_wider(names_from = "sample", 
                values_from = "peptide") %>%
    dplyr::rename(sample_1 = `HIV EC 9`, sample_2 = `replicate of HIV EC 9`) %>%
    arrange(rank)

edgeR_conc <- sapply(1:nrow(edgeR_ranks), function(rank){
    length(intersect(edgeR_ranks$sample_1[edgeR_ranks$rank <= rank], 
                     edgeR_ranks$sample_2[edgeR_ranks$rank <= rank]))/rank
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
         y = "concordance", 
         x = "rank") +
    coord_cartesian(ylim = c(0, 1), xlim = c(0, 200)) +
    scale_x_continuous(breaks = seq(0, 200, by = 25)) +
    scale_y_continuous(breaks = seq(0, 1, by = 0.2)) +
    scale_color_manual(values = c("red", "black")) +
    theme_bw() + 
    theme(aspect.ratio = 1, 
          title = element_text(size = 10),
          legend.background = element_rect(color = "black", size = 0.3), 
          legend.position = c(0.8, 0.2))

ggarrange(prop_reads, cat_plot, align = "hv", 
          nrow = 1, ncol = 2)

ggsave("figures/hiv_replicates.png", units = "in", width = 8, height = 4)

# Figure: hiv_ranked_prob.png -----------
post_prob <- hiv_tidy %>%
    filter(sample %in% c("HIV EC 9", "replicate of HIV EC 9")) %>%
    select(sample, peptide, beer_prob) %>%
    group_by(sample) %>%
    arrange(desc(beer_prob), .by_group = TRUE) %>%
    mutate(rank = 1:n()) %>%
    ggplot(aes(x = rank, y = beer_prob, group = sample, color = sample)) + 
    geom_point(size = 1) +
    geom_line() +
    coord_cartesian(xlim = c(0, 400), ylim = c(0, 1)) +
    labs(title = "BEER", 
         x = "rank", 
         y = "posterior probabilities", 
         color = "sample") +
    scale_color_manual(values = c("red", "black")) +
    theme_bw() + 
    theme(aspect.ratio = 1, 
          title = element_text(size = 10), 
          legend.background = element_rect(color = "black", size = 0.3), 
          legend.position = "none")

p_values <- hiv_tidy %>%
    filter(sample %in% c("HIV EC 9", "replicate of HIV EC 9")) %>%
    select(sample, peptide, edgeR_prob) %>%
    group_by(sample) %>%
    arrange(desc(edgeR_prob), .by_group = TRUE) %>%
    mutate(rank = 1:n()) %>%
    ggplot(aes(x = rank, y = edgeR_prob, group = sample, color = sample)) + 
    geom_point(size = 1) +
    geom_line() +
    labs(title = "edgeR", 
         x = "rank", 
         y = "-log10(p-values)", 
         color = "sample") +
    coord_cartesian(xlim = c(0, 400), ylim = c(0, 6)) +
    scale_color_manual(values = c("red", "black")) +
    theme_bw() +
    theme(aspect.ratio = 1, 
          title = element_text(size = 10), 
          legend.background = element_rect(color = "black", size = 0.3), 
          legend.position = c(0.725, 0.825))

ggarrange(post_prob, p_values, nrow = 1)

ggsave("figures/hiv_ranked_prob.png", units = "in", width = 8, height = 4, 
       bg = "white")

# Table: hiv_replicates -----------
hiv_tidy %>%
    filter(sample %in% c("HIV EC 9", "replicate of HIV EC 9")) %>%
    mutate(subtype_a = ifelse(taxon_species == "HIV type 1 group M subtype A", 1, 0)) %>%
    select(subtype_a, sample, peptide, beer_hits, edgeR_hits) %>%
    pivot_longer(cols = contains("hits"), 
                 names_to = "method", 
                 names_pattern = "([a-zA-Z]*)_hits", 
                 values_to = "hits") %>%
    pivot_wider(names_from = sample, values_from = hits) %>%
    dplyr::rename(sample_1 = `HIV EC 9`, 
                  sample_2 = `replicate of HIV EC 9`) %>%
    mutate(concordance = case_when(sample_1 == 1 & sample_2 == 1 ~ "both enriched", 
                                   sample_1 == 0 & sample_2 == 0 ~ "both not enriched", 
                                   TRUE ~ "discordant")) %>%
    group_by(subtype_a, method, concordance) %>%
    summarize(n = n(), .groups = "drop_last") %>%
    mutate(total_peps = sum(n), 
           proportion = n/total_peps) %>%
    pivot_longer(cols = c("n", "proportion"), 
                 names_to = "param_name", 
                 values_to = "value") %>%
    unite("param", method, param_name) %>%
    pivot_wider(names_from = "param", values_from = "value") %>%
    arrange(desc(subtype_a))
