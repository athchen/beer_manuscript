#' figure_hiv.R
#' 
#' Code to generate figures:
#' - hiv_protein.png
#' - hiv_replicates.png
#' - hiv_postprob.png
#' - hiv_pvalues.png

# Set-up --------------
source(file.path("R", "load_packages.R"))
source(file.path("R", "helper_functions.R"))

# Read in data
hiv <- readRDS(file.path("data_processed", "hiv_2e4.rds"))

sample_names <- c(paste0("HIV EC ", 1:(ncol(hiv) -1)), 
                  paste0("replicate of HIV EC ", ncol(hiv) - 1))
names(sample_names) <- colnames(hiv)

# Add peptide information, add hits for each approach
hiv_tidy <- as_df(hiv) %>%
    left_join(as_tibble(mcols(peptideInfo(hiv))), 
              by = c("peptide" = "pep_id")) %>%
    group_by(sample) %>%
    mutate(is_se = ifelse(sample != "beads" & is.na(beer_prob), TRUE, FALSE), 
           sample_id = sample_names[sample], 
           sample_id = factor(sample_id, levels = sample_names),
           hiv_species = gsub("Human immunodeficiency virus", 
                              "HIV", taxon_species),
           beer_hits = ifelse(beer_prob > 0.5 | is_se, 1, 0), 
           edgeR_bh = p.adjust(10^(-edgeR_prob), method = "BH"), 
           edgeR_hits = ifelse(edgeR_bh < 0.05, 1, 0)) %>%
    ungroup()

# Figure: hiv_protein.png ------------
hiv_tidy %>%
    filter(!grepl("BEADS", sample)) %>%
    select(sample, sample_id, peptide, UniProt_acc, hiv_species, beer_hits, edgeR_hits) %>%
    group_by(sample, sample_id, UniProt_acc, hiv_species) %>%
    summarize(prot_prop_Bayes = mean(beer_hits),
              prot_prop_edgeR = mean(edgeR_hits),
              num_peptides = n(), .groups = "drop") %>%
    arrange(desc(num_peptides)) %>%
    ggplot(aes(x = 0.5*(prot_prop_Bayes + prot_prop_edgeR), 
               y = prot_prop_Bayes - prot_prop_edgeR,
               color = hiv_species,
               size = num_peptides)) +
    facet_wrap(sample_id ~., ncol = 5) +
    geom_hline(aes(yintercept = 0), size = 0.5, color = "grey50") +
    geom_vline(aes(xintercept = 0), size = 0.5, color = "grey50") +
    geom_point() +
    scale_x_continuous(trans = "mysqrt", 
                       limits = c(0, 1), 
                       breaks = seq(0, 1, by = 0.2)) +
    labs(x = TeX("$\\frac{1}{2}$ (BEER + edgeR)"),
         y = "BEER - edgeR",
         color = "HIV strain",
         size = "# peptides") +
    scale_y_continuous(breaks = seq(-0.25, 1, by = 0.25), 
                       limits = c(-0.25, 1)) +
    theme_bw() +
    theme(legend.title = element_text(size = 10), 
          legend.text = element_text(size = 8), 
          legend.key.size = unit(0.75, "lines"), 
          legend.position = "bottom") +
    guides(size = guide_legend(ncol = 1), 
           color = guide_legend(ncol = 2))

ggsave("figures/hiv_protein.png", units = "in", width = 10.5, height = 6.5)

# Figure: hiv_replicates.png -----------
# Plot for proportion of reads
prop_reads <- hiv_tidy %>%
    mutate(prop_reads = counts/n) %>%
    filter(sample_id %in% c("HIV EC 15", "replicate of HIV EC 15")) %>%
    select(sample_id, peptide, prop_reads) %>%
    pivot_wider(names_from = "sample_id", 
                values_from = "prop_reads") %>%
    dplyr::rename(sample_1 = `HIV EC 15`, sample_2 = `replicate of HIV EC 15`) %>%
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
    filter(is_se & grepl("replicate", sample_id)) %>%
    mutate(rank = 1:n())

#   Ties are broken by the posterior marginal estimates of phi
beer_ranks <- hiv_tidy %>%
    filter(sample_id %in% c("HIV EC 15", "replicate of HIV EC 15")) %>%
    group_by(sample_id) %>%
    arrange(desc(beer_prob), desc(beer_logfc), peptide, .by_group = TRUE) %>%
    mutate(rank_se = nrow(se_peps) + cumsum(!peptide %in% se_peps$peptide)) %>%
    left_join(se_peps %>% select(peptide, rank), by = c("peptide")) %>%
    mutate(rank = ifelse(peptide %in% se_peps$peptide, rank, rank_se)) %>%
    ungroup() %>%
    select(sample_id, peptide, rank) %>%
    pivot_wider(names_from = "sample_id", 
                values_from = "peptide") %>%
    dplyr::rename(sample_1 = `HIV EC 15`, sample_2 = `replicate of HIV EC 15`) %>%
    arrange(rank)

beer_conc <- sapply(1:nrow(beer_ranks), function(rank){
    length(intersect(beer_ranks$sample_1[beer_ranks$rank <= rank], 
                     beer_ranks$sample_2[beer_ranks$rank <= rank]))/rank
})

# CAT curve for edgeR p-values
edgeR_ranks <- hiv_tidy %>%
    filter(sample_id %in% c("HIV EC 15", "replicate of HIV EC 15")) %>%
    group_by(sample_id) %>%
    arrange(edgeR_bh, desc(edgeR_logfc), .by_group = TRUE) %>%
    mutate(rank_se = 8 + cumsum(!peptide %in% se_peps$peptide)) %>%
    left_join(se_peps %>% select(peptide, rank), by = c("peptide")) %>%
    mutate(rank = ifelse(peptide %in% se_peps$peptide, rank, rank_se)) %>%
    ungroup() %>%
    select(sample_id, peptide, rank) %>%
    pivot_wider(names_from = "sample_id", 
                values_from = "peptide") %>%
    dplyr::rename(sample_1 = `HIV EC 15`, sample_2 = `replicate of HIV EC 15`) %>%
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

# Figure: hiv_postprob.png -----------
hiv_tidy %>%
    filter(sample_id %in% c("HIV EC 15", "replicate of HIV EC 15")) %>%
    select(sample_id, peptide, beer_prob) %>%
    group_by(sample_id) %>%
    arrange(desc(beer_prob), .by_group = TRUE) %>%
    mutate(rank = 1:n()) %>%
    ggplot(aes(x = rank, y = beer_prob, group = sample_id, color = sample_id)) + 
    geom_point(size = 1) +
    geom_line() +
    coord_cartesian(xlim = c(0, 400), ylim = c(0, 1)) +
    labs(x = "rank", 
         y = "posterior probabilities", 
         color = "sample") +
    scale_color_manual(values = c("red", "black")) +
    theme_bw() + 
    theme(aspect.ratio = 1, 
          title = element_text(size = 10), 
          legend.background = element_rect(color = "black", size = 0.3), 
          legend.position = c(0.725, 0.825))

ggsave("figures/hiv_postprob.png", units = "in", width = 4, height = 4)

# Figure: hiv_pvalues.png -----------
hiv_tidy %>%
    filter(sample_id %in% c("HIV EC 15", "replicate of HIV EC 15")) %>%
    select(sample_id, peptide, edgeR_prob) %>%
    group_by(sample_id) %>%
    arrange(desc(edgeR_prob), .by_group = TRUE) %>%
    mutate(rank = 1:n()) %>%
    ggplot(aes(x = rank, y = edgeR_prob, group = sample_id, color = sample_id)) + 
    geom_point(size = 1) +
    geom_line() +
    labs(x = "rank", 
         y = "-log10(p-values)", 
         color = "sample") +
    coord_cartesian(xlim = c(0, 400)) +
    scale_color_manual(values = c("red", "black")) +
    theme_bw() +
    theme(aspect.ratio = 1, 
          title = element_text(size = 10), 
          legend.background = element_rect(color = "black", size = 0.3), 
          legend.position = c(0.725, 0.825))

ggsave("figures/hiv_pvalues.png", units = "in", width = 4, height = 4)
