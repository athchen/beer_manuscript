#' figure_hiv.R
#' 
#' Code to generate figures:
#' - hiv_protein.png
#' - hiv_replicates.png

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
    facet_wrap(sample_id ~., ncol = 4) +
    geom_hline(aes(yintercept = 0), size = 0.5, color = "grey50") +
    geom_vline(aes(xintercept = 0), size = 0.5, color = "grey50") +
    geom_point() +
    labs(title = "Proportion of enriched peptides by protein",
         x = TeX("$\\frac{1}{2}$ (BEER + edgeR)"),
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

ggsave("figures/hiv_protein.png", units = "in", width = 7.5, height = 7.5)

# Figure: hiv_replicates.png -----------
# Plot for proportion of reads
prop_reads <- hiv_tidy %>%
    mutate(prop_reads = counts/n) %>%
    filter(sample_id %in% c("HIV EC 15", "replicate of HIV EC 15")) %>%
    select(sample_id, peptide, prop_reads) %>%
    pivot_wider(names_from = "sample_id", 
                values_from = "prop_reads") %>%
    rename(sample_1 = `HIV EC 15`, sample_2 = `replicate of HIV EC 15`) %>%
    ggplot(aes(x = 0.5*(sample_1 + sample_2), 
               y = sample_1 - sample_2)) + 
    geom_point() +
    labs(title = "Agreement in the proportion of reads pulled", 
         x = "average proportion of reads pulled", 
         y = "difference in proportion of reads pulled") +
    scale_x_continuous(breaks = seq(0, 0.03, by = 0.005), 
                       minor_breaks = seq(0, 0.03, by = 0.001)) +
    scale_y_continuous(breaks = seq(-0.006, 0.001, by = 0.001), 
                       minor_breaks = seq(-0.006, 0.001, by = 0.0005)) +
    theme_bw() +
    theme(aspect.ratio = 1, 
          title = element_text(size = 10))

# CAT curve for BEER posterior probabilities
# Ties are broken by the posterior marginal estimates of phi
beer_ranks <- hiv_tidy %>%
    filter(sample_id %in% c("HIV EC 15", "replicate of HIV EC 15")) %>%
    mutate(beer_prob = ifelse(is.na(beer_prob), 1, beer_prob), 
           beer_logfc = ifelse(is.na(beer_prob), Inf, beer_logfc)) %>%
    group_by(sample_id) %>%
    arrange(desc(beer_prob), desc(beer_logfc), .by_group = TRUE) %>%
    mutate(rank = 1:n()) %>%
    ungroup() %>%
    select(sample_id, peptide, rank) %>%
    pivot_wider(names_from = "sample_id", 
                values_from = "peptide") %>%
    rename(sample_1 = `HIV EC 15`, sample_2 = `replicate of HIV EC 15`) 

beer_conc <- sapply(1:nrow(beer_ranks), function(rank){
    length(intersect(beer_ranks$sample_1[1:rank], beer_ranks$sample_2[1:rank]))/rank
})

# CAT curve for edgeR p-values
edgeR_ranks <- hiv_tidy %>%
    filter(sample_id %in% c("HIV EC 15", "replicate of HIV EC 15")) %>%
    group_by(sample_id) %>%
    arrange(edgeR_bh, desc(edgeR_logfc), .by_group = TRUE) %>%
    mutate(rank = 1:n()) %>%
    ungroup() %>%
    select(sample_id, peptide, rank) %>%
    pivot_wider(names_from = "sample_id", 
                values_from = "peptide") %>%
    rename(sample_1 = `HIV EC 15`, sample_2 = `replicate of HIV EC 15`) 

edgeR_conc <- sapply(1:nrow(edgeR_ranks), function(rank){
    length(intersect(edgeR_ranks$sample_1[1:rank], edgeR_ranks$sample_2[1:rank]))/rank
})

cat_plot <- data.frame(rank = 1:length(beer_conc), 
           BEER = beer_conc, 
           edgeR = edgeR_conc) %>%
    pivot_longer(cols = c("BEER", "edgeR"), 
                 names_to = "approach", 
                 values_to = "concordance") %>%
    ggplot(aes(x = rank, y = concordance, color = approach)) +
    geom_line() +
    labs(title = "Concordance between technical replicates", 
         y = "Concordance", 
         x = "Rank") +
    coord_cartesian(ylim = c(0, 1)) +
    # scale_x_continuous(breaks = seq(0, 150, by = 25)) +
    scale_y_continuous(breaks = seq(0, 1, by = 0.2)) +
    scale_color_manual(values = c("red", "black")) +
    theme_bw() + 
    theme(aspect.ratio = 1)

ggarrange(prop_reads, cat_plot, nrow = 1, widths = c(1, 1.25))

ggsave("figures/hiv_replicates.png", units = "in", width = 9, height = 4)
