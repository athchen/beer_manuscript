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

# Plot for BEER posterior probabilities
concordance_beer <- hiv_tidy %>%
    filter(sample_id %in% c("HIV EC 15", "replicate of HIV EC 15")) %>%
    select(sample_id, peptide, beer_hits) %>%
    pivot_wider(names_from = "sample_id", 
                values_from = "beer_hits") %>%
    rename(sample_1 = `HIV EC 15`, sample_2 = `replicate of HIV EC 15`) %>%
    group_by(sample_1, sample_2) %>%
    summarize(num_peps = n(), .groups = "drop") %>%
    bind_cols(xpos = c(-0.03, -0.03, 1.03, 1.03),
              ypos =  c(-0.05, 1.05,-0.05, 1.05),
              hjustvar = c(0, 0, 1, 1),
              vjustvar = c(0, 1, 0, 1)) %>%
    mutate(annotation = paste0(c("Not enriched in both (", 
                                 "Not enriched in one (", 
                                 "Not enriched in one (", 
                                 "Enriched in both ("),
                                 num_peps, rep(")", 4)))
    
post_prob <- hiv_tidy %>%
    filter(sample_id %in% c("HIV EC 15", "replicate of HIV EC 15")) %>%
    select(sample_id, peptide, beer_prob) %>%
    pivot_wider(names_from = "sample_id", 
                values_from = "beer_prob") %>%
    rename(sample_1 = `HIV EC 15`, sample_2 = `replicate of HIV EC 15`) %>%
    ggplot(aes(x = sample_1, y = sample_2)) + 
    geom_point() +
    labs(title = "BEER posterior probabilities", 
         x = "HIV EC 15", 
         y = "replicate of HIV EC 15") +
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
    theme(title = element_text(size = 10))

# Plot for edgeR p-values
concordance_edgeR <- hiv_tidy %>%
    filter(sample_id %in% c("HIV EC 15", "replicate of HIV EC 15")) %>%
    select(sample_id, peptide, edgeR_hits) %>%
    pivot_wider(names_from = "sample_id", 
                values_from = "edgeR_hits") %>%
    rename(sample_1 = `HIV EC 15`, sample_2 = `replicate of HIV EC 15`) %>%
    group_by(sample_1, sample_2) %>%
    summarize(num_peps = n(), .groups = "drop") %>%
    bind_cols(xpos = c(-0.25, -0.25, 6.25, 6.25),
              ypos =  c(-0.25, 6.25, -0.25, 6.25),
              hjustvar = c(0, 0, 1, 1),
              vjustvar = c(0, 1, 0, 1)) %>%
    mutate(annotation = paste0(c("Not enriched in both (", 
                                 "Not enriched in one (", 
                                 "Not enriched in one (", 
                                 "Enriched in both ("),
                               num_peps, rep(")", 4)))
bh_cutoffs <- hiv_tidy %>%
    filter(sample_id %in% c("HIV EC 15", "replicate of HIV EC 15")) %>%
    group_by(sample_id) %>%
    filter(edgeR_bh == max(edgeR_bh[edgeR_bh < 0.05])) %>%
    select(sample_id, edgeR_prob)
    
edgeR_pval <- hiv_tidy %>%
    filter(sample_id %in% c("HIV EC 15", "replicate of HIV EC 15")) %>%
    select(sample_id, peptide, edgeR_prob) %>%
    pivot_wider(names_from = "sample_id", 
                values_from = "edgeR_prob") %>%
    rename(sample_1 = `HIV EC 15`, sample_2 = `replicate of HIV EC 15`) %>%
    ggplot(aes(x = sample_1, y = sample_2)) + 
    geom_point() +
    labs(title = "edgeR -log10(p-values)", 
         x = "HIV EC 15", 
         y = "replicate of HIV EC 15") +
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
    theme(title = element_text(size = 10))

ggarrange(prop_reads, post_prob, edgeR_pval, nrow = 1, widths = c(1.035, 1, 1))

ggsave("figures/hiv_replicates.png", units = "in", width = 12, height = 4)
