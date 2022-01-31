#' figure_attnconstant.R
#' 
#' Code to generate figures:
#' - Figure S16: attnconstant.png

# Set-up --------------
if(!"here" %in% installed.packages()){
    install.packages(here)
}
source(here::here("R", "load_packages.R"))

# Read in data
hiv <- readRDS(here::here("data_raw", "hiv_virscan.rds"))

# Convert to tidy format at add expected prop/counts
hiv_tidy <- as(hiv[, hiv$plate == 3], "DataFrame") %>%
    as_tibble() %>%
    group_by(sample) %>%
    mutate(n = sum(counts)) %>%
    ungroup() %>%
    select(sample, group, n, peptide, counts) %>%
    mutate(prop_reads = counts/n) %>%
    group_by(peptide) %>%
    mutate(expected_prop = mean(prop_reads[group == "beads"]), 
           expected_rc = expected_prop*n, 
           hit_2x = ifelse(counts > 2*expected_rc, TRUE, FALSE), 
           hit_5x = ifelse(counts > 5*expected_rc, TRUE, FALSE)) %>%
    ungroup()

## Identify one random sample on plate 3 for visualization purposes
set.seed(123)
sample_id <- sample(colnames(hiv[, hiv$group != "beads" & hiv$plate == 3]), 1)

# Figure S16: attnconstant.png ----------
hiv_tidy %>%
    filter(counts <= 250 & expected_rc <= 250 & sample == sample_id) %>%
    pivot_longer(cols = contains("hit_"), 
                 names_to = "threshold", 
                 values_to = "hit", 
                 names_pattern = "hit_([0-9])x") %>%
    ggplot(aes(x = expected_rc, y = counts, group = threshold)) + 
    facet_wrap(.~threshold) +
    geom_point(aes(color = hit), size = 0.5, alpha = 0.75) + 
    geom_abline(aes(intercept = 0, slope = 1), color = "black", size = 0.5) + 
    geom_smooth(formula = y ~ x - 1,
                data = hiv_tidy %>% 
                    filter(sample == sample_id) %>%
                    pivot_longer(cols = contains("hit_"), 
                                 names_to = "threshold", 
                                 values_to = "hit", 
                                 names_pattern = "hit_([0-9])x") %>%
                    filter(!hit),
                color = "blue", method = "lm", size = 0.5) +
    geom_text(aes(label = label), 
              x = 235, y = 10, size = 6, 
              #x = 57.75, y = 244, size = 5,
              data = data.frame(threshold = c("2", "5"), 
                                label = c("2X", "5X"))) +
    labs(x = "Expected # of reads", 
         y = "Observed # of reads", 
         color = "") +
    scale_color_manual(breaks = c("TRUE", "FALSE"), 
                       values = c("red", "grey70"), 
                       labels = c("enriched", "not enriched")) +
    coord_cartesian(xlim = c(0, 250), ylim = c(0, 250)) + 
    theme_bw() + 
    theme(aspect.ratio = 1, 
          strip.background = element_blank(), 
          strip.text = element_blank(),
          legend.position = "top", 
          legend.box.margin=margin(-10,-10,-10,-10))

ggsave("figures/attnconstant.png", units = "in", width = 6, height = 3.5,
       dpi = 600, bg = "white")
