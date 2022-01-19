#' figure_corstructure.R
#' 
#' Code to generate figure:
#' - Figure S13: corstructure.png
#' - Figure S14: figdatastructurebinom.png

# Set-up --------------
if(!"here" %in% installed.packages()){
    install.packages(here)
}
source(here::here("R", "load_packages.R"))

# Read in data
hiv_virscan <- readRDS(here::here("data_raw", "hiv_virscan.rds"))
hiv_beads <- hiv_virscan[, hiv_virscan$group == "beads"]

# Figure S13: corstructure.png ----------
## Left panel: peptide effect in beads-only samples
cpm <- propReads(hiv_beads)*1e6
set.seed(1)
breaks <- c(1, 2, 5, 10, 20, 50, 100)

plot_pepeff <- data.frame(beads_s1p3 = cpm[, "beads_3_C9"], 
                          beads_s1p4 = cpm[, "beads_4_A5"]) %>%
    ggplot(aes(x = log10(beads_s1p3 + 1), 
               y = log10(beads_s1p4 + 1))) + 
    geom_jitter(width = 0.01, height = 0.01, size = 0.1) +
    geom_abline(aes(slope = 1, intercept = 0), color = "red") +
    scale_x_continuous(breaks = log10(breaks), 
                       labels = breaks) +
    scale_y_continuous(breaks = log10(breaks), 
                       labels = breaks) +
    coord_fixed(xlim = c(0, log10(100)), ylim = c(0, log10(100))) +
    labs(x = "reads per million + 1 [sample 1 plate 3]", 
         y = "reads per million + 1 [sample 1 plate 4]") +
    theme_bw() + 
    theme(panel.grid = element_blank(), 
          aspect.ratio = 1)

# Spearman correlation 
cor(cpm[, "beads_3_C9"], cpm[, "beads_4_A5"], 
    method = "spearman")

## Middle panel: average cpm
avg_cpm <- vapply(paste0("plate", 1:5), function(plate){
    plate_num <- str_match(plate, "plate([0-9])")[, 2]
    rowMeans(cpm[, grepl(paste0("beads_", plate_num), colnames(cpm))])
}, numeric(nrow(cpm)))

plot_avgcpm <- avg_cpm %>%
    as_tibble() %>%
    ggplot(aes(x = log10(plate2 + 1), y = log10(plate3 + 1))) + 
    geom_jitter(width = 0.01, height = 0.01, size = 0.1) +
    geom_abline(aes(slope = 1, intercept = 0), color = "red") +
    scale_x_continuous(breaks = log10(breaks), 
                       labels = breaks) +
    scale_y_continuous(breaks = log10(breaks), 
                       labels = breaks) +
    coord_fixed(xlim = c(0, log10(100)), ylim = c(0, log10(100))) +
    labs(x = "average reads per million + 1 [plate 2]", 
         y = "average reads per million + 1 [plate 3]") +
    theme_bw() + 
    theme(panel.grid = element_blank(), 
          aspect.ratio = 1)

# Spearman correlation
cor(avg_cpm[, "plate2"], avg_cpm[, "plate3"], method = "spearman")   

## Right panel: correlation by read depth
# correlations between samples within plates
within_pl <- cor(cpm, method = "spearman") 
within_pl[lower.tri(within_pl, diag = TRUE)] <- NA
within_df <- within_pl %>%
    as_tibble(rownames = "sample_1") %>%
    pivot_longer(cols = -sample_1, 
                 names_to = "sample_2", 
                 values_to = "corr") %>%
    mutate(pl_1 = str_match(sample_1, "beads_([0-9])")[, 2], 
           pl_2 = str_match(sample_2, "beads_([0-9])")[, 2]) %>%
    left_join(data.frame(sample_id = colnames(hiv_beads), 
                         rd_1 = librarySize(hiv_beads)),
              by = c("sample_1"= "sample_id")) %>%
    left_join(data.frame(sample_id = colnames(hiv_beads), 
                         rd_2 = librarySize(hiv_beads)),
              by = c("sample_2"= "sample_id")) %>%
    filter(!is.na(corr) & pl_1 == pl_2) %>%
    rowwise() %>%
    mutate(min_rd = min(rd_1, rd_2)) %>%
    ungroup()

# correlations between plate averages
med_rd <- vapply(paste0("plate", 1:5), function(plate){
    plate_num <- str_match(plate, "plate([0-9])")[, 2]
    median(librarySize(hiv_beads)[grepl(paste0("beads_", plate_num), colnames(cpm))])
}, numeric(1))

between_pl <- cor(avg_cpm, method = "spearman")
between_pl[lower.tri(between_pl, diag = TRUE)] <- NA
between_df <- between_pl %>%
    as_tibble(rownames = "plate_1") %>%
    pivot_longer(cols = -plate_1, 
                 names_to = "plate_2", 
                 values_to = "corr") %>%
    left_join(data.frame(plate = names(med_rd), 
                         med_rd_1 = med_rd),
              by = c("plate_1"= "plate")) %>%
    left_join(data.frame(plate = names(med_rd), 
                         med_rd_2 = med_rd),
              by = c("plate_2"= "plate")) %>%
    filter(!is.na(corr)) %>%
    rowwise() %>%
    mutate(med_rd = min(med_rd_1, med_rd_2)) %>%
    ungroup() 

plot_corr <- ggplot() +
    geom_jitter(aes(x = min_rd/1e6, y = corr), 
                data = within_df,  width = 0.005, color = "red") +
    geom_jitter(aes(x = med_rd/1e6, y = corr),
               data = between_df, width = 0.005, color = "blue") +
    scale_x_continuous(breaks = seq(1.5, 4.0, by = 0.5)) +
    scale_y_continuous(breaks = seq(0.5, 1.0, by = 0.1)) +
    coord_cartesian(xlim = c(1.5, 4.0), ylim = c(0.5, 1.0)) +
    labs(x = "sample minimum read depth [million]", 
         y = "spearman correlation") +
    theme_bw() + 
    theme(aspect.ratio = 1)

ggarrange(plot_pepeff, plot_avgcpm, plot_corr, ncol = 3, 
          align = "v")    

ggsave("figures/corstructure.png", units = "in", width = 12, height = 4.5, 
       bg = "white")

# Figure S14: figdatastructurebinom.png ------------
# Plot library size
lib_size <- data.frame(libsize = librarySize(hiv_beads)) %>%
    rownames_to_column("sample") %>%
    mutate(plate = str_match(sample, "beads_([0-9])_")[, 2], 
           plate = factor(plate, levels = 1:5)) %>%
    arrange(plate, sample) %>%
    ggplot(aes(x = sample, y = libsize, fill = plate)) + 
    geom_bar(stat = "identity", color = "black", show.legend = FALSE) + 
    geom_point(shape = 21, size = 2, color = "black", 
               data = data.frame(sample = colnames(hiv_beads), 
                                 plate = factor(hiv_beads$plate, levels = 1:5), 
                                 libsize = rep(NA_integer_, ncol(hiv_beads)))) +
    scale_fill_manual(values = brewer.pal(5,"Accent"), 
                      labels = paste("Plate", 1:5)) +
    scale_y_continuous(breaks = (0:4)*1e6, 
                       labels = 0:4, 
                       expand = expansion(mult = c(0, 0.02))) +
    scale_x_discrete(labels = NULL, 
                     expand = expansion(add = c(1, 0))) +
    labs(y = "reads [ millions ]", 
         x = "") +
    guides(fill = guide_legend(keywidth = 0.1, keyheight = 0.75)) +
    theme_classic() +
    theme(aspect.ratio = 1, 
          axis.ticks.x = element_blank(), 
          panel.grid = element_blank(), 
          legend.title = element_blank(), 
          legend.position = c(0.92, 0.85))

# Plot variability for a peptide with large expected read counts
pep_large <- "pep_10160"
binom_large <- data.frame(counts = as.vector(counts(hiv_beads[pep_large, ])), 
                          n = librarySize(hiv_beads)) %>%
    rownames_to_column("sample_id") %>%
    rowwise() %>%
    mutate(cpm = counts/n*1e6,
           cpm_025 = binom.test(counts, n)$conf.int[1]*1e6, 
           cpm_975 = binom.test(counts, n)$conf.int[2]*1e6) %>%
    ungroup() %>%
    mutate(plate = str_match(sample_id, "beads_([0-9])")[, 2], 
           plate = factor(plate, levels = 1:5), 
           size = ifelse(sample_id %in% c("beads_1_D12", "beads_1_E8"), 2, 1)) %>%
    mutate(sample_id = factor(sample_id)) %>%
    ggplot(aes(x = as.numeric(sample_id))) +
    geom_linerange(aes(ymin = cpm_025, ymax = cpm_975)) +
    geom_point(aes(y = cpm, size = size, fill = plate), shape = 21) +
    scale_fill_manual(values = brewer.pal(5,"Accent"), 
                      labels = paste("Plate", 1:5)) +
    scale_y_continuous(breaks = c(40, 60, 80, 100), 
                       limits = c(30, 110)) +
    scale_x_continuous(labels = NULL, 
                       breaks =  cumsum(table(hiv_beads$plate)[-5]) + 0.5,
                       sec.axis = dup_axis(), 
                       expand = expansion(add = c(1, 1))) +
    geom_text(aes(x = xpos, y = ypos, label = label), 
              hjust = 0.5, nudge_x = 0.5, vjust = -0.3, size = 3, 
              data = data.frame(xpos = c(4, 12, 19, 24.5, 31.5), 
                                ypos = rep(110, 5), 
                                label = paste("Plate", 1:5))) +
    coord_cartesian(ylim = c(30, 110)) +
    scale_size(breaks = c(1, 2), 
               range = c(1.5, 3)) +
    labs(y = "expected reads per million", 
         x = "") +
    theme_bw() +
    theme(
        aspect.ratio = 1,
        panel.grid.major = element_line(linetype = "dotted", color = "lightgrey"),
        panel.grid.minor.y = element_line(linetype = "dotted", color = "lightgrey"), 
        panel.grid.minor.x = element_blank(), 
        axis.ticks.length.x = unit(-11, "pt"),
        legend.position = "none", 
        plot.margin = margin(0, 0, 0, 0)
    )

# Plot variability for a peptide with smaller expected read counts
pep_small <- "pep_58858"
binom_small <- data.frame(counts = as.vector(counts(hiv_beads[pep_small, ])), 
                          n = librarySize(hiv_beads)) %>%
    rownames_to_column("sample_id") %>%
    rowwise() %>%
    mutate(cpm = counts/n*1e6,
           cpm_025 = binom.test(counts, n)$conf.int[1]*1e6, 
           cpm_975 = binom.test(counts, n)$conf.int[2]*1e6) %>%
    ungroup() %>%
    mutate(plate = str_match(sample_id, "beads_([0-9])")[, 2], 
           plate = factor(plate, levels = 1:5), 
           size = ifelse(sample_id %in% c("beads_2_B5", "beads_2_D5"), 2, 1)) %>%
    mutate(sample_id = factor(sample_id)) %>%
    ggplot(aes(x = as.numeric(sample_id))) +
    geom_linerange(aes(ymin = cpm_025, ymax = cpm_975)) +
    geom_point(aes(y = cpm, size = size, fill = plate), shape = 21) +
    scale_fill_manual(values = brewer.pal(5,"Accent"), 
                      labels = paste("Plate", 1:5)) +
    scale_y_continuous(breaks = seq(0, 20, 5), 
                       limits = c(0, 20)) +
    scale_x_continuous(labels = NULL, 
                       breaks =  cumsum(table(hiv_beads$plate)[-5]) + 0.5,
                       sec.axis = dup_axis(), 
                       expand = expansion(add = c(1, 1))) +
    geom_text(aes(x = xpos, y = ypos, label = label), 
              hjust = 0.5, nudge_x = 0.5, vjust = -0.3, size = 3, 
              data = data.frame(xpos = c(4, 12, 19, 24.5, 31.5), 
                                ypos = rep(20, 5), 
                                label = paste("Plate", 1:5))) +
    coord_cartesian(ylim = c(0, 20)) +
    scale_size(breaks = c(1, 2), 
               range = c(1.5, 3)) +
    labs(y = "expected reads per million", 
         x = "") +
    theme_bw() +
    theme(
        aspect.ratio = 1,
        panel.grid.major = element_line(linetype = "dotted", color = "lightgrey"),
        panel.grid.minor.y = element_line(linetype = "dotted", color = "lightgrey"), 
        panel.grid.minor.x = element_blank(), 
        axis.ticks.length.x = unit(-11, "pt"),
        legend.position = "none", 
        plot.margin = margin(0, 0, 0, 0)
    )


ggarrange(lib_size, binom_large, binom_small, nrow = 1, 
          align = "v")

ggsave("figures/figdatastructurebinom.png", units = "in", 
       width = 12, height = 4, bg = "white")
