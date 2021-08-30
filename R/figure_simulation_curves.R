#' figure_simulation_curves.R
#' 
#' Code to generate tables:
#' - simulation_auc_roc_interp
#' - simulation_ppv_by_sens
#' 
#' Code to generate figures:
#' - simulation_roc_prc_interp.png
#' - simulation_roc_prc_interp_all.png
#' - simulation_roc_prc.png
#' - simulation_roc_bybeads.png
#' - simulation_prc_bybeads.png

# Set-up --------------
if(file.exists(here("data_processed", "simulation_curves.rda"))){
    load(here("data_processed", "simulation_curves.rda"))
} else {
    # Takes a while to run - only need to run once. 
    source(here("R", "load_curves.R"))
}

curves_summary <- interpolate %>%
    group_by(group, group_lab, method, num_beads, x) %>%
    summarize(n_point = n(), 
              sens_na = sum(is.na(sens)), 
              ppv_na = sum(is.na(ppv)), 
              mean_sens = mean(sens, na.rm = TRUE), 
              mean_ppv = mean(ppv, na.rm = TRUE), 
              .groups = "drop")

# Table: simulation_auc_roc_interp ------------
curves_summary %>%
    filter(!grepl("truth", method)) %>%
    group_by(method, group, num_beads) %>%
    arrange(x, .by_group = TRUE) %>%
    summarize(auc_roc = integrate_vector(x, mean_sens),
              .groups = "drop") %>%
    pivot_wider(names_from = "num_beads", values_from = "auc_roc") %>%
    arrange(method, group)

# Table: ppv_by_sens ----------
curves_summary %>%
    filter(round(x, 2) %in% c(0.5, 0.75, 0.9, 0.95) &
               grepl("edgeR", method)) %>%
    select(method, num_beads, group, x, mean_ppv) %>%
    mutate(mean_ppv = round(mean_ppv, 3)) %>%
    pivot_wider(names_from = x, values_from = mean_ppv) %>%
    group_by(num_beads) %>%
    arrange(method, group) %>%
    group_split() 

# Figure: simulation_roc_prc_interp.png ------------
roc_interp <- curves_summary %>%
    filter(method %in% c("BEER_truth", "BEER_edgeR", "edgeR")) %>% 
    ggplot(aes(x = x, y = mean_sens, group = method)) +
    geom_line(aes(color = method), size = 0.5) +
    facet_grid(num_beads ~ group_lab, 
               labeller = labeller(group_lab = label_parsed, 
                                   num_beads = c("2" = "2 beads",
                                                 "4" = "4 beads", 
                                                 "8" = "8 beads"))) +
    coord_fixed(xlim = c(0, 1), ylim = c(0, 1)) +
    labs(x = "1 - specificity", y = "sensitivity",
         color = "method") +
    scale_color_manual(breaks = c("BEER_truth", "BEER_mom", "BEER_mle", 
                                  "BEER_edgeR", "edgeR"),
                       labels = c("BEER, truth", "BEER, MOM", 
                                  "BEER, MLE", "BEER, edgeR", 
                                  "edgeR"), 
                       values = c(brewer.pal(5, "Reds")[2:5],
                                  brewer.pal(3, "Greys")[c(3)])) +
    scale_y_continuous(breaks = seq(0, 1, 0.2)) +
    scale_x_continuous(breaks = seq(0, 1, 0.2)) +
    theme_bw() +
    theme(legend.title = element_text(size = 6), 
          legend.text = element_text(size = 6), 
          legend.key.size = unit(0.6, "lines"), 
          legend.background = element_rect(color = "black", size = 0.3),  
          legend.position = c(0.92, 0.78), 
          axis.title = element_text(size = 8), 
          axis.text = element_text(size = 8), 
          plot.margin = unit(c(0, 0, 0.25, 0), "cm"))

prc_interp <- curves_summary %>%
    filter(method %in% c("BEER_truth", "BEER_edgeR", "edgeR")) %>% 
    ggplot(aes(x = x, y = mean_ppv, group = method)) +
    geom_line(aes(color = method), size = 0.5) +
    facet_grid(num_beads ~ group_lab, 
               labeller = labeller(group_lab = label_parsed, 
                                   num_beads = c("2" = "2 beads",
                                                 "4" = "4 beads", 
                                                 "8" = "8 beads"))) +
    coord_fixed(xlim = c(0, 1), ylim = c(0, 1)) +
    labs(x = "sensitivity", y = "positive predictive value (PPV)",
         color = "method") +
    scale_color_manual(breaks = c("BEER_truth", "BEER_mom", "BEER_mle", 
                                  "BEER_edgeR", "edgeR"),
                       labels = c("BEER, truth", "BEER, MOM", 
                                  "BEER, MLE", "BEER, edgeR", 
                                  "edgeR"), 
                       values = c(brewer.pal(5, "Reds")[2:5],
                                  brewer.pal(3, "Greys")[c(3)])) +
    scale_y_continuous(breaks = seq(0, 1, 0.2)) +
    scale_x_continuous(breaks = seq(0, 1, 0.2)) +
    theme_bw() +
    theme(axis.text = element_text(size = 8), 
          axis.title = element_text(size = 8), 
          plot.margin = unit(c(0, 0, 0.25, 0), "cm"))

roc_prc_interp <- ggarrange(roc_interp, 
                            prc_interp + theme(legend.position = "none"), 
                            ncol = 1, nrow = 2)

ggsave("figures/simulation_roc_prc_interp.png", roc_prc_interp, 
       units = "in", height = 9, width = 7)

# Figure: simulation_roc_prc_interp_all.png ------------
roc_interp <- curves_summary %>%
    ggplot(aes(x = x, y = mean_sens, group = method)) +
    geom_line(aes(color = method), size = 0.5) +
    facet_grid(num_beads ~ group_lab, 
               labeller = labeller(group_lab = label_parsed, 
                                   num_beads = c("2" = "2 beads",
                                                 "4" = "4 beads", 
                                                 "8" = "8 beads"))) +
    coord_fixed(xlim = c(0, 1), ylim = c(0, 1)) +
    labs(x = "1 - specificity", y = "sensitivity",
         color = "method") +
    scale_color_manual(breaks = c("BEER_truth", "BEER_mom", "BEER_mle", 
                                  "BEER_edgeR", "edgeR"),
                       labels = c("BEER, truth", "BEER, MOM", 
                                  "BEER, MLE", "BEER, edgeR", 
                                  "edgeR"), 
                       values = c(brewer.pal(5, "Reds")[2:5],
                                  brewer.pal(3, "Greys")[c(3)])) +
    scale_y_continuous(breaks = seq(0, 1, 0.2)) +
    scale_x_continuous(breaks = seq(0, 1, 0.2)) +
    theme_bw() +
    theme(legend.title = element_text(size = 6), 
          legend.text = element_text(size = 6), 
          legend.key.size = unit(0.6, "lines"), 
          legend.background = element_rect(color = "black", size = 0.3),  
          legend.position = c(0.92, 0.81),
          axis.title = element_text(size = 8), 
          axis.text = element_text(size = 8), 
          plot.margin = unit(c(0, 0, 0.25, 0), "cm"))

prc_interp <- curves_summary %>%
    ggplot(aes(x = x, y = mean_ppv, group = method)) +
    geom_line(aes(color = method), size = 0.5) +
    facet_grid(num_beads ~ group_lab, 
               labeller = labeller(group_lab = label_parsed, 
                                   num_beads = c("2" = "2 beads",
                                                 "4" = "4 beads", 
                                                 "8" = "8 beads"))) +
    coord_fixed(xlim = c(0, 1), ylim = c(0, 1)) +
    labs(x = "sensitivity", y = "positive predictive value (PPV)",
         color = "method") +
    scale_color_manual(breaks = c("BEER_truth", "BEER_mom", "BEER_mle", 
                                  "BEER_edgeR", "edgeR"),
                       labels = c("BEER, truth", "BEER, MOM", 
                                  "BEER, MLE", "BEER, edgeR", 
                                  "edgeR"), 
                       values = c(brewer.pal(5, "Reds")[2:5],
                                  brewer.pal(3, "Greys")[c(3)])) +
    scale_y_continuous(breaks = seq(0, 1, 0.2)) +
    scale_x_continuous(breaks = seq(0, 1, 0.2)) +
    theme_bw() +
    theme(axis.text = element_text(size = 8), 
          axis.title = element_text(size = 8), 
          plot.margin = unit(c(0, 0, 0.25, 0), "cm"))

roc_prc_interp <- ggarrange(roc_interp, 
                            prc_interp + theme(legend.position = "none"), 
                            ncol = 1, nrow = 2)

ggsave("figures/simulation_roc_prc_interp_all.png", roc_prc_interp, 
       units = "in", height = 9, width = 7)

# Figure: simulation_roc_prc.png ------------
roc <- roc_by_fc %>%
    ggplot(aes(x = 1-spec, y = sens, group = paste0(sim_num, method))) +
    geom_line(aes(color = method), size = 0.2) +
    facet_grid(num_beads ~ group_lab, 
               labeller = labeller(group_lab = label_parsed, 
                                   num_beads = c("2" = "2 beads",
                                                 "4" = "4 beads", 
                                                 "8" = "8 beads"))) +
    coord_fixed(xlim = c(0, 1), ylim = c(0, 1)) +
    labs(x = "1 - specificity", y = "sensitivity", color = "method") +
    scale_color_manual(breaks = c("BEER_truth", "BEER_mom", "BEER_mle", 
                                  "BEER_edgeR", "edgeR"),
                       labels = c("BEER, truth", "BEER, MOM", 
                                  "BEER, MLE", "BEER, edgeR", 
                                  "edgeR"), 
                       values = c(brewer.pal(5, "Reds")[2:5],
                                  brewer.pal(3, "Greys")[c(3)])) +
    scale_y_continuous(breaks = seq(0, 1, 0.2)) +
    scale_x_continuous(breaks = seq(0, 1, 0.2)) +
    theme_bw() +
    theme(legend.title = element_text(size = 6), 
          legend.text = element_text(size = 6), 
          legend.key.size = unit(0.6, "lines"), 
          legend.background = element_rect(color = "black", size = 0.3),  
          legend.position = c(0.93, 0.81), 
          axis.title = element_text(size = 8), 
          axis.text = element_text(size = 8), 
          plot.margin = unit(c(0, 0, 0.25, 0), "cm"))

prc <- roc_by_fc %>%
    ggplot(aes(x = sens, y = ppv, group = paste0(sim_num, method))) +
    geom_line(aes(color = method), size = 0.2) +
    facet_grid(num_beads ~ group_lab, 
               labeller = labeller(group_lab = label_parsed, 
                                   num_beads = c("2" = "2 beads",
                                                 "4" = "4 beads", 
                                                 "8" = "8 beads"))) +
    coord_fixed(xlim = c(0, 1), ylim = c(0, 1)) +
    labs(x = "sensitivity", y = "PPV", color = "method") +
    scale_color_manual(breaks = c("BEER_truth", "BEER_mom", "BEER_mle", 
                                  "BEER_edgeR", "edgeR"),
                       labels = c("BEER, truth", "BEER, MOM", 
                                  "BEER, MLE", "BEER, edgeR", 
                                  "edgeR"), 
                       values = c(brewer.pal(5, "Reds")[2:5],
                                  brewer.pal(3, "Greys")[c(3)])) +
    scale_y_continuous(breaks = seq(0, 1, 0.2)) +
    scale_x_continuous(breaks = seq(0, 1, 0.2)) +
    theme_bw() +
    theme(axis.text = element_text(size = 8), 
          axis.title = element_text(size = 8), 
          plot.margin = unit(c(0, 0, 0.25, 0), "cm"))

roc_prc <- ggarrange(roc, prc + theme(legend.position = "none"), ncol = 1, nrow = 2)

ggsave("figures/simulation_roc_prc.png", roc_prc, 
       units = "in", height = 9, width = 8.5)

# Figure: simulation_roc_bybeads.png ------------
curves_summary %>%
    filter(method !="BEER_truth") %>% 
    ggplot(aes(x = x, y = mean_sens, group = num_beads)) +
    geom_line(aes(color = factor(num_beads)), size = 0.5) +
    facet_grid(method ~ group_lab, 
               labeller = labeller(group_lab = label_parsed, 
                                   method = c("BEER_truth" = "BEER, truth",
                                              "BEER_mom" = "BEER, MOM", 
                                              "BEER_mle" = "BEER, MLE", 
                                              "BEER_edgeR" = "BEER, edgeR", 
                                              "edgeR" = "edgeR"))) +
    coord_fixed(xlim = c(0, 1), ylim = c(0, 1)) +
    labs(x = "1 - specificity", y = "sensitivity",
         color = "# beads") +
    scale_color_manual(breaks = c("2", "4", "8"),
                       values = c("red", "blue", "black")) +
    scale_y_continuous(breaks = seq(0, 1, 0.2)) +
    scale_x_continuous(breaks = seq(0, 1, 0.2)) +
    theme_bw() +
    theme(legend.title = element_text(size = 6), 
          legend.text = element_text(size = 6), 
          legend.key.size = unit(0.6, "lines"), 
          legend.background = element_rect(color = "black", size = 0.3),  
          legend.position = c(0.94, 0.08), 
          axis.title = element_text(size = 8), 
          axis.text = element_text(size = 8), 
          plot.margin = unit(c(0, 0, 0.25, 0), "cm"))

ggsave("figures/simulation_roc_bybeads.png", units = "in", 
       height = 5, width = 6)

# Figure: simulation_prc_bybeads.png ------------
curves_summary %>%
    filter(method !="BEER_truth") %>% 
    ggplot(aes(x = x, y = mean_ppv, group = num_beads)) +
    geom_line(aes(color = factor(num_beads)), size = 0.5) +
    facet_grid(method ~ group_lab, 
               labeller = labeller(group_lab = label_parsed)) +
    coord_fixed(xlim = c(0, 1), ylim = c(0, 1)) +
    labs(x = "sensitivity", y = "positive predictive value (PPV)",
         color = "# beads") +
    scale_color_manual(breaks = c("2", "4", "8"),
                       values = c("red", "blue", "black")) +
    scale_y_continuous(breaks = seq(0, 1, 0.2)) +
    scale_x_continuous(breaks = seq(0, 1, 0.2)) +
    theme_bw() +
    theme(legend.title = element_text(size = 6), 
          legend.text = element_text(size = 6), 
          legend.key.size = unit(0.6, "lines"), 
          legend.background = element_rect(color = "black", size = 0.3),  
          legend.position = c(0.94, 0.08), 
          axis.title = element_text(size = 8), 
          axis.text = element_text(size = 8), 
          plot.margin = unit(c(0, 0, 0.25, 0), "cm"))

ggsave("figures/simulation_prc_bybeads.png", units = "in", 
       height = 5, width = 6)
