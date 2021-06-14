#' figure_simulation_curves.R
#' 
#' Code to generate figures:
#' - simulation_roc_prc_interp.png
#' - simulation_roc_prc.png

# Set-up --------------
source(file.path("R", "load_curves.R"))

# Figure: simulation_roc_prc_interp.png ------------
roc_interp <- interpolate %>%
    group_by(group, group_lab, method, num_beads, x) %>%
    summarize(n_point = n(), 
              sens_na = sum(is.na(sens)), 
              ppv_na = sum(is.na(ppv)), 
              mean_sens = mean(sens, na.rm = TRUE), 
              mean_ppv = mean(ppv, na.rm = TRUE)) %>%
    filter(method %in% c("BEER_truth", "BEER_edgeR", "edgeR")) %>% 
    ggplot(aes(x = x, y = mean_sens, group = method)) +
    geom_line(aes(color = method), size = 0.5) +
    facet_grid(num_beads ~ group_lab, 
               labeller = labeller(group_lab = label_parsed, 
                                   num_beads = c("4" = "4 beads", 
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
          legend.position = c(0.92, 0.675), 
          axis.title = element_text(size = 8), 
          axis.text = element_text(size = 8), 
          plot.margin = unit(c(0, 0, 0.25, 0), "cm"))

prc_interp <- interpolate %>%
    group_by(group, group_lab, method, num_beads, x) %>%
    summarize(n_point = n(), 
              sens_na = sum(is.na(sens)), 
              ppv_na = sum(is.na(ppv)), 
              mean_sens = mean(sens, na.rm = TRUE), 
              mean_ppv = mean(ppv, na.rm = TRUE), 
              .groups = "drop") %>%
    filter(method %in% c("BEER_truth", "BEER_edgeR", "edgeR")) %>% 
    ggplot(aes(x = x, y = mean_ppv, group = method)) +
    geom_line(aes(color = method), size = 0.5) +
    facet_grid(num_beads ~ group_lab, 
               labeller = labeller(group_lab = label_parsed, 
                                   num_beads = c("4" = "4 beads", 
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
       units = "in", height = 6, width = 6.6)

# Figure: simulation_roc_prc_interp_all.png ------------
roc_interp <- interpolate %>%
    group_by(group, group_lab, method, num_beads, x) %>%
    summarize(n_point = n(), 
              sens_na = sum(is.na(sens)), 
              ppv_na = sum(is.na(ppv)), 
              mean_sens = mean(sens, na.rm = TRUE), 
              mean_ppv = mean(ppv, na.rm = TRUE)) %>%
    ggplot(aes(x = x, y = mean_sens, group = method)) +
    geom_line(aes(color = method), size = 0.5) +
    facet_grid(num_beads ~ group_lab, 
               labeller = labeller(group_lab = label_parsed, 
                                   num_beads = c("4" = "4 beads", 
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
          legend.position = c(0.92, 0.725),
          axis.title = element_text(size = 8), 
          axis.text = element_text(size = 8), 
          plot.margin = unit(c(0, 0, 0.25, 0), "cm"))

prc_interp <- interpolate %>%
    group_by(group, group_lab, method, num_beads, x) %>%
    summarize(n_point = n(), 
              sens_na = sum(is.na(sens)), 
              ppv_na = sum(is.na(ppv)), 
              mean_sens = mean(sens, na.rm = TRUE), 
              mean_ppv = mean(ppv, na.rm = TRUE), 
              .groups = "drop") %>%
    ggplot(aes(x = x, y = mean_ppv, group = method)) +
    geom_line(aes(color = method), size = 0.5) +
    facet_grid(num_beads ~ group_lab, 
               labeller = labeller(group_lab = label_parsed, 
                                   num_beads = c("4" = "4 beads", 
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
       units = "in", height = 6, width = 6.6)

# Figure: simulation_roc_prc.png ------------
roc <- roc_by_fc %>%
    ggplot(aes(x = 1-spec, y = sens, group = paste0(sim_num, method))) +
    geom_line(aes(color = method), size = 0.2) +
    facet_grid(num_beads ~ group_lab, 
               labeller = labeller(group_lab = label_parsed, 
                                   num_beads = c("4" = "4 beads", 
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
          legend.position = c(0.93, 0.74), 
          axis.title = element_text(size = 8), 
          axis.text = element_text(size = 8), 
          plot.margin = unit(c(0, 0, 0.25, 0), "cm"))

prc <- roc_by_fc %>%
    ggplot(aes(x = sens, y = ppv, group = paste0(sim_num, method))) +
    geom_line(aes(color = method), size = 0.2) +
    facet_grid(num_beads ~ group_lab, 
               labeller = labeller(group_lab = label_parsed, 
                                   num_beads = c("4" = "4 beads", 
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
       units = "in", height = 6)
