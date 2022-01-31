#' figure_postpred.R
#' 
#' Code to generate figures:
#' - Figure S13: postpred.png

# Set-up --------------
if(!"here" %in% installed.packages()){
    install.packages(here)
}
source(here::here("R", "load_packages.R"))

hiv <- readRDS(here("data_processed", "hiv_results.rds"))
hiv_samples <- readRDS(here("data_processed", "hiv_samples.rds"))

thetas <- as.matrix(hiv_samples)
thetas <- thetas[, grepl("theta", colnames(thetas))]

sample_num <- which(hiv$group != "beads")[1]
n <- librarySize(hiv)[sample_num]

# Sample from the posterior predictive distribution
set.seed(1234)
postpred_samples <- apply(thetas, 2, function(chain){
    as.vector(sapply(chain, function(x){rbinom(10, n, x)}))
}) 

# Get 95% credible intervals
ci <- map_dfr(1:nrow(hiv), function(row_number){
    row <- postpred_samples[, row_number]
    obs_count <- counts(hiv)[row_number, sample_num]
    
    output <- c(obs_count, 
      quantile(row, c(0.025, 0.5, 0.975)), 
      min(sum(row <= obs_count)/length(row), sum(row > obs_count)/length(row)))
    names(output) <-  c("counts","low", "med", "upper", "p_value")
    output
})

p_hist <- ci %>%
    ggplot(aes(x = p_value)) + 
    geom_histogram(binwidth = 0.025, boundary = 0, color = "white") +
    labs(x = "posterior predictive p-value", 
         y = "# peptides") +
    theme_bw() +
    theme(aspect.ratio = 1)

p_interval <- ci %>%
    arrange(counts) %>%
    mutate(peptide = 1:n()) %>%
    filter(peptide %% 50 == 0) %>%
    ggplot(aes(y = peptide)) +
    geom_errorbarh(aes(xmin = low, xmax = upper, color = "95% interval")) +
    geom_point(aes(x = counts, color = "observed count")) + 
    geom_point(aes(x = med, color = "median"), fill = NA, shape = 21) +
    coord_cartesian(xlim = c(0, 160), ylim = c(1, 3300)) +
    scale_x_continuous(breaks = seq(0, 160, by = 20)) +
    labs(x = "read counts", y = "", color = "") +
    scale_color_manual(breaks = c("95% interval", "median", "observed count"), 
                       values = c("black", "black", "blue")) +
    guides(colour = guide_legend(
        override.aes = list(linetype = c("solid", "blank", "blank"), 
                            shape = c(NA, 21, 16)))) +
    theme_bw() + 
    theme(aspect.ratio = 1, 
          legend.title = element_blank(), 
          panel.grid.major.y = element_blank(),
          panel.grid.minor.y = element_blank(), 
          axis.ticks.y = element_blank(), 
          axis.text.y = element_blank(), 
          legend.background = element_rect(color = "black"), 
          legend.position = c(0.75, 0.17)) 

ggarrange(p_hist, p_interval, nrow = 1, align = "v")

ggsave("figures/postpred.png", width = 8, height = 4)
