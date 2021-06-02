#' Code to load required packages to reproduce the results and figures in 
#' the manuscript. 
required_packages <- c('tidyverse', 'ggpubr', 'gridExtra', 'latex2exp', 
                       'kableExtra', 'RColorBrewer')
for (pkg in required_packages) {
    if (!(pkg %in% rownames(installed.packages()))) {
        install.packages(pkg)
    }
    library(pkg, character.only = TRUE)
}

# required_bioc_packages <- c("edgeR", "beer")
library(edgeR)
library(beer)

rm(list = c("required_packages", "pkg"))

#' Define global variables for plotting
hot_cold_cols <- c("navy", "blue", "deepskyblue", "cyan", "lightcyan",
                   "yellow", "orange", "red")