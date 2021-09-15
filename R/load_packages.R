#' Code to load required packages to reproduce the results and figures in 
#' the manuscript. 
required_packages <- c('plyr', 'tidyverse', 'here', 'ggpubr', 'gridExtra', 
                       'latex2exp', 'kableExtra', 'RColorBrewer', 'BiocManager')
for (pkg in required_packages) {
    if (!(pkg %in% rownames(installed.packages()))) {
        install.packages(pkg)
    }
    library(pkg, character.only = TRUE)
}

bioc_packages <- c("beer")
for(pkg in bioc_packages){
    if(!(pkg %in% rownames(installed.packages()))) {
        BiocManager::install(pkg)
    }
    library(pkg, character.only = TRUE)
}

rm(list = c("required_packages", "bioc_packages", "pkg"))

#' Define global variables for plotting
hot_cold_cols <- c("navy", "blue", "deepskyblue", "cyan", "lightcyan",
                   "yellow", "orange", "red")
