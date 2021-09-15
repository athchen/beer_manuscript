
<!-- README.md is generated from README.Rmd. Please edit that file -->

# BEER Manuscript

This repository contains all of the data and code to reproduce the work
presented in, *Detecting Enriched Antibody Peptides in
Phage-Immunoprecipitation Sequencing Data*.

The session info, including all packages used in the analyses and their
respective versions can be found in the [Session Info](#session-info)
section of this README.md. All code used to generate the results of the
manuscript are mapped in the [section
below](#mapping-of-manuscript-content-to-repository-files). A detailed
description of the contents of each file in the directory is provided in
the [Directory Structure](#directory-structure) section.

In particular, these results were generated using the `R` package
[`beer`](https://github.com/athchen/beer) version `0.99.0`.

# TO DO

-   In `README.Rmd`:
    -   Change `beer` hyperlink to Bioc URL once accepted
    -   Add figure and table numbers
    -   Add hyperlink to paper once it’s been uploaded/submitted.
-   In `workflowr` page:
    -   Update `beer` hyperlink to Bioc URL once accepted
    -   Review of figure/table legends and captions on each page
    -   Update CS data description
    -   Update HIV data description

# Directory structure

### **`analysis`**

Contains `.Rmd` files that are used to build
<https://athchen.github.io/beer_manuscript/>.

### **`data_processed`**

Generated output used in later parts of the analyses.

-   **`coronascan_results.rds`**: PhIPData object with edgeR and BEER
    results for the CoronaScan data set.
-   **`hiv_results.rds`**: PhIPData object with edgeR and BEER results
    for the HIV EC data set.
-   **`simulation_xbeads_zzzzz`**: folders for the simulation results
    with `x` beads-only samples and `zzzzz` method of estimating prior
    parameters for non-enriched peptides in BEER. Each folder has 10 RDS
    files, each corresponding to one simulated data set.
-   **`simulation_curves.rda`**: objects generated by `R/load_curves.R`.

### **`data_raw`**

Contains `PhIPData` objects of normalized read counts and relevant
sample and peptide annotation for example data sets and code to generate
simulated data sets.

-   **`coronascan.rds`**: coronascan data set.
-   **`hiv.rds`**: HIV elite controller data set.
-   **`simulations.R`**: script to generate simulated data sets.

### **`docs`**

`.html` files that are used to build
<https://athchen.github.io/beer_manuscript/>. `analysis/xxx.Rmd`
generates `docs/xxx.html`.

### **`figures`**

Contains png/pdfs of the manuscript figures.

### **`R`**

R scripts used to generate the objects in `data_processed` and the
figure files in `figures`.

-   **`figure_coronascan.R`** generates figures:
    -   `coronascan_protein.png`
    -   `coronascan_replicates.png`
    -   `coronascan_ranked_prob.png`
-   **`figure_hiv.R`** generates figures:
    -   `hiv_protein.png`
    -   `hiv_protein_noB.png`
    -   `hiv_protein_A.png`
    -   `hiv_replicates.png`
    -   `hiv_ranked_prob.png`
-   **`figure_simulation_curves.R`** generates tables
    `simulation_auc_roc_interp` and `simulation_ppv_by_sens` and
    figures:
    -   `simulation_roc_prc_interp.png`
    -   `simulation_roc_prc_interp_all.png`
    -   `simulation_roc_prc.png`
    -   `simulation_roc_bybeads.png`
    -   `simulation_prc_bybeads.png`
-   **`figure_simulation_enrichments.R`** generates figures:
    -   `simulation_enrichments_e.png`
    -   `simulation_enrichments_ne.png`
-   **`figure_simulation_fc.R`** generates figure `simulation_fc.png`.
-   **`figure_simulation_logistic.R`** generates figures:
    -   `simulation_logistic_fc.png`
    -   `simulation_logistic_rc.png`
-   **`figure_simulation_postprob.R`** generates figure
    `simulation_postprob.png`.
-   **`helper_functions.R`** contains helper functions for processing
    and analyzing output.
-   **`load_curves.R`** contains code to process and calculate
    sensitivity, specificity, and ppv from the simulation results.
-   **`load_packages.R`** code to load required packages and define
    global color variables.
-   **`load_simulations.R`** contains code to read in and convert
    simulation results to a tidy dataframe.
-   **`run_coronascan.R`** runs edgeR and BEER on CoronaScan data.
-   **`run_hiv.R`** runs edgeR and BEER on HIV EC data.
-   **`run_simulations.R`** run edgeR and BEER on simulated data sets
    with various beads-only samples and different methods of estimating
    prior parameters.

# Mapping of manuscript content to repository files

## Methods

## Mapping of Figures and Tables

### Main Text

### Supplementary Material

# Session Info

``` r
devtools::session_info()
#> ─ Session info ───────────────────────────────────────────────────────────────
#>  setting  value                       
#>  version  R version 4.1.1 (2021-08-10)
#>  os       macOS Big Sur 10.16         
#>  system   x86_64, darwin17.0          
#>  ui       X11                         
#>  language (EN)                        
#>  collate  en_US.UTF-8                 
#>  ctype    en_US.UTF-8                 
#>  tz       America/New_York            
#>  date     2021-09-15                  
#> 
#> ─ Packages ───────────────────────────────────────────────────────────────────
#>  package     * version date       lib source        
#>  cachem        1.0.6   2021-08-19 [1] CRAN (R 4.1.0)
#>  callr         3.7.0   2021-04-20 [1] CRAN (R 4.1.0)
#>  cli           3.0.1   2021-07-17 [1] CRAN (R 4.1.0)
#>  crayon        1.4.1   2021-02-08 [1] CRAN (R 4.1.0)
#>  desc          1.3.0   2021-03-05 [1] CRAN (R 4.1.0)
#>  devtools      2.4.2   2021-06-07 [1] CRAN (R 4.1.0)
#>  digest        0.6.27  2020-10-24 [1] CRAN (R 4.1.0)
#>  ellipsis      0.3.2   2021-04-29 [1] CRAN (R 4.1.0)
#>  evaluate      0.14    2019-05-28 [1] CRAN (R 4.1.0)
#>  fansi         0.5.0   2021-05-25 [1] CRAN (R 4.1.0)
#>  fastmap       1.1.0   2021-01-25 [1] CRAN (R 4.1.0)
#>  fs            1.5.0   2020-07-31 [1] CRAN (R 4.1.0)
#>  glue          1.4.2   2020-08-27 [1] CRAN (R 4.1.0)
#>  htmltools     0.5.2   2021-08-25 [1] CRAN (R 4.1.0)
#>  httpuv        1.6.2   2021-08-18 [1] CRAN (R 4.1.0)
#>  knitr         1.33    2021-04-24 [1] CRAN (R 4.1.0)
#>  later         1.3.0   2021-08-18 [1] CRAN (R 4.1.0)
#>  lifecycle     1.0.0   2021-02-15 [1] CRAN (R 4.1.0)
#>  magrittr      2.0.1   2020-11-17 [1] CRAN (R 4.1.0)
#>  memoise       2.0.0   2021-01-26 [1] CRAN (R 4.1.0)
#>  pillar        1.6.2   2021-07-29 [1] CRAN (R 4.1.0)
#>  pkgbuild      1.2.0   2020-12-15 [1] CRAN (R 4.1.0)
#>  pkgconfig     2.0.3   2019-09-22 [1] CRAN (R 4.1.0)
#>  pkgload       1.2.1   2021-04-06 [1] CRAN (R 4.1.0)
#>  prettyunits   1.1.1   2020-01-24 [1] CRAN (R 4.1.0)
#>  processx      3.5.2   2021-04-30 [1] CRAN (R 4.1.0)
#>  promises      1.2.0.1 2021-02-11 [1] CRAN (R 4.1.0)
#>  ps            1.6.0   2021-02-28 [1] CRAN (R 4.1.0)
#>  purrr         0.3.4   2020-04-17 [1] CRAN (R 4.1.0)
#>  R6            2.5.1   2021-08-19 [1] CRAN (R 4.1.0)
#>  Rcpp          1.0.7   2021-07-07 [1] CRAN (R 4.1.0)
#>  remotes       2.4.0   2021-06-02 [1] CRAN (R 4.1.0)
#>  rlang         0.4.11  2021-04-30 [1] CRAN (R 4.1.0)
#>  rmarkdown     2.10    2021-08-06 [1] CRAN (R 4.1.0)
#>  rprojroot     2.0.2   2020-11-15 [1] CRAN (R 4.1.0)
#>  rstudioapi    0.13    2020-11-12 [1] CRAN (R 4.1.0)
#>  sessioninfo   1.1.1   2018-11-05 [1] CRAN (R 4.1.0)
#>  stringi       1.7.4   2021-08-25 [1] CRAN (R 4.1.0)
#>  stringr       1.4.0   2019-02-10 [1] CRAN (R 4.1.0)
#>  testthat      3.0.4   2021-07-01 [1] CRAN (R 4.1.0)
#>  tibble        3.1.4   2021-08-25 [1] CRAN (R 4.1.0)
#>  usethis       2.0.1   2021-02-10 [1] CRAN (R 4.1.0)
#>  utf8          1.2.2   2021-07-24 [1] CRAN (R 4.1.0)
#>  vctrs         0.3.8   2021-04-29 [1] CRAN (R 4.1.0)
#>  withr         2.4.2   2021-04-18 [1] CRAN (R 4.1.0)
#>  workflowr   * 1.6.2   2020-04-30 [1] CRAN (R 4.1.0)
#>  xfun          0.25    2021-08-06 [1] CRAN (R 4.1.0)
#>  yaml          2.2.1   2020-02-01 [1] CRAN (R 4.1.0)
#> 
#> [1] /Library/Frameworks/R.framework/Versions/4.1/Resources/library
```
