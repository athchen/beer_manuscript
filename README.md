
<!-- README.md is generated from README.Rmd. Please edit that file -->

# BEER Manuscript

This repository contains all of the data and code to reproduce the work
presented in, [*Detecting Antibody Reactivities in Phage
Immunoprecipitation Sequencing
Data*](https://www.biorxiv.org/content/10.1101/2022.01.19.476926v1).

The session info, including all packages used in the analyses and their
respective versions can be found in the [Session Info](#session-info)
section of this README.md. All code used to generate the figures and
tables of the manuscript are mapped in the [section
below](#mapping-of-figures-and-tables). The contents of each file in the
repository are provided in the [Directory
Structure](#directory-structure) section.

The results in the manuscript were generated using the `R` package
[`beer`](https://github.com/athchen/beer) version `0.99.0`.

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
-   **`hiv_virscan.rds`**: additional HIV data set referenced in the
    Discussion
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
-   **`run_simulations.R`** runs edgeR and BEER on simulated data sets
    with various beads-only samples and different methods of estimating
    prior parameters.

## Mapping of Figures and Tables

### Main Text

-   Figure 1 - in `figures/simulation_roc_prc_interp.png`; generated
    with `R/figure_simulation_curves.R`
-   Figure 2 - in `figures/simulation_logistic_fc.png`; generated with
    `R/figure_simulation_fc.R`
-   Figure 3 - in `figures/hiv_protein.png`; generated with
    `R/figure_hiv.R`

### Supplementary Material

-   Table S1 - generated with `R/figure_simulation_curves.R`
-   Table S2 - generated with `R/figure_simulation_curves.R`
-   Table S3 - generated with `R/figure_hiv.R`
-   Table S4 - generated with `R/figure_coronascan.R`
-   Figure S1 - in `figures/simulation_roc_prc.png`; generated with
    `R/figure_simulation_curves.R`
-   Figure S2 - in `figures/simulation_roc_bybeads.png`; generated with
    `R/figure_simulation_curves.R`
-   Figure S3 - in `figures/simulation_prc_bybeads.png`; generated with
    `R/figure_simulation_curves.R`
-   Figure S4 - in `figures/simulation_postprob.png`; generated with
    `R/figure_simulation_postprob.R`
-   Figure S5 - in `figures/simulation_fc.png`; generated with
    `R/figure_simulation_fc.R`
-   Figure S6 - in `figures/hiv_protein_noB.png`; generated with
    `R/figure_hiv.R`
-   Figure S7 - in `figures/hiv_protein_A.png`; generated with
    `R/figure_hiv.R`
-   Figure S8 - in `figures/hiv_replicates.png`; generated with
    `R/figure_hiv.R`
-   Figure S9 - in `figures/hiv_ranked_prob.png`; generated with
    `R/figure_hiv.R`
-   Figure S10 - in `figures/coronascan_protein.png`; generated with
    `R/figure_coronascan.R`
-   Figure S11 - in `figures/coronascan_replicates.png`; generated with
    `R/figure_coronascan.R`
-   Figure S12 - in `figures/coronascan_ranked_prob.png`; generated with
    `R/figure_coronascan.R`
-   Figure S13 - in `figures/corstructure.png`; generated with
    `R/figure_corstructure.R`
-   Figure S14 - in `figures/figdatastructurebinom.png`; generated with
    `R/figure_corstructure.R`
-   Figure S15 - in `figures/attnconstant.png`; generated with
    `R/figure_attnconstant.R`
-   Figure 16 - in `figures/prior.pdf`; generated with
    `R/figure_priors.R`

# Session Info

``` r
devtools::session_info()
#> ─ Session info ───────────────────────────────────────────────────────────────
#>  setting  value
#>  version  R version 4.1.2 (2021-11-01)
#>  os       macOS Big Sur 10.16
#>  system   x86_64, darwin17.0
#>  ui       X11
#>  language (EN)
#>  collate  en_US.UTF-8
#>  ctype    en_US.UTF-8
#>  tz       America/New_York
#>  date     2022-01-23
#>  pandoc   2.16.2 @ /Applications/RStudio.app/Contents/MacOS/quarto/bin/ (via rmarkdown)
#> 
#> ─ Packages ───────────────────────────────────────────────────────────────────
#>  package     * version date (UTC) lib source
#>  brio          1.1.3   2021-11-30 [1] CRAN (R 4.1.0)
#>  cachem        1.0.6   2021-08-19 [1] CRAN (R 4.1.0)
#>  callr         3.7.0   2021-04-20 [1] CRAN (R 4.1.0)
#>  cli           3.1.1   2022-01-20 [1] CRAN (R 4.1.2)
#>  crayon        1.4.2   2021-10-29 [1] CRAN (R 4.1.0)
#>  desc          1.4.0   2021-09-28 [1] CRAN (R 4.1.0)
#>  devtools      2.4.3   2021-11-30 [1] CRAN (R 4.1.0)
#>  digest        0.6.29  2021-12-01 [1] CRAN (R 4.1.0)
#>  ellipsis      0.3.2   2021-04-29 [1] CRAN (R 4.1.0)
#>  evaluate      0.14    2019-05-28 [1] CRAN (R 4.1.0)
#>  fastmap       1.1.0   2021-01-25 [1] CRAN (R 4.1.0)
#>  fs            1.5.2   2021-12-08 [1] CRAN (R 4.1.0)
#>  getPass       0.2-2   2017-07-21 [1] CRAN (R 4.1.0)
#>  git2r         0.29.0  2021-11-22 [1] CRAN (R 4.1.0)
#>  glue          1.6.0   2021-12-17 [1] CRAN (R 4.1.0)
#>  htmltools     0.5.2   2021-08-25 [1] CRAN (R 4.1.0)
#>  httpuv        1.6.5   2022-01-05 [1] CRAN (R 4.1.2)
#>  httr          1.4.2   2020-07-20 [1] CRAN (R 4.1.0)
#>  knitr         1.37    2021-12-16 [1] CRAN (R 4.1.0)
#>  later         1.3.0   2021-08-18 [1] CRAN (R 4.1.0)
#>  lifecycle     1.0.1   2021-09-24 [1] CRAN (R 4.1.0)
#>  magrittr      2.0.1   2020-11-17 [1] CRAN (R 4.1.0)
#>  memoise       2.0.1   2021-11-26 [1] CRAN (R 4.1.0)
#>  pkgbuild      1.3.1   2021-12-20 [1] CRAN (R 4.1.0)
#>  pkgload       1.2.4   2021-11-30 [1] CRAN (R 4.1.0)
#>  prettyunits   1.1.1   2020-01-24 [1] CRAN (R 4.1.0)
#>  processx      3.5.2   2021-04-30 [1] CRAN (R 4.1.0)
#>  promises      1.2.0.1 2021-02-11 [1] CRAN (R 4.1.0)
#>  ps            1.6.0   2021-02-28 [1] CRAN (R 4.1.0)
#>  purrr         0.3.4   2020-04-17 [1] CRAN (R 4.1.0)
#>  R6            2.5.1   2021-08-19 [1] CRAN (R 4.1.0)
#>  Rcpp          1.0.8   2022-01-13 [1] CRAN (R 4.1.2)
#>  remotes       2.4.2   2021-11-30 [1] CRAN (R 4.1.0)
#>  rlang         0.4.12  2021-10-18 [1] CRAN (R 4.1.0)
#>  rmarkdown     2.11    2021-09-14 [1] CRAN (R 4.1.0)
#>  rprojroot     2.0.2   2020-11-15 [1] CRAN (R 4.1.0)
#>  rstudioapi    0.13    2020-11-12 [1] CRAN (R 4.1.0)
#>  sessioninfo   1.2.2   2021-12-06 [1] CRAN (R 4.1.0)
#>  stringi       1.7.6   2021-11-29 [1] CRAN (R 4.1.0)
#>  stringr       1.4.0   2019-02-10 [1] CRAN (R 4.1.0)
#>  testthat      3.1.2   2022-01-20 [1] CRAN (R 4.1.2)
#>  usethis       2.1.5   2021-12-09 [1] CRAN (R 4.1.0)
#>  whisker       0.4     2019-08-28 [1] CRAN (R 4.1.0)
#>  withr         2.4.3   2021-11-30 [1] CRAN (R 4.1.0)
#>  workflowr   * 1.7.0   2021-12-21 [1] CRAN (R 4.1.0)
#>  xfun          0.29    2021-12-14 [1] CRAN (R 4.1.0)
#>  yaml          2.2.1   2020-02-01 [1] CRAN (R 4.1.0)
#> 
#>  [1] /Library/Frameworks/R.framework/Versions/4.1/Resources/library
#> 
#> ──────────────────────────────────────────────────────────────────────────────
```
