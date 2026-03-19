# Discounting Individuals Power Prior (DIPP) Project

## Overview

This repository contains the programs used to conduct the simulation studies and data analysis for the paper, "Leveraging External Data in Rare Disease Trials Through Individualized Discounting Within the Power Prior: A Case Study in Hereditary Angioedema."
The project evaluates a newly developed Bayesian methodology (DIPP) by comparing it against established approaches, including LEAP, the standard power prior, the normalized power prior, and the propensity-score integrated power prior (PSIPP). The primary goal of this study is to demonstrate that DIPP provides a robust framework for individualized data discounting, offering a more nuanced and valid alternative to traditional methods that rely on blanket discounting.

## Repository Structure

Below is an outline of the directory structure and a brief description of the contents within each folder:

```text
├── Data_Analysis/
│   ├── COMPACT Data Analysis.R
│   ├── LEAP data analysis.R
│   ├── partial dependence plots.R
│   └── propensity score analysis.R
├── DummyData/
│   ├── curr_dummy.rds
│   └── hist_dummy.rds
└── Simulation_Study/
    ├── Batch_Jobs/
    │   ├── dipp_shell.sh
    │   ├── leap_shell.sh
    │   ├── ndipp_shell.sh
    │   ├── pp_shell.sh
    │   ├── psipp_shell.sh
    │   └── sims_add_shell.sh
    ├── Compile_Results/
    │   ├── compile_additivemodels_results.R
    │   ├── compile_dipp_sim_results.R
    │   ├── compile_leap_sim_results.R
    │   ├── compile_ndipp_sim_results.R
    │   ├── compile_pp_sim_results.R
    │   └── compile_psipp_sim_results.R
    ├── Data_Generation/
    │   ├── create grid.R
    │   ├── generate data additive.R
    │   └── generate data mult.R
    └── Programs/
        ├── leap_sim.R
        ├── nimble_dipp_simcode.R
        ├── nimble_ndipp_simcode.R
        ├── nimble_powerprior_simcode.R
        ├── nimble_psipp_simcode.R
        └── sims_additivenonexch.R
```

## Dependencies and Reproducibility

The analysis was conducted using **R version 4.4.0** on a Red Hat Enterprise Linux system. For full transparency and version control of all packages used (including `nimble`, `tidyverse`, and `cmdstanr`), please see the session information below.

<details>
<summary>Click to expand R Session Info</summary>

```text
 R version 4.4.0 (2024-04-24)
Platform: x86_64-pc-linux-gnu
Running under: Red Hat Enterprise Linux 9.7 (Plow)

Matrix products: default
BLAS/LAPACK: /nas/longleaf/rhel9/apps/r/4.4.0/lib/libopenblas_haswellp-r0.3.27.so;  LAPACK version 3.12.0

locale:
 [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C               LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8     LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
 [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                  LC_ADDRESS=C               LC_TELEPHONE=C             LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       

time zone: America/New_York
tzcode source: system (glibc)

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] ggpubr_0.6.1    mvtnorm_1.3-3   sandwich_3.1-1  lmtest_0.9-40   zoo_1.8-14      MatchIt_4.7.2   patchwork_1.3.0 cmdstanr_0.9.0  hdbayes_0.1.1   haven_2.5.4     posterior_1.6.1
[12] lubridate_1.9.4 forcats_1.0.0   stringr_1.5.1   dplyr_1.1.4     purrr_1.0.4     readr_2.1.5     tidyr_1.3.1     tibble_3.2.1    ggplot2_3.5.2   tidyverse_2.0.0 nimble_1.3.0   

loaded via a namespace (and not attached):
 [1] gtable_0.3.6         tensorA_0.36.2.1     xfun_0.52            formula.tools_1.7.1  rstatix_0.7.2        processx_3.8.6       lattice_0.22-7       callr_3.7.6         
 [9] tzdb_0.5.0           numDeriv_2016.8-1.1  vctrs_0.6.5          tools_4.4.0          ps_1.9.1             generics_0.1.4       parallel_4.4.0       enrichwith_0.4.0    
[17] pkgconfig_2.0.3      Matrix_1.7-3         checkmate_2.3.2      RColorBrewer_1.1-3   distributional_0.5.0 lifecycle_1.0.4      compiler_4.4.0       farver_2.1.2        
[25] Brobdingnag_1.2-9    carData_3.0-5        Formula_1.2-5        pracma_2.4.4         car_3.1-3            pillar_1.10.2        bridgesampling_1.1-2 abind_1.4-8         
[33] mclust_6.1.1         tidyselect_1.2.1     stringi_1.8.7        operator.tools_1.6.3 grid_4.4.0           cli_3.6.5            magrittr_2.0.3       dichromat_2.0-0.1   
[41] broom_1.0.8          withr_3.0.2          scales_1.4.0         backports_1.5.0      timechange_0.3.0     igraph_2.1.4         instantiate_0.2.3    ggsignif_0.6.4      
[49] hms_1.1.3            coda_0.19-4.1        evaluate_1.0.3       knitr_1.50           rlang_1.1.6          Rcpp_1.0.14          glue_1.8.0           rstudioapi_0.17.1   
[57] jsonlite_2.0.0       R6_2.6.1             fs_1.6.6
```
