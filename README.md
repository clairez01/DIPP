# Discounted Individuals Power Prior (DIPP) Project

## Overview

This repository contains the programs used to conduct the simulation studies and data analysis for the paper, "Leveraging External Data in Rare Disease Trials Through Individualized Discounting Within the Power Prior: A Case Study in Hereditary Angioedema." >
The project evaluates a newly developed Bayesian methodology—DIPP—by comparing it against established approaches, including LEAP, the standard power prior, and the propensity-score integrated power prior (PSIPP). The primary goal of this study is to demonstrate that DIPP provides a robust framework for individualized data discounting, offering a more nuanced and valid alternative to traditional methods that rely on blanket discounting.

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
