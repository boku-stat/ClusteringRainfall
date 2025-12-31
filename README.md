# Improving stochastic rainfall models through event clustering - Supplementary code

In this repository, we have collected all code and scripts necessary to (re-)run the analyses from "Improving stochastic rainfall models through event clustering" by Ortega Menjivar et al. (2026+).

## Package organization

This repository is organized to "work like an `R` package", that is to say we follow all necessary conventions to simulate it as an `R` package via `devtools::load_all()`:

    - Ee make use of a dummy DESCRIPTION file,
    - All necessary functions are stored in the `R/` folder,
    - All analysis script are placed in a subfolder of `inst/` (in our case we have opted for `inst/scripts/`), where they are named and numbered by their specific step in the workflow.
    - Running the analyses will create a `data/` folder, 
Note that function documentations are written in `roxygen2` syntax. This means that you can create minimalistic help pages (stored in a `man/` folder) to them by running `roxygen2::roxygenise()` or `devtools::document()` from within this repository.

## How to use this repository

 1. Clone, or download and unzip.
 2. Apply for access to the lightning data, and store it in a `data/` folder within this repository.
 3. Simulate loading this repo as a package by running `devtools::load_all()`. I have added statements explaining which repo functions are used in each script, in case a more base `R` approach via `source()` is preferred.
 4. Run through the analysis scripts in `inst/scripts` in the order of their numbering, starting from `001_download_station_data.R` up to `006_aggregate_distribution_calc.R`. **Note that a few select analysis steps are computationally heavy.** These steps are marked as such by setting them in a `DON'T RUN` section. Rerun this steps either by decreasing relevant parameters, such as number of repetitions, or run them only on an adequate computational setup with sufficient cores and GPU available.
 5. The robustness analyses, using different Minimum Interevent Time levels are stored in `101_robustness_MIT6.Rmd` to `103_robustness_MIT24.Rmd`.

## TODO:

 -[ ] Lena: Add note on how to get access to the data. (Will do that once I know how to do that for the lightning data)
 -[ ] Lena: Review all filepaths to new locations
 -[ ] Lena: Review whether I have correctly transferred all necessary functions from `EROSA-Stat`
 -[ ] Lena: Run robustness analyses for different MIT values (finally!), and add these scripts as 101...R - 103...R
 -[ ] Others: Please review whether your scripts and functions have valid `roxygen2` syntax, and `DON'T RUN`-statements for computationally heavy operations.
 
