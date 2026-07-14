# Process‑aware mixture modelling of precipitation-event distributions using event‑type clustering - Supplementary code

In this repository, we have collected all code and scripts necessary to (re-)run the analyses from "Process-aware mixture modelling of precipitation-event distributions using event-type clustering" by Ortega Menjivar et al. (2026+).

## Package organization

This repository is organized to "work like an `R` package", that is to say we follow all necessary conventions to simulate it as an `R` package via `devtools::load_all()`:

 - We make use of a dummy DESCRIPTION file,
 - All necessary functions are stored in the `R/` folder,
 - All analysis scripts are placed in a subfolder of `inst/` (in our case we have opted for `inst/scripts/`), where they are named and numbered by their specific step in the workflow.
- Running the analyses will populate the `data/` folder. Note that `data/` is git-ignored except for the provided aggregate lightning counts in `data/flash_events/`, so generated files will not appear in `git status`.

Note that function documentations are written in `roxygen2` syntax. This means that you can create minimalistic help pages (stored in a `man/` folder) to them by running `roxygen2::roxygenise()` or `devtools::document()` from within this repository.

## How to use this repository

 1. Clone, or download and unzip.
 2. Apply for access to the lightning data, and store it in a `data/ALDIS/` folder within this repository. For more information on this process, see [ALDIS/BLIDS](https://www.aldis.at/en/). (In order to allow for reproduction of the paper results in a review process, without the full application process and following data access, we have provided the aggregated and processed lightning counts for our selected stations in `data/flash_events/`.)
 3. Simulate loading this repo as a package by running `devtools::load_all()`. We have added comments listing which repo functions are used in each script, in case a more base `R` approach via `source()` is preferred.
 4. Run through the analysis scripts in `inst/scripts` in the order of their numbering, starting from `001_download_station_data.R` up to `006_aggregate_distribution_calc.R`. **Note that a few select analysis steps are computationally heavy.** These steps are marked as such by setting them in a `DON'T RUN` section. Rerun these steps either by decreasing relevant parameters, such as number of repetitions, or run them only on an adequate computational setup with sufficient cores and memory available.
 5. The robustness analyses for different Minimum Interevent Time (MIT) values (3, 6, 12, 18, 24 h) are orchestrated by the master script `100_robustnessAnalysis_iet_values.R`, which re-runs the full pipeline per MIT value and renders a diagnostic report per value from the template `101_analysis_report_template.Rmd` (output written to `inst/reports/`). By default, reports use the ks chosen for the 4 h baseline as starting values; after inspecting the ARI plots in a report, rerun `step6_generate_report()` with visually chosen ks.
 
## How to cite

If you use this code, please cite this repository via its archived version: [Zenodo DOI to be added upon release],  
or the accompanying paper:

> Ortega Menjivar, L., Özcelik, N. B., Laimighofer, J., Fischer, S., and Laaha, G. (2026+): Process-aware mixture modelling of precipitation-event distributions using event-type clustering, Hydrol. Earth Syst. Sci. [DOI to be added upon publication]

This code is released under the GPL-3 license (see `LICENSE`).