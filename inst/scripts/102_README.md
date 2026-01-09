# Precipitation Event Analysis Pipeline

This repository contains a master R script and report template for analyzing precipitation events with different Inter-Event Time (IET) thresholds.

## Overview

The pipeline performs:
1. **Data Preparation**: Combines flash, precipitation, and event data
2. **Partitioning Clustering**: Identifies event clusters with stability analysis
3. **Clusterwise Regression**: Fits Gamma GLM models to cluster-specific precipitation patterns
4. **Report Generation**: Creates detailed HTML reports with visualizations

## Files

### `master_analysis.R`
Main orchestration script that runs the entire analysis pipeline for multiple IET values.

**Key Features:**
- Processes IET values: 3, 6, 12, 18, and 24 hours
- Saves all models to `data/` folder
- Generates reports for each IET value
- Configurable parameters for stations, clustering, and regression

### `analysis_report_template.Rmd`
R Markdown template that generates comprehensive reports including:
- Adjusted Rand Index (ARI) visualizations
- Kullback-Leibler divergence comparison tables
- Factor Analysis of Mixed Data (FAMD) plots
- Partial Dependence Plots (PDP)
- Empirical Cumulative Distribution Functions (ECDF)
- CDF comparison plots (empirical vs. fitted)

## Prerequisites

### Required R Packages

```r
# Core packages
install.packages(c("tidyverse", "sf", "flexclust", "flexmix"))

# For reporting
install.packages(c("rmarkdown", "knitr", "kableExtra"))

# For visualization
install.packages(c("FactoMineR", "factoextra", "pdp", "patchwork"))

# For parallel processing
install.packages("parallel")
```

### Custom Functions

The scripts use several custom functions that should be available in your development package:
- `get_precip()`: Extract precipitation data for time windows
- `get_flash_data()`: Check flash occurrence in time windows
- `do_stabAn()`: Perform stability analysis on clustering
- `bootEventmix()`: Bootstrap analysis for flexmix models

**IMPORTANT:** Update these function definitions in `master_analysis.R` with your actual implementations.

## Directory Structure

```
project/
├── data/
│   ├── selected_stations_metadata.csv
│   ├── events/
│   ├── raw_ts/
│   ├── flash_events/           # Created by script
│   ├── events_combined/        # Created by script
│   └── clusres/                # Model outputs
├── reports/                    # Generated reports
├── master_analysis.R
└── analysis_report_template.Rmd
```

## Configuration

### Modify Settings in `master_analysis.R`

#### 1. IET Values
```r
IET_VALUES <- c(3, 6, 12, 18, 24)  # Change as needed
```

#### 2. Station Selection
```r
STATIONS <- c(30, 171)  # Add/remove station IDs
```

#### 3. Clustering Parameters
```r
K_VALUES <- 2:8  # Range of cluster numbers to test
CLUS_VARS <- c('severity', 'magnitude', 'duration',
               'intensity', 'time_to_peak', 'flash')
```

#### 4. Regression Settings
```r
REGRESSION_FORMULA <- as.formula(
  precipitation + 10 ~ -1 + (log(magnitude+0.01) + log(duration+0.01) +
                              time_to_peak)*flash
)
FILTER_VALUES <- c(0, 0.21)  # Minimum precipitation thresholds
```

#### 5. Computational Settings
```r
N_CORES <- max(1, floor(parallel::detectCores() * 0.8))
N_BOOTSTRAP <- 100  # Number of bootstrap iterations
```

#### 6. File Paths
Update these paths to match your directory structure:
```r
PATH_STATION_METADATA <- "data/selected_stations_metadata.csv"
PATH_FLASH_RAW <- "../ALDIS/ALDIS_csv/lightning.csv"
PATH_EVENTS <- "data/events/"
PATH_RAW_TS <- "data/raw_ts/"
```

## Usage

### Running the Full Pipeline

```r
# Source the master script
source("master_analysis.R")

# Run with default settings (Steps 2-4, skip Step 1)
main()
```

### Running Specific Steps

```r
# Run all steps including data preparation
main(
  run_step1 = TRUE,      # Generate combined event data
  run_step2 = TRUE,      # Clustering analysis
  run_step3 = TRUE,      # Regression models
  run_step4 = TRUE,      # Generate reports
  run_bootstrap = FALSE  # Skip bootstrap (faster)
)

# Run only clustering and reports (skip Step 1 and 3)
main(
  run_step1 = FALSE,
  run_step2 = TRUE,
  run_step3 = FALSE,
  run_step4 = TRUE,
  run_bootstrap = FALSE
)
```

### Individual Functions

```r
# Run individual steps for specific IET
step2_partitioning_clustering(iet = 6)
step3_clusterwise_regression(iet = 6, run_bootstrap = FALSE)
step4_generate_report(iet = 6)
```

### Manual Report Generation

```r
library(rmarkdown)

render(
  input = "analysis_report_template.Rmd",
  output_file = "report_iet_6.html",
  output_dir = "reports/",
  params = list(
    iet = 6,
    stations = c(30, 171),
    path_results = "data/clusres/"
  )
)
```

## Output Files

### Model Files (saved in `data/clusres/`)

**Clustering Results:**
- `part_iet_[IET]_station_[STATION].RDA`
  - Contains: clustering object, global stability ARIs, segment-level ARIs

**Regression Models:**
- `mod_iet_[IET]_filt_[FILTER]_station_[STATION].RDA`
  - Contains: stepFlexmix model object

**Stability Analysis (if run):**
- `mod_iet_[IET]_filt_[FILTER]_station_[STATION]_stabAn.RDA`
  - Contains: bootstrap ARI results (100+ bootstrap samples)

### Report Files (saved in `reports/`)

- `report_iet_3.html`
- `report_iet_6.html`
- `report_iet_12.html`
- `report_iet_18.html`
- `report_iet_24.html`

Each report includes:
- Global and segment-level ARI plots
- KL divergence comparisons
- FAMD visualizations
- Partial dependence plots
- ECDF plots
- CDF comparison plots
- Model selection criteria (AIC/BIC)

## Important Notes

### Performance Considerations

1. **Bootstrap Analysis:**
   - Set `run_bootstrap = FALSE` for faster testing
   - Bootstrap with 100+ iterations can take hours to days
   - Consider running on HPC for large datasets

2. **Memory Usage:**
   - Large datasets may require substantial RAM
   - Monitor memory usage during clustering and regression steps

3. **Parallel Processing:**
   - Adjust `N_CORES` based on your system
   - Default uses 80% of available cores

### Data Requirements

The script expects the following input files:
- `selected_stations_metadata.csv`: Station metadata with id, lon, lat columns
- `lightning.csv`: Flash/lightning data with longitude, latitude, dt columns
- Event files in `data/events/`: Named as `events_iet_[IET]h_[STATION].csv`
- Raw precipitation files in `data/raw_ts/`: Time series data

### Customization

#### Adding New Visualizations

Edit `analysis_report_template.Rmd` to add sections:

```r
# 10. Your New Section

## 10.1 Your Subsection

```{r your-chunk}
# Your analysis code here
```
```

#### Modifying Clustering Variables

```r
# In master_analysis.R, update:
CLUS_VARS <- c('severity', 'magnitude', 'duration',
               'intensity', 'time_to_peak', 'flash',
               'your_new_variable')
```

#### Changing Regression Formula

```r
REGRESSION_FORMULA <- as.formula(
  precipitation ~ your_predictors + interactions
)
```

## Troubleshooting

### Common Issues

1. **"No event data found"**
   - Check that `data/events_combined/` contains files
   - Verify file naming: `model_input_events_iet_[IET]h_[STATION].RDS`
   - Run Step 1 first if files are missing

2. **Custom function errors**
   - Update function definitions in `master_analysis.R`
   - Ensure your package is loaded with `devtools::load_all()`
   - Check function signatures match expected inputs

3. **Memory errors**
   - Reduce `N_BOOTSTRAP` value
   - Process fewer stations simultaneously
   - Increase system RAM or use HPC

4. **Report generation fails**
   - Check that all required packages are installed
   - Verify `path_results` contains model files
   - Look for error messages in R console

5. **Missing plots in reports**
   - Ensure model files exist for specified IET and stations
   - Check that stability analysis ran successfully
   - Verify data extraction code in report template

## Example Workflow

```r
# 1. First time setup: run data preparation
main(run_step1 = TRUE, run_step2 = FALSE, 
     run_step3 = FALSE, run_step4 = FALSE)

# 2. Quick analysis without bootstrap
main(run_step1 = FALSE, run_step2 = TRUE, 
     run_step3 = TRUE, run_step4 = TRUE, 
     run_bootstrap = FALSE)

# 3. Review reports in reports/ folder

# 4. (Optional) Run bootstrap on HPC for final analysis
main(run_step1 = FALSE, run_step2 = FALSE, 
     run_step3 = TRUE, run_step4 = FALSE, 
     run_bootstrap = TRUE)

# 5. Regenerate reports with stability results
main(run_step1 = FALSE, run_step2 = FALSE, 
     run_step3 = FALSE, run_step4 = TRUE)
```

## Citation

If you use this pipeline in your research, please cite:

```
[Your citation information here]
```

## Support

For questions or issues:
- Open an issue on [GitHub repository]
- Contact: [your email]

## License

[Your license information]

---

**Version:** 1.0  
**Last Updated:** January 2026  
**Author:** [Your name]
