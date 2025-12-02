# Review all-cause country-level mortality data

This repository contains code and data for a review and analysis of clustering methods applied to all-cause mortality data.

## Project Structure

### `data/`
Contains mortality data used for the analysis:
- **`life_tables_5x1.rds`**: Life table data at 5-year age intervals and 1-year time periods. This is the primary dataset used for clustering analyses.

### `plots/`
Output directory containing all generated visualizations from the analysis, including.

### `rscripts/`
R scripts for data processing and analysis:

#### `review-methods.R`
Main analysis script that:
- Loads and preprocesses life table data
- Implements and compares multiple clustering methodologies:
  - K-means clustering (standard and ILC-based)
  - Functional data clustering
  - Model-based clustering
  - Fuzzy clustering
  - Hierarchical clustering (Ward's method)
  - PCA-based approaches
- Generates visualizations for each method
- Compares clustering results across different distance metrics and approaches
- Tracks publication years of reviewed methods (median year: 2021)

Key dependencies: `tidyverse`, `NbClust`, `cluster`, `vegan`, `dtw`, `mclust`, `fda`, `vital`, `patchwork`, `ggraph`, `tidygraph`, `igraph`, `e1071`, `ggdendro`

#### `utils.R`
Utility functions supporting the main analysis:
- **`fit_nb_clust()`**: Wrapper function for the NbClust package that fits clustering models and evaluates optimal cluster numbers using multiple validation indices (30 different indices available including silhouette, gap statistic, Calinski-Harabasz, Davies-Bouldin, etc.)
- Supports both data matrices and dissimilarity matrices
- Handles k-means clustering with various distance metrics
- Includes error handling for robust execution across different indices

### `review-all-cause-mortality-clustering.Rproj`
RStudio project file for managing the workspace and project settings.

## Getting Started

1. Open the project in RStudio by double-clicking `review-all-cause-mortality-clustering.Rproj`
2. Ensure all required packages are installed (see dependencies in `review-methods.R`)
3. Run `rscripts/review-methods.R` to reproduce the analysis
4. Generated plots will be saved to the `plots/` directory

## Analysis Overview

This project reviews and implements various clustering methodologies applied to mortality data, comparing their effectiveness across different countries and time periods. 
