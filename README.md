# Review of all-cause mortality clustering methods

This repository accompanies the preprint
[*Clustering methods for all-cause mortality: a review* (arXiv:2512.04831)](https://arxiv.org/abs/2512.04831),
a review of clustering methods applied to country-level all-cause mortality
data. It bundles the reviewed methods into an R package (`clustmrev`) and
provides the scripts and data used to produce the analyses and figures.

## The `clustmrev` package

[clustmrev/](clustmrev/) is an R package that implements and applies four
clustering methods identified in the literature for grouping countries by
all-cause mortality patterns:

- **Hellinger distance with complete linkage** — `fit_hell_complete()`
- **ILC k-means** — `fit_ilc_kmeans()`
- **PCA-based fuzzy c-means** — `fit_pca_fuzzy()`
- **Functional k-means on life expectancy at birth** — `fit_func_kmeans()`

Each method ships with a matching `analyse_*()` function for downstream
visualisation and diagnostics. `fit_review_methods()` runs all four on a
common dataset and returns a combined internal-quality-metrics plot, and
`plot_clust_graph()` visualises cluster structure as a graph. A sample
dataset (`hmd_data`) of Human Mortality Database period life tables for 30
countries (1960–2010, ages 0–110 in 5-year groups) is included.

### Installation

From the repository root:

```r
# install.packages("remotes")
remotes::install_local("clustmrev")
```

A prebuilt source tarball (`clustmrev_0.0.0.9000.tar.gz`) is also provided.

## Repository structure

- [clustmrev/](clustmrev/) — The `clustmrev` R package (source, data, docs, tests).
- [data/](data/) — Preprocessed life-table datasets (`.rds`) for different
  sex, year range, and age range combinations, plus `life_tables_5x1.rds`.
- [rscripts/](rscripts/) — Analysis scripts:
  - [process-data.R](rscripts/process-data.R) — builds the datasets in
    [data/](data/) from the Human Mortality Database.
  - [review-methods.R](rscripts/review-methods.R) — runs the four reviewed
    clustering methods via `clustmrev` and produces the figures.
  - [utils.R](rscripts/utils.R) — helper functions used by the analysis
    scripts.
- [plots/](plots/) — Output directory for generated figures.
- [review-all-cause-mortality-clustering.Rproj](review-all-cause-mortality-clustering.Rproj) —
  RStudio project file.

## Getting started

1. Open the project in RStudio via `review-all-cause-mortality-clustering.Rproj`.
2. Install `clustmrev` (see above) and its dependencies.
3. Run [rscripts/review-methods.R](rscripts/review-methods.R) to reproduce
   the analysis; figures are written to [plots/](plots/).

