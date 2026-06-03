# clustmrev

`clustmrev` is an R package accompanying the preprint
[*Clustering methods for all-cause mortality: a review* (arXiv:2512.04831)](https://arxiv.org/abs/2512.04831).
It implements and applies four clustering methods identified in the
literature for grouping countries by all-cause mortality patterns.

## Methods

| Method | Fit | Analyse |
| --- | --- | --- |
| Hellinger distance with complete linkage | `fit_hell_complete()` | `analyse_hell_complete()` |
| ILC k-means | `fit_ilc_kmeans()` | `analyse_ilc_kmeans()` |
| PCA-based fuzzy c-means | `fit_pca_fuzzy()` | `analyse_pca_fuzzy()` |
| Functional k-means on life expectancy at birth | `fit_func_kmeans()` | `analyse_func_kmeans()` |

`fit_review_methods()` runs all four methods on a common dataset and returns
a combined internal-quality-metrics plot. `plot_clust_graph()` visualises
cluster structure as a graph.

## Data

The package bundles `hmd_data`, a sample of Human Mortality Database period
life tables for 30 countries over 1960–2010, ages 0–110 in 5-year groups.
See `?hmd_data` for details and sources.

## Installation

From the parent repository root:

```r
# install.packages("remotes")
remotes::install_local("clustmrev")
```

Or from the prebuilt tarball:

```r
install.packages("clustmrev_0.0.0.9000.tar.gz", repos = NULL, type = "source")
```

## Example

```r
library(clustmrev)

fit <- fit_review_methods(hmd_data)
fit$internal_metrics_plot
```

## License

MIT. See [LICENSE](LICENSE).
