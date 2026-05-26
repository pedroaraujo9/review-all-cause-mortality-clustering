#' clustmrev: Reviewing clustering methods for all-cause mortality
#'
#' @description
#' Provides implementations of four clustering methods reviewed in the
#' literature for grouping countries by all-cause mortality patterns:
#' Hellinger distance with complete linkage, ILC k-means, PCA-based fuzzy
#' c-means, and functional k-means on life expectancy at birth. The package
#' also includes analysis and visualisation functions for each method, along
#' with the Human Mortality Database (`hmd_data`) sample dataset.
#'
#' @details
#' The main entry point is [fit_review_methods()], which fits all four methods
#' simultaneously and returns combined internal quality metric plots.
#' Individual methods can also be fitted independently:
#' \describe{
#'   \item{[fit_hell_complete()]}{Hellinger distance with complete linkage
#'     hierarchical clustering on death distributions.}
#'   \item{[fit_ilc_kmeans()]}{k-means on country-specific Lee-Carter
#'     \eqn{\beta_x} age-sensitivity vectors.}
#'   \item{[fit_pca_fuzzy()]}{Fuzzy c-means on PCA scores derived from
#'     logit-transformed age-specific mortality rates \eqn{q_x}.}
#'   \item{[fit_func_kmeans()]}{Functional k-means on B-spline coefficients
#'     of life expectancy at birth \eqn{e_0}.}
#' }
#' After fitting, the corresponding `analyse_*()` functions extract cluster
#' assignments for a chosen number of clusters and produce diagnostic plots.
#'
#' @seealso [fit_review_methods()], [hmd_data]
#'
#' @importFrom magrittr %>%
#' @keywords package
#' @name clustmrev
"_PACKAGE"
