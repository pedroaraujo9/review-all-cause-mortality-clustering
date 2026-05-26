#' Fit functional k-means clustering on life expectancy at birth
#'
#' Smooths country-specific life expectancy at birth (\eqn{e_0}) curves over
#' time using B-spline basis functions, then applies k-means clustering on the
#' resulting spline coefficients. Internal quality metrics are evaluated for 2
#' to 10 clusters.
#'
#' @param data A data frame containing period life table data for multiple
#'   countries. Must include columns \code{country}, \code{year}, \code{age},
#'   and \code{ex}.
#'
#' @return A list with the following components:
#'   \describe{
#'     \item{\code{internal_metrics_plot}}{A \code{ggplot2} plot of internal
#'       quality metrics (silhouette and Calinski-Harabasz) across 2 to 10
#'       clusters.}
#'     \item{\code{class_matrix}}{An integer matrix of cluster assignments,
#'       with one column per number of clusters (1 to 10) and one row per
#'       country.}
#'   }
#'
#' @seealso \code{\link{fit_review_methods}}
#'
#' @export
fit_func_kmeans = function(data) {

  ex_tidy = data %>% dplyr::select(country, year, age, ex)
  ex = ex_tidy %>% tidyr::spread(age, ex)
  period_range = range(data$year)

  ex_0 = ex %>%
    dplyr::select(country, year , `0`) %>%
    tidyr::spread(year, `0`)

  ex_0_matrix = ex_0 %>%
    dplyr::select(-country) %>%
    as.matrix()

  years = data$year %>% unique()
  basis = fda::create.bspline.basis(
    rangeval = period_range,
    nbasis = floor(diff(range(period_range))/2),
    norder = 3
  )

  md_fd_obj = fda::smooth.basis(
    argvals = years,
    y = t(ex_0_matrix),
    fdParobj = basis
  )

  set.seed(1)
  func_kmeans_fit = fit_hik(
    data = t(md_fd_obj$fd$coefs),
    diss = stats::dist(t(md_fd_obj$fd$coefs)),
    G = 2:10,
    n_start = 100,
    n_iters = 1000,
    method = "kmeans",
    seed = 1
  )

  internal_metrics_plot = func_kmeans_fit$internal_metrics_plot
  class_matrix = cbind(1, func_kmeans_fit$class_matrix)

  out = list(
    internal_metrics_plot = internal_metrics_plot,
    class_matrix = class_matrix
  )

  return(out)

}
