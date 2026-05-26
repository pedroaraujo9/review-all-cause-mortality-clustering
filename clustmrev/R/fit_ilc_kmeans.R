#' Fit ILC k-means clustering
#'
#' Fits a Lee-Carter model to each country's mortality rates to extract the
#' age sensitivity vector (\eqn{\beta_x}), then applies k-means clustering on
#' those vectors. Internal quality metrics are evaluated for 2 to 10 clusters.
#'
#' @param data A data frame containing period life table data for multiple
#'   countries. Must include columns \code{country}, \code{year}, \code{age},
#'   and \code{mx}.
#'
#' @return A list with the following components:
#'   \describe{
#'     \item{\code{internal_metrics_plot}}{A \code{ggplot2} plot of internal
#'       quality metrics (silhouette and Calinski-Harabasz) across 2 to 10
#'       clusters.}
#'     \item{\code{class_matrix}}{An integer matrix of cluster assignments,
#'       with one column per number of clusters (1 to 10) and one row per
#'       country.}
#'     \item{\code{beta_x}}{A numeric matrix of estimated \eqn{\beta_x}
#'       vectors, with one row per country.}
#'   }
#'
#' @seealso \code{\link{fit_review_methods}}
#'
#' @export
fit_ilc_kmeans = function(data) {

  countries = data$country %>% unique() %>% as.character()
  mx_tidy = data %>% dplyr::select(country, year, age, mx)

  ILC = lapply(countries, function(country_name){
    fit = mx_tidy %>%
      dplyr::filter(country == country_name) %>%
      dplyr::select(year, age, mx) %>%
      dplyr::arrange(age, year) %>%
      tibble::as_tibble() %>%
      vital::as_vital(index = year, key = age, .age = "age") %>%
      vital::model(lee_carter = vital::LC(log(mx), scale = TRUE))

    list(kt = fit %>% vital::time_components() %>% .$kt,
         bx = fit %>% vital::age_components() %>% .$bx)
  })

  beta_x = ILC %>% purrr::map(~{.x$bx}) %>% do.call(rbind, .)

  ILC_k_means_fit = fit_hik(
    data = beta_x,
    diss = stats::dist(beta_x),
    G = 2:10,
    n_start = 100,
    n_iters = 1000,
    method = "kmeans",
    seed = 1
  )

  internal_metrics_plot = ILC_k_means_fit$internal_metrics_plot
  class_matrix = cbind(1, ILC_k_means_fit$class_matrix)

  out = list(
    internal_metrics_plot = internal_metrics_plot,
    class_matrix = class_matrix,
    beta_x = beta_x
  )

  return(out)

}
