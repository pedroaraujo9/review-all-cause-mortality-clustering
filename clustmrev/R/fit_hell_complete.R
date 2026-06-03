#' Fit Hellinger distance with complete linkage clustering
#'
#' Computes the average Hellinger distance matrix across all years between
#' countries' death distributions and applies hierarchical clustering with
#' complete linkage. Internal quality metrics (silhouette, Calinski-Harabasz)
#' are evaluated for 2 to 10 clusters.
#'
#' @param data A data frame containing period life table data for multiple
#'   countries. Must include columns \code{country}, \code{year}, \code{age},
#'   and \code{dx}.
#'
#' @return A list with the following components:
#'   \describe{
#'     \item{\code{dendrogram_plot}}{A \code{ggplot2} dendrogram of the
#'       complete-linkage hierarchical clustering solution.}
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
fit_hell_complete = function(data) {

  n = data$country %>% unique() %>% length()
  time_unique = data$year %>% unique()
  dH = matrix(0, nrow = n, ncol = n)

  for(i in seq_along(time_unique)) {

    d_time = data %>%
      dplyr::filter(year == time_unique[i]) %>%
      dplyr::mutate(dx_norm = dx/100000) %>%
      dplyr::select(country, age, dx_norm) %>%
      tidyr::spread(age, dx_norm) %>%
      dplyr::select(-country) %>%
      as.matrix() %>%
      sqrt() %>%
      stats::dist() %>%
      as.matrix() %>%
      `/`(sqrt(2))

    dH = dH + d_time/length(time_unique)

  }

  colnames(dH) = rownames(dH) = data$country %>% unique()
  dH = stats::as.dist(dH)

  ward_fit = stats::hclust(dH, method = "complete")

  dd = ggdendro::dendro_data(ward_fit)

  dendrogram_plot = ggplot2::ggplot() +
    ggplot2::geom_segment(data = dd$segments,
                          ggplot2::aes(x = x, y = y, xend = xend, yend = yend)) +
    ggplot2::geom_text(
      data = dd$labels,
      ggplot2::aes(x = x, y = y - 0.02 * max(dd$segments$y), label = label),
      hjust = 1, angle = 90
    ) +
    ggplot2::labs(
      x = "Countries",
      y = "Height (dissimilarity)",
    ) +
    ggplot2::scale_x_continuous(breaks = NULL) +
    ggplot2::theme_minimal()

  h_complete_fit = fit_hik(
    diss = dH, G = 2:10, seed = 1, method = "complete"
  )

  internal_metrics_plot = h_complete_fit$internal_metrics_plot
  class_matrix = cbind(1, h_complete_fit$class_matrix)

  out = list(
    dendrogram_plot = dendrogram_plot,
    internal_metrics_plot = internal_metrics_plot,
    class_matrix = class_matrix
  )

  return(out)

}
