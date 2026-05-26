#' Analyse Hellinger distance complete-linkage clustering results
#'
#' Extracts cluster assignments from a fitted Hellinger distance
#' complete-linkage hierarchical clustering solution and produces a plot of
#' the median age-at-death distribution (\eqn{d_{xi}^{*}}) per cluster,
#' together with 95 percent interval ribbons.
#'
#' @param review_fit A list returned by \code{\link{fit_review_methods}},
#'   containing the fitted clustering models and the original data.
#' @param clust_size An integer specifying the number of clusters to extract
#'   from the \code{class_matrix}.
#' @param new_label An optional integer vector used to relabel the clusters.
#'   If \code{NULL} (default), the original labels are kept.
#'
#' @return A list with the following components:
#'   \describe{
#'     \item{\code{d_cluster_plot}}{A \code{ggplot2} plot of the median
#'       normalised death distribution \eqn{d_{xi}^{*}} by age group for each
#'       cluster, with pointwise 2.5 to 97.5 percent ribbon bands.}
#'     \item{\code{dendrogram_plot}}{The \code{ggplot2} dendrogram produced
#'       during model fitting (passed through from \code{review_fit}).}
#'     \item{\code{cluster}}{A named integer vector of cluster assignments,
#'       one entry per country.}
#'   }
#'
#' @seealso \code{\link{fit_hell_complete}}, \code{\link{fit_review_methods}}
#'
#' @export
analyse_hell_complete = function(review_fit, clust_size, new_label = NULL) {

  cluster = review_fit$review_fit$hell_complete$class_matrix[, clust_size]

  if(!is.null(new_label)) {
    cluster = relabel_cluster(cluster, new_label)
  }

  data = review_fit$data

  d_cluster_plot = data %>%
    dplyr::mutate(dx_norm = sqrt(dx/100000)) %>%
    dplyr::select(country, year, age, dx_norm) %>%
    dplyr::mutate(id = paste0(country, "-", year)) %>%
    dplyr::mutate(class = cluster[country]) %>%
    tibble::as_tibble() %>%
    dplyr::group_by(class, age) %>%
    dplyr::summarise(med = stats::median(dx_norm),
                     li = stats::quantile(dx_norm, 0.025),
                     ui = stats::quantile(dx_norm, 0.975)) %>%
    ggplot2::ggplot(ggplot2::aes(x=age, y=med, color=factor(class))) +
    ggplot2::geom_line(linewidth = 1) +
    ggplot2::geom_ribbon(ggplot2::aes(x=age, ymin = li, ymax=ui, fill=factor(class)),
                         inherit.aes = F, alpha = 0.2) +
    ggplot2::labs(x="Age group", y=latex2exp::TeX("$d_{x i}^{*}$"),
                  fill = "Cluster", color = "Cluster") +
    ggplot2::theme_bw()

  out = list(
    d_cluster_plot = d_cluster_plot,
    dendrogram_plot = review_fit$review_fit$hell_complete$dendrogram_plot,
    cluster = cluster
  )

  return(out)

}
