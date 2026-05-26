#' Analyse PCA-based fuzzy c-means clustering results
#'
#' Extracts cluster assignments from a fitted PCA fuzzy c-means solution and
#' produces a stacked bar chart of country-level membership levels per
#' cluster.
#'
#' @param review_fit A list returned by \code{\link{fit_review_methods}},
#'   containing the fitted clustering models and the original data.
#' @param clust_size An integer specifying the number of clusters to extract
#'   from the \code{class_matrix}.
#' @param new_label An optional integer vector used to relabel the clusters.
#'   If \code{NULL} (default), the original labels are kept.
#' @param country_order An optional character vector giving the desired
#'   ordering of countries on the y-axis of the membership plot. If
#'   \code{NULL} (default), countries are ordered by cluster and then by
#'   descending maximum membership level.
#'
#' @return A list with the following components:
#'   \describe{
#'     \item{\code{prob_plot}}{A \code{ggplot2} stacked bar chart of
#'       membership levels by country and cluster.}
#'     \item{\code{cluster}}{A named integer vector of hard cluster
#'       assignments, one entry per country.}
#'   }
#'
#' @seealso \code{\link{fit_pca_fuzzy}}, \code{\link{fit_review_methods}}
#'
#' @export
analyse_pca_fuzzy = function(review_fit, clust_size, new_label = NULL, country_order = NULL) {

  cluster = review_fit$review_fit$pca_fuzzy$class_matrix[, clust_size]

  if(!is.null(new_label)) {
    cluster = relabel_cluster(cluster, new_label)
  }

  probs = review_fit$review_fit$pca_fuzzy$membership_level[[clust_size]]
  max_prob = apply(probs, 1, max)

  if(is.null(country_order)) {

    country_order = data.frame(
      country = review_fit$data$country %>% unique(),
      cluster = cluster,
      max_prob = max_prob
    ) %>%
      dplyr::arrange(cluster, dplyr::desc(max_prob)) %>%
      .$country

  }


  probs_df = probs %>%
    as.data.frame() %>%
    dplyr::mutate(country = unique(review_fit$data$country))

  prob_plot = probs_df %>%
    tidyr::gather(class, level, -country) %>%
    dplyr::mutate(country = factor(country, levels = country_order)) %>%
    ggplot2::ggplot(ggplot2::aes(x=level, y=country, fill=class)) +
    ggplot2::geom_bar(stat="identity") +
    ggplot2::labs(x="Membership level", y="Country", fill="Cluster") +
    ggplot2::theme_minimal() +
    ggplot2::geom_vline(xintercept = 0.5, linetype = "dashed")

  out = list(
    prob_plot = prob_plot,
    cluster = cluster
  )

  return(out)

}
