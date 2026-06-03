#' Analyse functional k-means clustering results
#'
#' Extracts cluster assignments from a fitted functional k-means solution and
#' produces a plot of life expectancy at birth (\eqn{e_0}) trajectories over
#' time, coloured by cluster.
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
#'     \item{\code{ex_plot}}{A \code{ggplot2} plot of \eqn{e_0} trajectories
#'       over time for each country, coloured by cluster assignment.}
#'     \item{\code{cluster}}{A named integer vector of cluster assignments,
#'       one entry per country.}
#'   }
#'
#' @seealso \code{\link{fit_func_kmeans}}, \code{\link{fit_review_methods}}
#'
#' @export
analyse_func_kmeans = function(review_fit, clust_size, new_label = NULL) {

  cluster = review_fit$review_fit$func_kmeans$class_matrix[, clust_size]

  if(!is.null(new_label)) {
    cluster = relabel_cluster(cluster, new_label)
  }

  ex_tidy = review_fit$data %>% dplyr::select(country, year, age, ex)
  ex = ex_tidy %>% tidyr::spread(age, ex)

  ex_0 = ex %>%
    dplyr::select(country, year , `0`) %>%
    tidyr::spread(year, `0`)

  ex_plot = ex_0 %>%
    dplyr::mutate(class = factor(cluster)) %>%
    tidyr::gather(time, e0, -country, -class) %>%
    tibble::as_tibble() %>%
    dplyr::mutate(time = stringr::str_extract(time, "\\d{1,10}") %>% as.numeric()) %>%
    ggplot2::ggplot(ggplot2::aes(x=time, y=e0, group=country, color=class)) +
    ggplot2::geom_line() +
    ggplot2::labs(x="Period", y=expression(e[0]), color="Cluster")  +
    ggplot2::theme(text = ggplot2::element_text(size = 13)) +
    ggplot2::theme_bw()

  out = list(
    ex_plot = ex_plot,
    cluster = cluster
  )

  return(out)

}
