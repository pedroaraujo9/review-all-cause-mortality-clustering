#' Analyse ILC k-means clustering results
#'
#' Extracts cluster assignments from a fitted ILC k-means solution and
#' produces three diagnostic plots: country-specific \eqn{\beta_x} curves
#' coloured by cluster, average \eqn{\log(m_x)} over time for selected age
#' groups, and average \eqn{\log(m_x)} age curves for selected periods.
#'
#' @param review_fit A list returned by \code{\link{fit_review_methods}},
#'   containing the fitted clustering models and the original data.
#' @param clust_size An integer specifying the number of clusters to extract
#'   from the \code{class_matrix}.
#' @param ages A numeric vector of age groups to display in the
#'   \code{mx_over_time_plot} (e.g. \code{c(0, 15, 45, 80)}).
#' @param periods A numeric vector of calendar years to display in the
#'   \code{mx_curve_plot} (e.g. \code{c(1960, 1980, 2000, 2010)}).
#' @param new_label An optional integer vector used to relabel the clusters.
#'   If \code{NULL} (default), the original labels are kept.
#'
#' @return A list with the following components:
#'   \describe{
#'     \item{\code{beta_plot}}{A \code{ggplot2} plot of country-specific
#'       \eqn{\beta_x} vectors by age group, coloured by cluster.}
#'     \item{\code{mx_over_time_plot}}{A \code{ggplot2} faceted plot of
#'       cluster-average \eqn{\log(m_x)} over time for age groups 0, 15, 45,
#'       and 80.}
#'     \item{\code{mx_curve_plot}}{A \code{ggplot2} faceted plot of
#'       cluster-average \eqn{\log(m_x)} age curves for the periods 1960,
#'       1980, 2000, and 2010.}
#'     \item{\code{cluster}}{A named integer vector of cluster assignments,
#'       one entry per country.}
#'   }
#'
#' @seealso \code{\link{fit_ilc_kmeans}}, \code{\link{fit_review_methods}}
#'
#' @export
analyse_ilc_kmeans = function(review_fit, clust_size, ages, periods, new_label = NULL) {

  cluster = review_fit$review_fit$ilc_kmeans$class_matrix[, clust_size]

  if(!is.null(new_label)) {
    cluster = relabel_cluster(cluster, new_label)
  }

  beta_x = review_fit$review_fit$ilc_kmeans$beta_x
  countries = review_fit$data$country %>% unique()
  mx_tidy = review_fit$data %>% dplyr::select(country, year, age, mx)
  colnames(beta_x) = review_fit$data$age %>% unique()

  beta_plot = beta_x %>%
    as.data.frame() %>%
    dplyr::mutate(country = countries , class = factor(cluster)) %>%
    tidyr::gather(age_group, bx, -country, -class) %>%
    tibble::as_tibble() %>%
    dplyr::mutate(age_group = age_group %>% as.numeric()) %>%
    ggplot2::ggplot(ggplot2::aes(x=age_group, y=bx, group=country, color=class)) +
    ggplot2::geom_line() +
    ggplot2::labs(x="Age group", y=latex2exp::TeX("$\\beta_{x i}$"), color="Cluster") +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed") +
    ggplot2::theme_bw()

  mx_ilc = mx_tidy %>%
    tibble::as_tibble() %>%
    dplyr::left_join(
      data.frame(
        country = countries,
        class = factor(cluster)
      ),
      by = "country"
    ) %>%
    dplyr::mutate(class = paste0("Cluster ", class))

  mx_over_time_plot = mx_ilc %>%
    dplyr::group_by(year, class, age) %>%
    dplyr::summarise(mx = mean(log(mx))) %>%
    dplyr::filter(age %in% ages) %>%
    dplyr::mutate(age = ifelse(age == 0, paste0("Age group [0, 1)"),
                               paste0("Age group [", age, ", ", age + 5, ")"))) %>%
    dplyr::mutate(class = stringr::str_remove(class, "Cluster ")) %>%
    ggplot2::ggplot(ggplot2::aes(x=year, y=(mx), color=class)) +
    ggplot2::geom_line(alpha = 0.8) +
    ggplot2::facet_wrap(. ~ age, scales = "free") +
    ggplot2::labs(x = "Period", y=latex2exp::TeX("Average $\\log(m_{x})$"),
                  color = "Cluster") +
    ggplot2::theme(text = ggplot2::element_text(size = 17)) +
    ggplot2::theme_bw()

  mx_curve_plot = mx_ilc %>%
    dplyr::group_by(year, class, age) %>%
    dplyr::summarise(mx = mean(log(mx))) %>%
    dplyr::filter(year %in% periods) %>%
    ggplot2::ggplot(ggplot2::aes(x=age, y=(mx), color=class)) +
    ggplot2::geom_line() +
    ggplot2::facet_wrap(. ~ year, scales = "free") +
    ggplot2::labs(x = "Age group", y=latex2exp::TeX("Average $\\log(m_{x})$"),
                  color = "Cluster") +
    ggplot2::theme_bw()

  out = list(
    beta_plot = beta_plot,
    mx_over_time_plot = mx_over_time_plot,
    mx_curve_plot = mx_curve_plot,
    cluster = cluster
  )

  return(out)

}
