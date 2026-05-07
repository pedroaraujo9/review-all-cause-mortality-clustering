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
      dist() %>%
      as.matrix() %>%
      `/`(sqrt(2))

    dH = dH + d_time/length(time_unique)

  }

  colnames(dH) = rownames(dH) = data$country %>% unique()
  dH = as.dist(dH)

  ward_fit = hclust(dH, method = "complete")

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
    ggplot2::scale_y_continuous(limits = c(-0.10, 0.25), breaks = c(0.1, 0.25)) +
    ggplot2::theme_minimal()

  h_complete_fit = fit_hik(
    diss = dH, G = 2:10, seed = 1, method = "complete"
  )

  metric_plot = h_complete_fit$metrics_plot
  class_matrix = h_complete_fit$class_matrix

  out = list(
    dendrogram_plot = dendrogram_plot,
    metric_plot = metric_plot,
    class_matrix = class_matrix
  )

  return(out)

}
