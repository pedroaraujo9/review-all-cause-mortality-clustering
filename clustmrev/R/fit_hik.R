fit_hik = function(data = NULL,
                   diss,
                   G,
                   n_start,
                   n_iters,
                   method = "kmeans",
                   seed) {

  set.seed(seed)

  if(method == "kmeans") {

    class_matrix = lapply(G, function(g){
      stats::kmeans(data, centers = g, iter.max = n_iters, nstart = n_start)$cluster
    }) %>% do.call(cbind, .)

  }else if(method == "ward") {

    ward_fit = stats::hclust(d = diss, method = "ward.D")

    class_matrix = lapply(G, function(g){
      stats::cutree(ward_fit, k = g)
    }) %>% do.call(cbind, .)

  }else if(method == "complete") {

    complete_fit = stats::hclust(d = diss, method = "complete")

    class_matrix = lapply(G, function(g){
      stats::cutree(complete_fit, k = g)
    }) %>% do.call(cbind, .)

  }

  colnames(class_matrix) = G

  internal_metrics = lapply(as.character(G), function(g){

    quality_stats = WeightedCluster::wcClusterQuality(
      diss = diss, clustering = class_matrix[, g]
    )$stats

    c("Silhouette" = quality_stats["ASW"] %>% as.numeric(),
      "CH" = quality_stats["CH"] %>% as.numeric(),
      "PBC" = quality_stats["PBC"] %>% as.numeric()
    )

  }) %>%
    do.call(rbind, .) %>%
    as.data.frame()

  internal_metrics$G = G

  internal_metrics_plot = internal_metrics %>%
    tidyr::gather(metric, value, -G) %>%
    ggplot2::ggplot(ggplot2::aes(x=G, y=value)) +
    ggplot2::geom_point() +
    ggplot2::geom_line() +
    ggplot2::facet_wrap(. ~ metric, scales = "free_y") +
    ggplot2::labs(x="Number of clusters", y="Metric") +
    ggplot2::scale_x_continuous(breaks = G)

  out = list(
    internal_metrics = internal_metrics,
    internal_metrics_plot = internal_metrics_plot,
    class_matrix = class_matrix
  )

  return(out)
}









