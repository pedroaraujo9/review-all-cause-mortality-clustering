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
      kmeans(data, centers = g, iter.max = n_iters, nstart = n_start)$cluster
    }) %>% do.call(cbind, .)

  }else if(method == "ward") {

    ward_fit = hclust(d = diss, method = "ward.D")

    class_matrix = lapply(G, function(g){
      cutree(ward_fit, k = g)
    }) %>% do.call(cbind, .)

  }else if(method == "complete") {

    complete_fit = hclust(d = diss, method = "complete")

    class_matrix = lapply(G, function(g){
      cutree(complete_fit, k = g)
    }) %>% do.call(cbind, .)

  }

  colnames(class_matrix) = G

  metrics = lapply(as.character(G), function(g){

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

  metrics$G = G

  metrics_plot = metrics %>%
    gather(metric, value, -G) %>%
    ggplot(aes(x=G, y=value)) +
    geom_point() +
    geom_line() +
    facet_wrap(. ~ metric, scales = "free_y") +
    labs(x="Number of clusters", y="Metric") +
    scale_x_continuous(breaks = G)

  out = list(
    metrics = metrics,
    metrics_plot = metrics_plot,
    class_matrix = class_matrix
  )

  return(out)
}









