fit_pca_fuzzy = function(data) {

  qx_tidy = data %>%
    dplyr::select(country, year, age, qx) %>%
    dplyr::filter(age < 110)

  qx = qx_tidy %>%
    tidyr::spread(age, qx)

  logitqx = qx_tidy %>%
    tibble::as_tibble() %>%
    dplyr::mutate(year_age = paste0("v_", year, "_", age)) %>%
    dplyr::select(country, year_age, qx) %>%
    dplyr::mutate(qx = log(qx/(1-qx))) %>%
    tidyr::spread(year_age, qx)

  qx_matrix = logitqx %>% dplyr::select(-country) %>% as.matrix() %>% scale()
  dim(qx_matrix)

  set.seed(1)
  qx_eigen_dec = qx_matrix %>% cor() %>% eigen()
  lambda_cumsum = cumsum(qx_eigen_dec$values)/sum((qx_eigen_dec$values))
  n_dim = length(lambda_cumsum[lambda_cumsum < 0.9]) + 1

  qxPCS = qx_matrix %*% qx_eigen_dec$vectors[, 1:n_dim]
  dim(qxPCS)

  PCA_fuzzy_fit = lapply(2:10, function(g){
    set.seed(1)
    e1071::cmeans(qxPCS, centers = g, m = 2)
  })

  PCA_fuzzy_class_matrix = lapply(PCA_fuzzy_fit, function(fit){
    fit$cluster
  }) %>% do.call(cbind, .)

  # metrics
  metrics = lapply(2:10, function(k){

    fcm_result = PCA_fuzzy_fit[[k-1]]

    sil = cluster::silhouette(fcm_result$cluster, dist(qxPCS), FUN = mean)
    avg_sil = mean(sil[, 3])

    partition_coefficient = sum(fcm_result$membership^2) / nrow(qxPCS)
    partition_entropy = -sum(fcm_result$membership * log(fcm_result$membership)) / nrow(qxPCS)

    min_intercluster_dist = min(dist(fcm_result$centers))^2

    xie_beni = sum(apply(fcm_result$membership^2 * rowSums((qxPCS - fcm_result$centers[fcm_result$cluster, ])^2), 1, sum)) /
      (nrow(qxPCS) * min_intercluster_dist)

    global_mean = colMeans(qxPCS)

    c(
      "sil" = avg_sil,
      "PC" = partition_coefficient,
      "xie_beni" = xie_beni
    )

  }) %>%
    do.call(rbind, .) %>%
    as.data.frame() %>%
    mutate(K = 2:10)

  metric_plot = metrics %>%
    gather(metric, val, -K) %>%
    mutate(metric = ifelse(metric == "xie_beni", "Xie-Beni", metric),
           metric = ifelse(metric == "sil", "Silhouette", metric)) %>%
    ggplot(aes(x=K, y=val)) +
    geom_point() +
    geom_line() +
    facet_wrap(. ~ metric, scales = "free") +
    scale_x_continuous(breaks = 2:10) +
    labs(x="Number of clusters", y="Metric")


  class_matrix = PCA_fuzzy_class_matrix

  out = list(
    metric_plot = metric_plot,
    class_matrix = class_matrix,
    fuzzy_fit = PCA_fuzzy_fit
  )

  return(out)
}
