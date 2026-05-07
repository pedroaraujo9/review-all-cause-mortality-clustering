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

  max_prob = lapply(PCA_fuzzy_fit, function(fit){
    fit$membership %>% apply(1, max)
  }) %>% do.call(cbind, .)

  # metrics
  metrics = lapply(2:10, function(k){

    fcm_result = PCA_fuzzy_fit[[k-1]]

    # Fuzzy silhouette: weight crisp silhouette widths by confidence gap
    # between the largest and second-largest memberships per observation.
    sil = cluster::silhouette(fcm_result$cluster, dist(qxPCS))
    sorted_membership = t(apply(fcm_result$membership, 1, sort, decreasing = TRUE))
    membership_gap = sorted_membership[, 1] - sorted_membership[, 2]
    if (sum(membership_gap) > 0) {
      avg_sil = sum(membership_gap * sil[, 3]) / sum(membership_gap)
    } else {
      avg_sil = mean(sil[, 3])
    }

    partition_coefficient = sum(fcm_result$membership^2) / nrow(qxPCS)
    partition_entropy = -sum(fcm_result$membership * log(fcm_result$membership)) / nrow(qxPCS)

    min_intercluster_dist = min(dist(fcm_result$centers))^2

    dist_sq_to_centers = sapply(seq_len(nrow(fcm_result$centers)), function(j) {
      rowSums((qxPCS - matrix(fcm_result$centers[j, ], nrow = nrow(qxPCS), ncol = ncol(qxPCS), byrow = TRUE))^2)
    })

    xie_beni = sum((fcm_result$membership^2) * dist_sq_to_centers) /
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
    dplyr::mutate(K = 2:10)

  metric_plot = metrics %>%
    tidyr::gather(metric, val, -K) %>%
    dplyr::mutate(metric = ifelse(metric == "xie_beni", "Xie-Beni", metric),
           metric = ifelse(metric == "sil", "Silhouette", metric)) %>%
    ggplot2::ggplot(ggplot2::aes(x=K, y=val)) +
    ggplot2::geom_point() +
    ggplot2::geom_line() +
    ggplot2::facet_wrap(. ~ metric, scales = "free") +
    ggplot2::scale_x_continuous(breaks = 2:10) +
    ggplot2::labs(x="Number of clusters", y="Metric")


  class_matrix = PCA_fuzzy_class_matrix

  out = list(
    metric_plot = metric_plot,
    class_matrix = class_matrix,
    max_prob = max_prob
  )

  return(out)
}
