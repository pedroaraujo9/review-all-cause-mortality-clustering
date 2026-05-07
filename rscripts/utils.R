fit_nb_clust = function(data = NULL, diss = NULL, method, indices = NULL) {
  
  if(is.null(indices)) {
    
    indices = c(
      "kl", "ch", "hartigan", "ccc", "scott", 
      "marriot", "trcovw", "tracew", "friedman", 
      "rubin", "cindex", "db", "silhouette", "duda", 
      "pseudot2", "beale", "ratkowsky", "ball", 
      "ptbiserial", "gap", "frey", "mcclain", "gamma", "gplus", "tau", 
      "dunn", "hubert", "sdindex", "dindex", "sdbw"
    )
    
  }
  
  
  if(is.null(diss) & method == "kmeans"){
    
    cluster_fit = lapply(indices, function(index){
      
      pdf(NULL) 
      
      fit = safely(NbClust, otherwise = NA)(
        data = data, 
        method = "kmeans", 
        index = index, 
        distance = "euclidean"
      )
      
      dev.off()  
      
      par(mfrow=c(1,1))
      return(fit)
    })
    
  }else if (method == "kmeans"){
    
    cluster_fit = lapply(indices, function(index){
      
      pdf(NULL) 
      
      fit = safely(NbClust, otherwise = NA)(
        data = data, 
        diss = diss, 
        method = "kmeans", 
        index = index, 
        distance = NULL
      )
      
      dev.off()  
      
      return(fit)
    })
    
  }else{
    
    cluster_fit = lapply(indices, function(index){
      
      pdf(NULL) 
      
      fit = safely(NbClust, otherwise = NA)(
        diss = diss, 
        method = method, 
        index = index, 
        distance = NULL
      )
      
      dev.off()  
      
      return(fit)
    })
    
  }
  
  number_cluster = lapply(cluster_fit, function(fit){
    if(is.null(fit$error)){
      if(is.null(fit$result$Best.partition)) {
        return(NA)
      }else{
        return(length(unique(fit$result$Best.partition)))
      }
    }else{
      return(NA)
    }
  }) %>%
    do.call(c, .)
  
  
  table_number_cluster = number_cluster %>% table() 
  opt_n = which.max(table_number_cluster) %>% names() %>% as.numeric()
  
  class = cluster_fit[[which(number_cluster == opt_n)[1]]]$result$Best.partition
  return(class)
}

fit_compare = function(data = NULL, 
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
  
  vote = G[apply(metrics, MARGIN = 2, FUN = which.max)]
  best = table(vote) %>% which.max() %>% names()
  
  out = list(
    metrics = metrics,
    metrics_plot = metrics_plot,
    class_matrix = class_matrix,
    best = class_matrix[, best]
    
  )
  
  return(out)
}


fit_hell_complete = function(data) {
  
  n = data$country %>% unique() %>% length()
  time_unique = data$year %>% unique()
  dH = matrix(0, nrow = n, ncol = n)
  
  for(i in seq_along(time_unique)) {
    
    d_time = data %>%
      filter(year == time_unique[i]) %>%
      mutate(dx_norm = dx/100000) %>% #divide by l0 = 100000
      select(country, age, dx_norm) %>%
      spread(age, dx_norm) %>%
      select(-country) %>%
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
  
  dd = dendro_data(ward_fit)
  
  dendrogram_plot = ggplot() +
    geom_segment(data = dd$segments,
                 aes(x = x, y = y, xend = xend, yend = yend)) +
    geom_text(
      data = dd$labels,
      aes(x = x, y = y - 0.02 * max(dd$segments$y), label = label),
      hjust = 1, angle = 90
    ) + 
    labs(
      x = "Countries",
      y = "Height (dissimilarity)",
    ) +
    scale_x_continuous(breaks = NULL) + 
    scale_y_continuous(limits = c(-0.10, 0.25), breaks = c(0.1, 0.25)) + 
    theme_minimal() 
  
  h_complete_fit = fit_compare(
    diss = dH, G = 2:10, seed = 1, method = "complete"
  )
  
  metric_plot = h_complete_fit$metrics_plot

  h_complete_class_matrix = h_complete_fit$class_matrix

  out = list(
    dendrogram_plot = dendrogram_plot,
    metric_plot = metric_plot,
    class_matrix = h_complete_class_matrix
  )
  
  return(out)
  
}

fit_ilc_kmeans = function(data) {
  
  country_analyzed = data$country %>% unique() %>% as.character()
  mx_tidy = data %>% dplyr::select(country, year, age, mx)
  
  ILC = lapply(countries_analyzed, function(country_name){
    fit = mx_tidy %>%
      filter(country == country_name) %>%
      dplyr::select(year, age, mx) %>%
      dplyr::arrange(age, year) %>%
      as_tibble() %>%
      as_vital(index = year, key = age, .age = "age") %>%
      model(lee_carter = LC(log(mx), scale = T)) 
    
    list(kt = fit %>% time_components() %>% .$kt, 
         bx = fit %>% age_components() %>% .$bx)
  })
  
  kappa_t = ILC %>% purrr::map(~{.x$kt}) %>% do.call(rbind, .) 
  beta_x = ILC %>% purrr::map(~{.x$bx}) %>% do.call(rbind, .)
  
  ILC_k_means_fit = fit_compare(
    data = beta_x, 
    diss = dist(beta_x), 
    G = 2:10, 
    n_start = 100, 
    n_iters = 1000, 
    method = "kmeans",
    seed = 1
  )
  
  metric_plot = ILC_k_means_fit$metrics_plot
  class_matrix = ILC_kmeans_class_matrix = ILC_k_means_fit$class_matrix
  
  out = list(
    metric_plot = metric_plot,
    class_matrix = class_matrix
  )
  
  return(out)
  
}

fit_pca_fuzzy = function(data) {
  
  qx_tidy = data %>%
    dplyr::select(country, year, age, qx) %>%
    dplyr::filter(age < 110)
  
  qx = qx_tidy %>%
    spread(age, qx)
  
  logitqx = qx_tidy %>%
    as_tibble() %>%
    mutate(year_age = paste0("v_", year, "_", age)) %>%
    dplyr::select(country, year_age, qx) %>%
    mutate(qx = log(qx/(1-qx))) %>%
    spread(year_age, qx) 
  
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
    cmeans(qxPCS, centers = g, m = 2)
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
  
  metrics_plot = metrics %>%
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
    class_matrix = h_complete_class_matrix,
    PCA_fuzzy_fit = PCA_fuzzy_fit
  )
  
  return(out)
}

fit_func_kmeans = function(data) {
  
  ex_tidy = data %>% dplyr::select(country, year, age, ex)
  ex = ex_tidy %>% spread(age, ex)
  
  ex_0 = ex %>% 
    dplyr::select(country, year , `0`) %>%
    spread(year, `0`)
  
  ex_0_matrix = ex_0 %>%
    select(-country) %>%
    as.matrix()
  
  years = lt$year %>% unique()
  basis = create.bspline.basis(
    rangeval = period_range, 
    nbasis = 25,  
    norder = 3
  )
  
  md_fd_obj = smooth.basis(
    argvals = years, 
    y = t(ex_0_matrix), 
    fdParobj = basis
  )
  
  plot(md_fd_obj)
  
  set.seed(1)
  func_kmeans_fit = fit_compare(
    data = t(md_fd_obj$fd$coefs), 
    diss = dist(t(md_fd_obj$fd$coefs)), 
    G = 2:10, 
    n_start = 100, 
    n_iters = 1000, 
    method = "kmeans", 
    seed = 1
  )
  
  metrics_plot = func_kmeans_fit$metrics_plot
  class_matrix = func_kmeans_fit$class_matrix
  
  out = list(
    metric_plot = metric_plot,
    class_matrix = h_complete_class_matrix,
    PCA_fuzzy_fit = PCA_fuzzy_fit
  )
  
  return(out)
  
}



plot_graph = function(class_df, seed = 1, show_legend = TRUE) {
  
  match_matrix = matrix(0, n_country, n_country)
  rownames(match_matrix) = colnames(match_matrix) = class_df$country
  
  for (i in 1:n_country) {
    for (j in 1:n_country) {
      match_matrix[i, j] = sum(class_df[i, -1] == class_df[j, -1])
    }
  }
  
  
  country_conn = match_matrix %>%
    as.data.frame() %>%
    mutate(to = rownames(.)) %>%
    gather(from, weight, -to) %>%
    filter(from != to) %>%
    filter(weight != 0)
  
  country_conn$id = NA
  for(i in 1:nrow(country_conn)) {
    country_conn$id[i] = sort(c(country_conn$from[i], country_conn$to[i])) %>% 
      str_flatten(collapse = "-")
  }
  
  country_conn = country_conn %>%
    distinct(id, .keep_all = T) %>%
    select(-id)
  
  country_graph = graph_from_data_frame(country_conn, directed = FALSE)
  country_tbl_graph = as_tbl_graph(country_graph)
  
  set.seed(seed)
  ggraph(country_tbl_graph, layout = "fr") +  # Fruchterman-Reingold layout
    geom_edge_link(aes(alpha = weight), color = "grey60", show.legend = show_legend) +  
    geom_node_point() + 
    geom_node_text(aes(label = name), repel = TRUE, family = "serif") + 
    scale_edge_alpha(name = "Weight") + 
    theme_graph()
  
} 

