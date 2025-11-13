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


plot_graph = function(class_df, seed = 1) {
  
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
    geom_edge_link(aes(alpha = weight), color = "grey60", show.legend = TRUE) +  
    geom_node_point() + 
    geom_node_text(aes(label = name), repel = TRUE, family = "serif") + 
    scale_edge_alpha(name = "Weight") + 
    theme_graph()
  
} 

