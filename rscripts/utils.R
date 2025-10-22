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
    
    fpc_stats = fpc::cluster.stats(
      d = diss, clustering = class_matrix[, g], silhouette = T, G2 = T, G3 = T
    )
    
    
    c("Silhouette" = fpc_stats$avg.silwidth, 
      "Dunn" = fpc_stats$dunn2, 
      "CH" = fpc_stats$ch,
      "G2" = fpc_stats$g2
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