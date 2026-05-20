analyse_hell_complete = function(review_fit, clust_size) {

  data = review_fit$data
  cluster = review_fit$review_fit$hell_complete$class_matrix[, clust_size]

  d_cluster_plot = data %>%
    mutate(dx_norm = sqrt(dx/100000)) %>%
    select(country, year, age, dx_norm) %>%
    mutate(id = paste0(country, "-", year)) %>%
    mutate(class = cluster[country]) %>%
    as_tibble() %>%
    group_by(class, age) %>%
    summarise(med = median(dx_norm),
              li = quantile(dx_norm, 0.025),
              ui = quantile(dx_norm, 0.975)) %>%
    ggplot(aes(x=age, y=med, color=factor(class))) +
    geom_line(linewidth = 1) +
    geom_ribbon(aes(x=age, ymin = li, ymax=ui, fill=factor(class)),
                inherit.aes = F, alpha = 0.2) +
    labs(x="Age group", y=latex2exp::TeX("$d_{x i}^{*}$"),
         fill = "Cluster", color = "Cluster") +
    theme_bw()

  out = list(
    d_cluster_plot = d_cluster_plot,
    dendrogram_plot = review_fit$review_fit$hell_complete$dendrogram_plot
  )

  return(out)

}
