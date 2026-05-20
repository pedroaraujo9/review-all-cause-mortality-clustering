analyse_pca_fuzzy = function(review_fit, clust_size) {

  probs = review_fit$review_fit$pca_fuzzy$membership_level[[clust_size]]
  cluster = review_fit$review_fit$pca_fuzzy$class_matrix[, clust_size]
  max_prob = apply(probs, 1, max)

  country_order = data.frame(
    country = review_fit$data$country %>% unique(),
    cluster = cluster,
    max_prob = max_prob
  ) %>%
    arrange(cluster, desc(max_prob)) %>%
    .$country

  probs_df = probs %>%
    as.data.frame() %>%
    mutate(country = unique(review_fit$data$country))

  probs_df %>%
    gather(class, level, -country) %>%
    mutate(country = factor(country, levels = country_order)) %>%
    ggplot(aes(x=level, y=country, fill=class)) +
    geom_bar(stat="identity") +
    labs(x="Membership level", y="Country", fill="Cluster") +
    theme_minimal() +
    geom_vline(xintercept = 0.5, linetype = "dashed")

}
