analyse_func_kmeans = function(review_fit, clust_size) {

  ex_tidy = review_fit$data %>% dplyr::select(country, year, age, ex)
  ex = ex_tidy %>% tidyr::spread(age, ex)
  cluster = review_fit$review_fit$func_kmeans$class_matrix[, clust_size]

  ex_0 = ex %>%
    dplyr::select(country, year , `0`) %>%
    tidyr::spread(year, `0`)

  ex_0 %>%
    mutate(class = factor(cluster)) %>%
    gather(time, e0, -country, -class) %>%
    as_tibble() %>%
    mutate(time = time %>% str_extract("\\d{1,10}") %>% as.numeric()) %>%
    ggplot(aes(x=time, y=e0, group=country, color=class)) +
    geom_line() +
    labs(x="Period", y=expression(e[0]), color="Cluster")  +
    theme(text = element_text(size = 13)) +
    theme_bw()
}
