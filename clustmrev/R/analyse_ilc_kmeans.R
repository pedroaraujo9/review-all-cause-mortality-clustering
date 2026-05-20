analyse_ilc_kmeans = function(review_fit, clust_size) {

  beta_x = review_fit$review_fit$ilc_kmeans$beta_x
  countries = review_fit$data$country %>% unique()
  cluster = review_fit$review_fit$ilc_kmeans$class_matrix[, clust_size]
  mx_tidy = review_fit$data %>% dplyr::select(country, year, age, mx)
  colnames(beta_x) = review_fit$data$age %>% unique()

  beta_plot = beta_x %>%
    as.data.frame() %>%
    mutate(country = countries , class = factor(cluster)) %>%
    gather(age_group, bx, -country, -class) %>%
    as_tibble() %>%
    mutate(age_group = age_group %>% as.numeric()) %>%
    ggplot(aes(x=age_group, y=bx, group=country, color=class)) +
    geom_line() +
    labs(x="Age group", y=latex2exp::TeX("$\\beta_{x i}$"), color="Cluster") +
    geom_hline(yintercept = 0, linetype = "dashed") +
    theme_bw()

  mx_ilc = mx_tidy %>%
    as_tibble() %>%
    left_join(
      data.frame(
        country = countries,
        class = factor(cluster)
      ),
      by = "country"
    ) %>%
    mutate(class = paste0("Cluster ", class))

  mx_over_time_plot = mx_ilc %>%
    group_by(year, class, age) %>%
    summarise(mx = mean(log(mx))) %>%
    filter(age %in% c(0, 15, 45, 80)) %>%
    mutate(age = ifelse(age == 0, paste0("Age group [0, 1)"),
                        paste0("Age group [", age, ", ", age + 5, ")"))) %>%
    mutate(class = class %>% str_remove("Cluster ")) %>%
    ggplot(aes(x=year, y=(mx), color=class)) +
    geom_line(alpha = 0.8) +
    facet_wrap(. ~ age, scales = "free") +
    labs(x = "Period", y=latex2exp::TeX("Average $\\log(m_{x})$"),
         color = "Cluster") +
    theme(text = element_text(size = 17)) +
    theme_bw()

  mx_curve_plot = mx_ilc %>%
    group_by(year, class, age) %>%
    summarise(mx = mean(log(mx))) %>%
    filter(year %in% c(1960, 1980, 2000, 2010)) %>%
    ggplot(aes(x=age, y=(mx), color=class)) +
    geom_line() +
    facet_wrap(. ~ year, scales = "free") +
    labs(x = "Age group", y=latex2exp::TeX("Average $\\log(m_{x})$"),
         color = "Cluster") +
    theme_bw()

  out = list(
    beta_plot = beta_plot,
    mx_over_time_plot = mx_over_time_plot,
    mx_curve_plot = mx_curve_plot
  )

  return(out)

}
