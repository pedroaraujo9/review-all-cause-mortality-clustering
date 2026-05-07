fit_ilc_kmeans = function(data) {

  countries = data$country %>% unique() %>% as.character()
  mx_tidy = data %>% dplyr::select(country, year, age, mx)

  ILC = lapply(countries, function(country_name){
    fit = mx_tidy %>%
      dplyr::filter(country == country_name) %>%
      dplyr::select(year, age, mx) %>%
      dplyr::arrange(age, year) %>%
      tibble::as_tibble() %>%
      vital::as_vital(index = year, key = age, .age = "age") %>%
      vital::model(lee_carter = vital::LC(log(mx), scale = TRUE))

    list(kt = fit %>% vital::time_components() %>% .$kt,
         bx = fit %>% vital::age_components() %>% .$bx)
  })

  kappa_t = ILC %>% purrr::map(~{.x$kt}) %>% do.call(rbind, .)
  beta_x = ILC %>% purrr::map(~{.x$bx}) %>% do.call(rbind, .)

  ILC_k_means_fit = fit_hik(
    data = beta_x,
    diss = dist(beta_x),
    G = 2:10,
    n_start = 100,
    n_iters = 1000,
    method = "kmeans",
    seed = 1
  )

  metric_plot = ILC_k_means_fit$metrics_plot
  class_matrix = ILC_k_means_fit$class_matrix

  out = list(
    metric_plot = metric_plot,
    class_matrix = class_matrix
  )

  return(out)

}
