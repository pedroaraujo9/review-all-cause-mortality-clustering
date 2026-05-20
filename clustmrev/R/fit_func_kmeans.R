fit_func_kmeans = function(data) {

  ex_tidy = data %>% dplyr::select(country, year, age, ex)
  ex = ex_tidy %>% tidyr::spread(age, ex)
  period_range = range(data$year)

  ex_0 = ex %>%
    dplyr::select(country, year , `0`) %>%
    tidyr::spread(year, `0`)

  ex_0_matrix = ex_0 %>%
    dplyr::select(-country) %>%
    as.matrix()

  years = data$year %>% unique()
  basis = fda::create.bspline.basis(
    rangeval = period_range,
    nbasis = 25,
    norder = 3
  )

  md_fd_obj = fda::smooth.basis(
    argvals = years,
    y = t(ex_0_matrix),
    fdParobj = basis
  )

  set.seed(1)
  func_kmeans_fit = fit_hik(
    data = t(md_fd_obj$fd$coefs),
    diss = stats::dist(t(md_fd_obj$fd$coefs)),
    G = 2:10,
    n_start = 100,
    n_iters = 1000,
    method = "kmeans",
    seed = 1
  )

  internal_metrics_plot = func_kmeans_fit$internal_metrics_plot
  class_matrix = cbind(1, func_kmeans_fit$class_matrix)

  out = list(
    internal_metrics_plot = internal_metrics_plot,
    class_matrix = class_matrix
  )

  return(out)

}
