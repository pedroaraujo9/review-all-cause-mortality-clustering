fit_methods = function(data) {

  fits = lapply(list(
    fit_hell_complete,
    fit_ilc_kmeans,
    fit_pca_fuzzy,
    fit_func_kmeans
  ), function(f) f(data))

}
