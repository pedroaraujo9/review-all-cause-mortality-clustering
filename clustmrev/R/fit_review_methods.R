fit_review_methods = function(data) {

  fits = lapply(list(
    fit_hell_complete,
    fit_ilc_kmeans,
    fit_pca_fuzzy,
    fit_func_kmeans
  ), function(f) f(data))

  names(fits) = c("hell_complete", "ilc_kmeans", "pca_fuzzy", "func_kmeans")

  return(fits)

}
