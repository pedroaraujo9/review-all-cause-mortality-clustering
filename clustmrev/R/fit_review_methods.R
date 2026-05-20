fit_review_methods = function(data) {

  fits = lapply(list(
    fit_hell_complete,
    fit_ilc_kmeans,
    fit_pca_fuzzy,
    fit_func_kmeans
  ), function(f) f(data))

  names(fits) = c("hell_complete", "ilc_kmeans", "pca_fuzzy", "func_kmeans")

  internal_metrics_plot = patchwork::wrap_plots(
    fits$hell_complete$internal_metrics_plot + ggplot2::ggtitle("Hellinger-Complete") + ggplot2::theme_bw(),
    fits$ilc_kmeans$internal_metrics_plot  + ggplot2::ggtitle("ILC-k-means") + ggplot2::theme_bw(),
    fits$pca_fuzzy$internal_metrics_plot + ggplot2::ggtitle("PCA-fuzzy") + ggplot2::theme_bw(),
    fits$func_kmeans$internal_metrics_plot + ggplot2::ggtitle("Func-k-means") + ggplot2::theme_bw(),
    ncol = 1,
    axis_titles = "collect"
  )


  out = list(
    review_fit = fits,
    data = data,
    internal_metrics_plot = internal_metrics_plot
  )

  class(out) = "review_fit"
  return(out)

}
