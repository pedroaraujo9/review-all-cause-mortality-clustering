#' Fit all review clustering methods
#'
#' Applies the four clustering methods reviewed in the paper to a common
#' mortality dataset: Hellinger distance with complete linkage
#' (\code{fit_hell_complete}), ILC k-means (\code{fit_ilc_kmeans}),
#' PCA-based fuzzy c-means (\code{fit_pca_fuzzy}), and functional k-means
#' on life expectancy at birth (\code{fit_func_kmeans}). A combined plot of
#' the internal quality metrics across all methods is also produced.
#'
#' @param data A data frame containing period life table data for multiple
#'   countries. Must include columns \code{country}, \code{year}, \code{age},
#'   \code{dx}, \code{mx}, \code{qx}, and \code{ex}.
#'
#' @return An object of class \code{"review_fit"}, which is a list with
#'   the following components:
#'   \describe{
#'     \item{\code{review_fit}}{A named list with elements \code{hell_complete},
#'       \code{ilc_kmeans}, \code{pca_fuzzy}, and \code{func_kmeans}, each
#'       being the output of the corresponding fitting function.}
#'     \item{\code{data}}{The input data frame passed to \code{data}.}
#'     \item{\code{internal_metrics_plot}}{A \code{patchwork} plot combining
#'       the internal quality metrics plots from all four methods.}
#'   }
#'
#' @seealso \code{\link{fit_hell_complete}}, \code{\link{fit_ilc_kmeans}},
#'   \code{\link{fit_pca_fuzzy}}, \code{\link{fit_func_kmeans}}
#'
#' @export
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
