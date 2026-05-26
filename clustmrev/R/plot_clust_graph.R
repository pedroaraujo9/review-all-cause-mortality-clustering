#' Plot a country co-clustering graph for a single cluster solution
#'
#' Constructs an undirected weighted graph where nodes are countries and edge
#' weights represent the number of methods (out of four) that assign the two
#' countries to the same cluster. The graph is rendered using a
#' Fruchterman-Reingold layout.
#'
#' @param review_fit An object of class \code{"review_fit"} returned by
#'   \code{\link{fit_review_methods}}.
#' @param clust_number A single integer specifying which column of each
#'   method's \code{class_matrix} to use (i.e., the number of clusters).
#' @param seed Integer random seed passed to \code{set.seed} before computing
#'   the graph layout. Default is \code{1}.
#' @param show_legend Logical; whether to display the edge-weight legend.
#'   Default is \code{TRUE}.
#'
#' @return A \code{ggraph} plot object.
#'
#' @keywords internal
plot_graph = function(review_fit, clust_number, seed = 1, show_legend = TRUE) {

  class_df = data.frame(
    country = review_fit$data$country %>% unique(),
    review_fit$review_fit$hell_complete$class_matrix[, clust_number],
    review_fit$review_fit$ilc_kmeans$class_matrix[, clust_number],
    review_fit$review_fit$pca_fuzzy$class_matrix[, clust_number],
    review_fit$review_fit$func_kmeans$class_matrix[, clust_number]
  )

  n_country = nrow(class_df)
  match_matrix = matrix(0, n_country, n_country)
  rownames(match_matrix) = colnames(match_matrix) = class_df$country

  for (i in 1:n_country) {
    for (j in 1:n_country) {
      match_matrix[i, j] = sum(class_df[i, -1] == class_df[j, -1])
    }
  }

  country_conn = match_matrix %>%
    as.data.frame() %>%
    dplyr::mutate(to = rownames(.)) %>%
    tidyr::gather(from, weight, -to) %>%
    dplyr::filter(from != to) %>%
    dplyr::filter(weight != 0)

  country_conn$id = NA
  for(i in 1:nrow(country_conn)) {
    country_conn$id[i] = sort(c(country_conn$from[i], country_conn$to[i])) %>%
      stringr::str_flatten(collapse = "-")
  }

  country_conn = country_conn %>%
    dplyr::distinct(id, .keep_all = T) %>%
    dplyr::select(-id)

  country_graph = igraph::graph_from_data_frame(country_conn, directed = FALSE)
  country_tbl_graph = tidygraph::as_tbl_graph(country_graph)

  set.seed(seed)
  out = ggraph::ggraph(country_tbl_graph, layout = "fr") +  # Fruchterman-Reingold layout
    ggraph::geom_edge_link(ggplot2::aes(alpha = weight), color = "grey60", show.legend = show_legend) +
    ggraph::geom_node_point() +
    ggraph::geom_node_text(ggplot2::aes(label = name), repel = TRUE, family = "serif") +
    ggraph::scale_edge_alpha(name = "Weight") +
    ggraph::theme_graph()

  return(out)

}

#' Plot country co-clustering graphs across multiple cluster solutions
#'
#' Calls \code{plot_graph} for each requested cluster size and combines the
#' resulting panels into a single \code{patchwork} figure. The edge weight
#' legend is shown only on the last panel.
#'
#' @param review_fit An object of class \code{"review_fit"} returned by
#'   \code{\link{fit_review_methods}}.
#' @param clust_number An integer vector of cluster sizes to display, e.g.
#'   \code{c(3, 4, 5)}.
#' @param seed An integer vector of random seeds (one per element of
#'   \code{clust_number}) passed to the Fruchterman-Reingold layout.
#'
#' @return A \code{patchwork} plot combining one \code{ggraph} panel per
#'   element of \code{clust_number}.
#'
#' @seealso \code{\link{fit_review_methods}}
#'
#' @export
plot_clust_graph = function(review_fit, clust_number, seed) {

  plots = lapply(seq_along(clust_number), function(g){
    plot_graph(
      review_fit,
      clust_number = clust_number[g],
      seed = seed[g],
      show_legend = ifelse(g == length(clust_number), TRUE, FALSE)
    ) + ggplot2::ggtitle(paste0(clust_number[g], " clusters"))
  })

  patchwork::wrap_plots(plots)

}


