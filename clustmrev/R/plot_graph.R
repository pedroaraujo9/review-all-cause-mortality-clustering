plot_graph = function(class_df, seed = 1, show_legend = TRUE) {

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
