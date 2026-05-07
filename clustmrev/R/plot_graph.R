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
    mutate(to = rownames(.)) %>%
    gather(from, weight, -to) %>%
    filter(from != to) %>%
    filter(weight != 0)

  country_conn$id = NA
  for(i in 1:nrow(country_conn)) {
    country_conn$id[i] = sort(c(country_conn$from[i], country_conn$to[i])) %>%
      str_flatten(collapse = "-")
  }

  country_conn = country_conn %>%
    distinct(id, .keep_all = T) %>%
    select(-id)

  country_graph = graph_from_data_frame(country_conn, directed = FALSE)
  country_tbl_graph = as_tbl_graph(country_graph)

  set.seed(seed)
  out = ggraph(country_tbl_graph, layout = "fr") +  # Fruchterman-Reingold layout
    geom_edge_link(aes(alpha = weight), color = "grey60", show.legend = show_legend) +
    geom_node_point() +
    geom_node_text(aes(label = name), repel = TRUE, family = "serif") +
    scale_edge_alpha(name = "Weight") +
    theme_graph()

  return(out)

}
