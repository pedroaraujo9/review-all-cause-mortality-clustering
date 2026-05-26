gv = c(
  ".", "year", "country", "age", "0", "dx", "dx_norm", "x", "y", "xend", "label",
  "value", "mx", "qx", "year_age", "K", "val", "metric", "from", "weight", "id",
  "yend", "n_country", "to", "name", "yend", "med", "li", "ui", "clust_size", "method",
  "time", "e0", "age_group", "bx", "level", "ari", ""
)

utils::globalVariables(gv)

relabel_cluster = function(cluster, new_labels) {
  cluster = as.character(cluster)
  names(new_labels) = as.character(1:length(new_labels))
  as.numeric(new_labels[cluster])
}



