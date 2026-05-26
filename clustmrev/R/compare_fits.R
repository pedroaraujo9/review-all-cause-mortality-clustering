compute_ari = function(class_matrix_1, class_matrix_2) {

  G = ncol(class_matrix_1)

  ARI = lapply(1:G, function(g){
    mclust::adjustedRandIndex(
      class_matrix_1[, g],
      class_matrix_2[, g]
    )
  }) %>% do.call(c, .)

  return(ARI)
}

compare_fits = function(review_fit_1, review_fit_2) {

  method_names = names(review_fit_1$review_fit)

  ARIs = lapply(method_names, function(name){
    compute_ari(
      class_matrix_1 = review_fit_1$review_fit[[name]]$class_matrix,
      class_matrix_2 = review_fit_2$review_fit[[name]]$class_matrix
    )
  }) %>% do.call(rbind, .)

  ARIs = tibble::tibble(as.data.frame(ARIs))
  colnames(ARIs) = paste0("G", 1:ncol(ARIs))

  ARIs = ARIs %>%
    dplyr::mutate(method = method_names) %>%
    dplyr::select(method, dplyr::everything()) %>%
    as.data.frame()

  ARI_plot = ARIs %>%
    tidyr::gather(clust_size, ari, -method) %>%
    dplyr::mutate(clust_size = as.numeric(gsub("G", "", clust_size))) %>%
    dplyr::mutate(method = method %>%
                    stringr::str_replace("func_kmeans", "Func-k-means") %>%
                    stringr::str_replace("hell_complete", "Hellinger-Complete") %>%
                    stringr::str_replace("ilc_kmeans", "ILC-k-means") %>%
                    stringr::str_replace("pca_fuzzy", "PCA-fuzzy")) %>%
    dplyr::mutate(method = factor(method, levels = c("Hellinger-Complete", "ILC-k-means", "PCA-fuzzy", "Func-k-means"))) %>%
    ggplot2::ggplot(ggplot2::aes(x=clust_size, y=ari, group=method)) +
    ggplot2::geom_line() +
    ggplot2::geom_point() +
    ggplot2::scale_x_continuous(breaks = 1:ncol(ARIs)) +
    ggplot2::facet_wrap(~method) +
    ggplot2::theme_bw()

  out = list(
    ARI = ARIs,
    ARI_plot = ARI_plot
  )

  return(out)

}

