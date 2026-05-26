test_that("analyse works", {

  fit = fit_review_methods(data = hmd_data)
  analyse_func_kmeans(review_fit = fit, clust_size = 2) %>% expect_no_error()

})
