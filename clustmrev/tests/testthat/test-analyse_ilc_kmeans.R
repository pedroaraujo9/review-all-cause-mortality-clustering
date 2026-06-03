test_that("analyse works", {

  fit = fit_review_methods(data = hmd_data)
  analyse_ilc_kmeans(review_fit = fit, clust_size = 2, ages = c(0, 15, 45, 80), periods = c(1960, 1980, 2000, 2010)) %>% expect_no_error()

})
