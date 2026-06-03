test_that("fit works", {

  fit = fit_ilc_kmeans(data = hmd_data) %>% expect_no_error()

})
