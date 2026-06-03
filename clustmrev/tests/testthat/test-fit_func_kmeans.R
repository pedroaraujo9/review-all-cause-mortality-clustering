test_that("fit works", {

  fit = fit_func_kmeans(data = hmd_data) %>% expect_no_error()

})
