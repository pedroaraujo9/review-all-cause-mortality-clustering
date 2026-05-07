test_that("fit works", {

  fit = fit_pca_fuzzy(data = hmd_data) %>% expect_no_error()

})
