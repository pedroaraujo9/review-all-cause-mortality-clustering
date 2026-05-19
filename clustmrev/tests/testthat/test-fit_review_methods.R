test_that("function works", {

  fit = fit_review_methods(data = hmd_data) %>% expect_no_error()

})
