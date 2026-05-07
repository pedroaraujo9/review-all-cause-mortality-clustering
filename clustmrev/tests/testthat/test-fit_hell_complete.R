test_that("fit works", {

  fit = fit_hell_complete(data = hmd_data) %>%
    expect_no_error()

})
