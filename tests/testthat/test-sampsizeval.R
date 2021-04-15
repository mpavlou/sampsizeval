test_that("Error message for missing prevalence", {
  expect_error(sampsizeval(c=0.77, se_c=0.025, se_cs =0.15, se_cl = 0.15))
})

test_that("Error message for missing C-statistic", {
  expect_error(sampsizeval(p=0.057, se_c=0.025, se_cs =0.15, se_cl = 0.15))
})

test_that("Error message for missing se_c", {
  expect_error(sampsizeval(p=0.057, c=0.77, se_cs =0.15, se_cl = 0.15))
})


test_that("Error message for missing se_cs", {
  expect_error(sampsizeval(p=0.057, c=0.77,  se_c=0.025, se_cl = 0.15))
})


test_that("Error message for missing se_cl", {
  expect_error(sampsizeval(p=0.057, c=0.77,  se_c=0.025, se_cs =0.15))
})


test_that("Output is a list of length 4", {
  a <- sampsizeval(p=0.057, c=0.77, se_c=0.025, se_cs =0.15, se_cl = 0.15)
  expect_equal(length(a), 4)
})


test_that("Sample size calculations are correct", {
  a <- sampsizeval(p=0.057, c=0.77, se_c=0.025, se_cs =0.15, se_cl = 0.15)
  expect_equal(a$size_c_statistic,1599)
  expect_equal(a$size_calibration_slope,934)
  expect_equal(a$size_calibration_large,897)
  expect_equal(a$size_recommended,1599)
})









