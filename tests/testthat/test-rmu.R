test_that("rho should be in [-1,1]", {
  expect_warning(rmu(0,1,
                     rep(0.25, 5),
                     rep(0.125, 4),
                     0.1, 0.25, 0.1, rho = -1.1, h=1),
                 "NaNs produced")
})
