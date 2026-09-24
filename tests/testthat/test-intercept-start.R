# The start a family reads off the data for a regression's intercepts.

test_that("the base method asks for nothing", {
  expect_identical(distrib_intercept_start(poisson_distrib(), c(0, 1, 3)), list())
  expect_identical(distrib_intercept_start(gaussian1_distrib(), rnorm(5)), list())
})

test_that("zero_inflated() starts its mixing weight at the proportion of zeros", {
  d <- zero_inflated(negbin2_distrib())
  y <- c(0, 0, 0, 1, 4, 2, 0, 7)
  s <- distrib_intercept_start(d, y)
  expect_named(s, "zi")
  expect_equal(s$zi, 0.5)
  # on the parameter scale, whatever the link: the layer maps it
  dp <- zero_inflated(negbin2_distrib(), link_zi = linkfunctions7::probit_link())
  expect_equal(distrib_intercept_start(dp, y)$zi, 0.5)
  # kept inside the open interval where there are no zeros or only zeros
  expect_equal(distrib_intercept_start(d, c(1, 2, 3, 4))$zi, 1 / 8)
  expect_equal(distrib_intercept_start(d, c(0, 0, 0, 0))$zi, 1 - 1 / 8)
  expect_identical(distrib_intercept_start(d, numeric(0)), list())
})
