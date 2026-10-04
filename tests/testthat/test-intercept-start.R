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

test_that("skewnormal1 starts at its moment estimate beyond the skewness bound", {
  d <- skewnormal1_distrib()
  # a mild skewness: the intercept-only fit is the start
  set.seed(3)
  y0 <- distrib_rng(d, n = 500, theta = list(mu = 0, sigma = 1, alpha = 1))
  expect_identical(distrib_intercept_start(d, y0), list())
  # a skewness beyond 0.9953: the moment estimate with the skewness held at 0.9
  y <- rexp(500)
  s <- distrib_intercept_start(d, y)
  expect_named(s, c("mu", "sigma", "alpha"))
  r <- (2 * 0.9 / (4 - pi))^(1 / 3)
  expect_equal(s$alpha, r / sqrt(2 / pi + (2 / pi - 1) * r^2))
  expect_equal(s$sigma, sd(y) * sqrt(1 + r^2))
  expect_equal(s$mu, mean(y) - sd(y) * r)
  # and with the sign of a left skewness
  expect_equal(distrib_intercept_start(d, -y)$alpha, -s$alpha)
  expect_identical(distrib_intercept_start(d, c(1, 2)), list())
})

test_that("pseudohuber starts at nu = 10 below gaussian kurtosis", {
  d <- pseudohuber_distrib()
  # heavy tails: the intercept-only fit is the start
  set.seed(4)
  expect_identical(distrib_intercept_start(d, rt(500, df = 4)), list())
  # a uniform sample has excess kurtosis -1.2, below the family's range
  y <- runif(500, 2, 6)
  s <- distrib_intercept_start(d, y)
  expect_named(s, c("mu", "sigma", "nu"))
  expect_equal(s$nu, 10)
  expect_equal(s$mu, mean(y))
  # the variance at the start is the sample variance
  expect_equal(variance(d, theta = s), var(y))
  expect_identical(distrib_intercept_start(d, c(1, 2, 3)), list())
})
