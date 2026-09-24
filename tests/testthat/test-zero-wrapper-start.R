# The zero wrappers start from the data, and fit_distrib() does not keep the
# first converged run when a later start reaches a higher maximum.

za_pig1_sample <- function(seed = 1) {
  set.seed(seed)
  distrib_rng(zero_adjusted(pig1_distrib()), 1500,
              list(mu = 2, sigma = 0.7, za = 0.5))
}

test_that("a zero-adjusted start reads the proportion of zeros and the positives", {
  y <- za_pig1_sample()
  d <- zero_adjusted(pig1_distrib())
  set.seed(1)
  st <- distrib_start(d, y, n_start = 4)
  expect_length(st, 4L)
  for (s in st) expect_identical(names(s), d@params)
  expect_equal(st[[1]]$za, mean(y == 0))
  set.seed(1)
  pst <- distrib_start(pig1_distrib(), y[y > 0], n_start = 4)
  expect_equal(st[[1]][c("mu", "sigma")], pst[[1]])
  # The start the fallback used to give read the median of the whole sample,
  # which is zero when half the values are.
  expect_gt(st[[1]]$mu, 1)
})

test_that("a continuous zero-adjusted start drops only the exact zeros", {
  set.seed(3)
  y <- rgamma(800, shape = 2, rate = 1)
  y[runif(800) < 0.3] <- 0
  d <- zero_adjusted(gamma1_distrib())
  st <- distrib_start(d, y, n_start = 2)
  expect_equal(st[[1]]$za, mean(y == 0))
  expect_equal(st[[1]][c("mu", "phi")],
               distrib_start(gamma1_distrib(), y[y != 0], n_start = 2)[[1]])
})

test_that("a zero-inflated start reads the parent on every observation", {
  y <- za_pig1_sample(2)
  d <- zero_inflated(pig1_distrib())
  st <- distrib_start(d, y, n_start = 3)
  expect_identical(names(st[[1]]), d@params)
  expect_equal(st[[1]]$zi, mean(y == 0))
  expect_equal(st[[1]][c("mu", "sigma")],
               distrib_start(pig1_distrib(), y, n_start = 3)[[1]])
})

test_that("the proportion of zeros is kept inside the open interval", {
  d <- zero_adjusted(poisson_distrib())
  expect_equal(distrib_start(d, c(1, 2, 3, 4), 1)[[1]]$za, 1 / 8)
  expect_equal(distrib_start(d, c(0, 0, 0, 5), 1)[[1]]$za, 3 / 4)
})

test_that("fit_distrib keeps the best of the converged runs, not the first", {
  y <- za_pig1_sample(1)
  d <- zero_adjusted(pig1_distrib())
  # The starts the univariate fallback produced before the wrappers had a
  # method: the first converges to a degenerate point, the others to the
  # maximum.
  set.seed(20260810L)
  old <- start_from_moments(d, y, 5L)
  lls <- vapply(old, function(s) {
    f <- suppressWarnings(fit_distrib(d, y, start = s))
    if (f@converged) as.numeric(logLik(f)) else NA_real_
  }, 0)
  expect_true(isTRUE(lls[[1]] < max(lls, na.rm = TRUE) - 50))
  local_mocked_bindings(distrib_start = function(distrib, y, n_start = 5L, ...) old)
  f <- fit_distrib(d, y)
  expect_true(f@converged)
  expect_equal(as.numeric(logLik(f)), max(lls, na.rm = TRUE), tolerance = 1e-6)
})
