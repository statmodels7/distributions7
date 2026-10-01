test_that("constant_to_scalar collapses only a constant vector", {
  expect_identical(constant_to_scalar(rep(2, 5)), 2)
  expect_identical(constant_to_scalar(c(2, 3)), c(2, 3))
  expect_identical(constant_to_scalar(4), 4)
  expect_identical(constant_to_scalar(c(1, NA)), c(1, NA))
  expect_identical(constant_to_scalar(c(NA_real_, NA_real_)), c(NA_real_, NA_real_))
})

test_that("the negative binomial expected information does not depend on collapsing", {
  d <- negbin2_distrib()
  y <- c(0, 1, 3, 7, 2, 0)
  th <- list(mu = rep(4, 6), theta = rep(0.7, 6))
  h <- distrib_expected_hessian(d, y, th)
  # the kernel called with the full vectors, one series per observation
  k <- negbin_expected_hessian_cpp(y, th$mu, th$theta, 1L)
  for (nm in names(k)) expect_identical(h[[nm]], k[[nm]])
  # one entry per observation, whichever way it was computed
  expect_length(h$theta_theta, 6L)
  expect_true(all(h$theta_theta == h$theta_theta[1]))
})

test_that("the collapse leaves parameters that vary by observation alone", {
  d <- negbin2_distrib()
  y <- c(0, 1, 3, 7)
  th <- list(mu = c(1, 2, 3, 4), theta = rep(2, 4))
  h <- distrib_expected_hessian(d, y, th)
  k <- negbin_expected_hessian_cpp(y, th$mu, th$theta, 1L)
  for (nm in names(k)) expect_identical(h[[nm]], k[[nm]])
})

test_that("a dispersion at the edge of the space is finite and does not run for minutes", {
  d <- negbin2_distrib()
  th <- list(mu = rep(0.0257, 5), theta = rep(1e-13, 5))
  t0 <- system.time(h <- distrib_expected_hessian(d, c(0, 1, 2, 3, 0), th))[["elapsed"]]
  expect_true(all(is.finite(unlist(h))))
  expect_lt(t0, 30)
})

test_that("a small dispersion, where a series is long, gives the same values", {
  d <- negbin2_distrib()
  n <- 50L
  y <- seq_len(n) %% 5
  th <- list(mu = rep(0.05, n), theta = rep(1e-4, n))
  h <- distrib_expected_hessian(d, y, th)
  k <- negbin_expected_hessian_cpp(y, th$mu, th$theta, 1L)
  for (nm in names(k)) expect_identical(h[[nm]], k[[nm]])
})
