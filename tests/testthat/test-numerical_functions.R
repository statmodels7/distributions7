test_that("numerical_series calculates standard convergent series", {
  # Basel problem: 1/x^2 converges to pi^2 / 6
  f_basel <- function(x) 1 / x^2
  expect_equal(numerical_series(f_basel, start = 1, end = Inf), pi^2 / 6, tolerance = 1e-5)
})

test_that("numerical_series handles alternating series correctly", {
  # We use a fast-converging alternating series to prevent the test from taking too long.
  # Geometric alternating series: sum_{x=0}^{Inf} (-0.5)^x = 1 / (1 - (-0.5)) = 2/3
  f_alternating <- function(x) (-0.5)^x
  expect_equal(numerical_series(f_alternating, start = 0, end = Inf), 2 / 3, tolerance = 1e-5)
})

test_that("numerical_series protects against underflow with delayed mass", {
  # Distribution shifted far to the right
  f_delayed <- function(x) stats::dnorm(x, mean = 250, sd = 5)
  # The sum of this probability mass over the integer grid should be ~1
  expect_equal(numerical_series(f_delayed, start = 0, end = Inf, step = 100), 1, tolerance = 1e-5)
})

test_that("numerical_series handles doubly infinite domains", {
  # Exponential decay in both directions
  f_doubly_inf <- function(x) exp(-abs(x))
  expected_sum <- (exp(1) + 1) / (exp(1) - 1)
  expect_equal(numerical_series(f_doubly_inf, start = -Inf, end = Inf), expected_sum, tolerance = 1e-5)
})

test_that("numerical_series catches divergent series and exits early", {
  f_divergent <- function(x) x / 1000
  expect_warning(
    res <- numerical_series(f_divergent, start = 1, end = Inf),
    "divergent"
  )
  expect_true(is.numeric(res))
})

test_that("expectation calculates correct values for discrete distributions", {
  d <- poisson_distrib()
  f_mean <- function(y, theta) y
  
  # Expected value of Poisson(mu) is mu
  res <- expectation(d, f_mean, theta = list(mu = 5))
  expect_equal(res, 5, tolerance = 1e-5)
})

test_that("expectation calculates correct values for continuous distributions", {
  d <- gaussian_distrib()
  f_mean <- function(y, theta) y
  f_var <- function(y, theta) (y - theta$mu)^2
  
  # Expected value of Gaussian is mu, variance is sigma^2
  res_mean <- expectation(d, f_mean, theta = list(mu = 10, sigma = 2))
  res_var <- expectation(d, f_var, theta = list(mu = 10, sigma = 2))
  
  expect_equal(res_mean, 10, tolerance = 1e-5)
  expect_equal(res_var, 4, tolerance = 1e-5) # sigma^2 = 2^2 = 4
})

test_that("expectation handles vectorization correctly", {
  d <- poisson_distrib()
  f_pow <- function(y, theta, gamma = 1) y^gamma
  
  # Vectorized over extra arguments (gamma)
  # For Poisson(2): E[Y] = 2, E[Y^2] = Var(Y) + E[Y]^2 = 2 + 4 = 6
  res_gamma <- expectation(d, f_pow, theta = list(mu = 2), gamma = c(1, 2))
  expect_equal(res_gamma, c(2, 6), tolerance = 1e-5)
  
  # Vectorized over distribution parameters (mu)
  res_mu <- expectation(d, f_pow, theta = list(mu = c(2, 3)), gamma = 1)
  expect_equal(res_mu, c(2, 3), tolerance = 1e-5)
})

test_that("expectation throws error on name collision", {
  d <- poisson_distrib()
  f_dummy <- function(y, theta, mu) y
  expect_error(
    expectation(d, f_dummy, theta = list(mu = 2), mu = 5),
    "cannot have the same names"
  )
})
