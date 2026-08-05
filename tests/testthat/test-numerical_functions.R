# The expectation engines, rewired onto numericals7. References are closed
# forms and direct finite sums that share no code with the engines.

test_that("discrete_support_sum handles the four support shapes", {
  # finite: an exact one-shot sum
  got <- discrete_support_sum(function(k, i) dbinom(k, 10, 0.3), 0, 10, 1L)
  expect_equal(got, 1, tolerance = 1e-12)

  # unbounded above, with the mass far from the start (the tail-guard case)
  got <- discrete_support_sum(function(k, i) dpois(k, 250), 0, Inf, 1L)
  expect_equal(got, 1, tolerance = 1e-8)

  # unbounded below, by reflection
  got <- discrete_support_sum(function(k, i) dpois(-k, 3), -Inf, 0, 1L)
  expect_equal(got, 1, tolerance = 1e-9)

  # unbounded on both sides, folded around zero: sum of exp(-|k|)
  got <- discrete_support_sum(function(k, i) exp(-abs(k)), -Inf, Inf, 1L)
  expect_equal(got, (exp(1) + 1) / (exp(1) - 1), tolerance = 1e-9)

  # several rows retire independently
  lam <- c(0.5, 4, 60)
  got <- discrete_support_sum(function(k, i) dpois(k, lam[i]), 0, Inf, 3L)
  expect_equal(got, rep(1, 3), tolerance = 1e-9)
})


test_that("a divergent series is refused with an error, not returned", {
  expect_error(
    discrete_support_sum(function(k, i) 1 / (k + 1), 0, Inf, 1L),
    "did not converge"
  )
})


test_that("expectation calculates correct values for discrete distributions", {
  d <- poisson_distrib()
  f_mean <- function(y, theta) y

  # Expected value of Poisson(mu) is mu
  res <- expectation(d, f_mean, theta = list(mu = 5))
  expect_equal(res, 5, tolerance = 1e-5)
})

test_that("expectation calculates correct values for continuous distributions", {
  d <- gaussian1_distrib()
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

test_that("a vector theta is one batched call, continuous case", {
  d <- gamma2_distrib()
  mu <- c(0.5, 2, 7, 30)
  s2 <- c(0.2, 1, 4, 10)
  got <- expectation(d, function(y, theta) y, list(mu = mu, sigma2 = s2))
  expect_equal(got, mu, tolerance = 1e-7)

  # and a distribution whose mass sits far from the origin
  dn <- gaussian1_distrib()
  got <- expectation(dn, function(y, theta) y, list(mu = c(-300, 1000), sigma = 2))
  expect_equal(got, c(-300, 1000), tolerance = 1e-6)
})

test_that("expectation throws error on name collision", {
  d <- poisson_distrib()
  f_dummy <- function(y, theta, mu) y
  expect_error(
    expectation(d, f_dummy, theta = list(mu = 2), mu = 5),
    "cannot have the same names"
  )
})


test_that("a combination the batch refuses is rescued by one scalar integrate", {
  # gamma2 at shape mu^2/sigma2 = 0.49: the mu_mu Hessian component times
  # the density behaves like y^(shape-2+k) at zero, an integrable
  # singularity too harsh for bisection at the default depth; found by a
  # 200-row benchmark whose row 22 this is. The batch refuses, the scalar
  # rescue reaches it, and the answer must satisfy the second Bartlett
  # identity within quadrature accuracy.
  d <- gamma2_distrib()
  th <- list(mu = 1.2916372833, sigma2 = 3.3876254116)
  eh <- distrib_expected_hessian(d, 1, th)          # closed form, reference
  got <- expectation(d, function(y, theta)
    distrib_hessian(d, y, theta)[["mu_mu"]], th)
  expect_equal(got, eh[["mu_mu"]], tolerance = 1e-6)
})
