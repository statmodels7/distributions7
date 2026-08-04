# The beta-binomial: a binomial whose success probability is drawn from a
# Beta. Not reachable from anything else here -- the mixing is over the
# probability rather than over the outcome, so no wrapper produces it.

test_that("the mass function is the one written out by hand", {
  n <- 10
  d <- betabinom_distrib(size = n)
  th <- list(mu = 0.3, sigma = 0.5)
  k <- 0:n
  a <- th$mu / th$sigma
  b <- (1 - th$mu) / th$sigma

  expect_equal(distrib_pdf(d, k, th),
               choose(n, k) * beta(k + a, n - k + b) / beta(a, b))
  expect_equal(sum(distrib_pdf(d, k, th)), 1)

  # outside the support, and off the lattice
  expect_identical(distrib_pdf(d, c(-1, n + 1), th), c(0, 0))
  expect_identical(distrib_pdf(d, 2.5, th), 0)
})


test_that("the mean and the variance are the ones the family promises", {
  n <- 12
  d <- betabinom_distrib(size = n)
  for (th in list(list(mu = 0.3, sigma = 0.5), list(mu = 0.7, sigma = 0.1))) {
    k <- 0:n
    p <- distrib_pdf(d, k, th)
    expect_equal(sum(k * p), n * th$mu, tolerance = 1e-10)
    expect_equal(sum((k - n * th$mu)^2 * p),
                 n * th$mu * (1 - th$mu) *
                   (1 + (n - 1) * th$sigma / (1 + th$sigma)),
                 tolerance = 1e-10)
  }
})


test_that("it is overdispersed and reaches the binomial as sigma vanishes", {
  n <- 10
  d <- betabinom_distrib(size = n)
  db <- binomial_distrib(size = n)
  k <- 0:n

  expect_gt(variance(d, list(mu = 0.3, sigma = 0.5)),
            variance(db, list(mu = 0.3)))

  # the gap closes linearly in sigma, so a tenfold decrease is a tenfold gap
  g1 <- max(abs(distrib_pdf(d, k, list(mu = .3, sigma = 1e-4)) -
                distrib_pdf(db, k, list(mu = .3))))
  g2 <- max(abs(distrib_pdf(d, k, list(mu = .3, sigma = 1e-6)) -
                distrib_pdf(db, k, list(mu = .3))))
  expect_lt(g2, g1 / 50)
})


test_that("the analytical derivatives match one Richardson differentiation", {
  skip_if_not_installed("numDeriv")
  d <- betabinom_distrib(size = 10)
  th <- list(mu = 0.3, sigma = 0.5)
  p0 <- c(0.3, 0.5)
  at <- function(p) list(mu = p[1], sigma = p[2])
  set.seed(1)
  y <- distrib_rng(d, 40, th)

  g <- distrib_gradient(d, y, th)
  expect_equal(vapply(g, sum, numeric(1)),
    numDeriv::grad(function(p) sum(distrib_pdf(d, y, at(p), log = TRUE)), p0),
    tolerance = 1e-6, ignore_attr = TRUE)

  J <- numDeriv::jacobian(function(p) {
    vapply(distrib_gradient(d, y, at(p)), sum, numeric(1))
  }, p0)
  h <- distrib_hessian(d, y, th)
  expect_equal(sum(h$mu_mu), J[1, 1], tolerance = 1e-6)
  expect_equal(sum(h$mu_sigma), J[1, 2], tolerance = 1e-6)
  expect_equal(sum(h$sigma_sigma), J[2, 2], tolerance = 1e-6)
})


test_that("the expected information is an exact sum over the support", {
  d <- betabinom_distrib(size = 10)
  th <- list(mu = 0.3, sigma = 0.5)

  # written out here rather than taken from the kernel: the mass times the
  # observed Hessian, summed over the finite support
  k <- 0:10
  p <- distrib_pdf(d, k, th)
  hk <- distrib_hessian(d, k, th)
  want <- vapply(hk, function(v) sum(p * v), numeric(1))

  got <- distrib_expected_hessian(d, 0, th)
  expect_equal(vapply(got, function(v) v[1], numeric(1)), want,
               tolerance = 1e-12)

  # it does not depend on the data, only on the parameters
  e2 <- distrib_expected_hessian(d, c(1, 5, 9), th)
  expect_equal(length(unique(e2$mu_mu)), 1L)
})


test_that("the validator passes and a fit recovers the parameters", {
  d <- betabinom_distrib(size = 10)
  set.seed(3)
  res <- check_distrib(d, verbose = FALSE)
  expect_true(all(res$status == "OK"),
    label = paste(res$check[res$status != "OK"], collapse = ", "))

  set.seed(5)
  y <- distrib_rng(d, 4000, list(mu = 0.3, sigma = 0.5))
  f <- fit_distrib(d, y)
  expect_true(f@converged, info = fit_report(f, d, list(mu = 0.3, sigma = 0.5)))
  expect_equal(unname(coef(f)), c(0.3, 0.5), tolerance = 0.06)
})


test_that("size is a constant of the distribution and is validated", {
  expect_error(betabinom_distrib(size = 0), "positive integer")
  expect_error(betabinom_distrib(size = 2.5), "positive integer")
  expect_error(betabinom_distrib(size = c(2, 3)), "positive integer")

  d <- betabinom_distrib(size = 7)
  expect_identical(d@size, 7)
  expect_identical(d@bounds, c(0, 7))
  expect_identical(d@params, c("mu", "sigma"))
})
