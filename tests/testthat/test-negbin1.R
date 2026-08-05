# NB1: the negative binomial whose variance is LINEAR in the mean. It is a
# different family from negbin2_distrib(), not a reparametrization: the size is
# mu/theta, so the mean sits inside the gamma functions rather than outside
# them, and the two have different variance functions at every (mu, theta).

test_that("the variance is linear in the mean, unlike the quadratic form", {
  d <- negbin1_distrib()
  d2 <- negbin2_distrib()
  for (th in list(list(mu = 4, theta = 4), list(mu = 2, theta = 0.5),
                  list(mu = 9, theta = 3))) {
    expect_equal(mean(d, th), th$mu, tolerance = 1e-6)
    expect_equal(variance(d, th), th$mu * (1 + th$theta), tolerance = 1e-5)
    expect_equal(variance(d2, th), th$mu + th$mu^2 / th$theta, tolerance = 1e-5)
  }

  # the ratio to the Poisson variance is constant here and grows there, which
  # is the whole difference between the two families
  r1 <- vapply(c(2, 8, 32), function(m)
    variance(d, list(mu = m, theta = 3)) / m, numeric(1))
  r2 <- vapply(c(2, 8, 32), function(m)
    variance(d2, list(mu = m, theta = 3)) / m, numeric(1))
  expect_equal(diff(range(r1)), 0, tolerance = 1e-4)
  expect_gt(diff(range(r2)), 5)
})


test_that("the mass is the negative binomial at size mu/theta", {
  d <- negbin1_distrib()
  th <- list(mu = 4, theta = 2)
  k <- 0:30
  expect_equal(distrib_pdf(d, k, th),
               stats::dnbinom(k, size = th$mu / th$theta,
                              prob = 1 / (1 + th$theta)))
  expect_equal(sum(distrib_pdf(d, 0:400, th)), 1, tolerance = 1e-10)
})


test_that("the analytical derivatives match one Richardson differentiation", {
  skip_if_not_installed("numDeriv")
  d <- negbin1_distrib()
  th <- list(mu = 4, theta = 2)
  q0 <- c(4, 2)
  at <- function(q) list(mu = q[1], theta = q[2])
  set.seed(1)
  y <- distrib_rng(d, 50, th)

  g <- distrib_gradient(d, y, th)
  expect_equal(vapply(g, sum, numeric(1)),
    numDeriv::grad(function(q) sum(distrib_pdf(d, y, at(q), log = TRUE)), q0),
    tolerance = 1e-6, ignore_attr = TRUE)

  J <- numDeriv::jacobian(function(q) {
    vapply(distrib_gradient(d, y, at(q)), sum, numeric(1))
  }, q0)
  h <- distrib_hessian(d, y, th)
  expect_equal(sum(h$mu_mu), J[1, 1], tolerance = 1e-6)
  expect_equal(sum(h$mu_theta), J[1, 2], tolerance = 1e-6)
  expect_equal(sum(h$theta_theta), J[2, 2], tolerance = 1e-6)
})


test_that("the expected information agrees with a Monte Carlo average", {
  # The series against the exact mass is checked against a route that shares
  # nothing with it: the observed Hessian averaged over draws.
  d <- negbin1_distrib()
  th <- list(mu = 4, theta = 2)
  eh <- distrib_expected_hessian(d, 0, th)
  set.seed(2)
  ys <- distrib_rng(d, 4e5, th)
  hs <- distrib_hessian(d, ys, th)
  expect_equal(vapply(eh, function(v) v[1], numeric(1)),
               vapply(hs, mean, numeric(1)), tolerance = 5e-3)
})


test_that("it approaches the Poisson as the dispersion vanishes", {
  d <- negbin1_distrib()
  dp <- poisson_distrib()
  k <- 0:20
  g1 <- max(abs(distrib_pdf(d, k, list(mu = 4, theta = 1e-3)) -
                distrib_pdf(dp, k, list(mu = 4))))
  g2 <- max(abs(distrib_pdf(d, k, list(mu = 4, theta = 1e-5)) -
                distrib_pdf(dp, k, list(mu = 4))))
  expect_lt(g2, g1 / 50)
})


test_that("the validator passes and a fit recovers the parameters", {
  d <- negbin1_distrib()
  set.seed(3)
  res <- check_distrib(d, verbose = FALSE)
  expect_true(all(res$status == "OK"),
    label = paste(res$check[res$status != "OK"], collapse = ", "))

  set.seed(5)
  th <- list(mu = 4, theta = 2)
  y <- distrib_rng(d, 5000, th)
  f <- fit_distrib(d, y)
  expect_true(f@converged, info = fit_report(f, d, th))
  expect_equal(unname(coef(f)), c(4, 2), tolerance = 0.25)
})
