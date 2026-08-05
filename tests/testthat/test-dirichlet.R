# The Dirichlet: the first multivariate family here that is not elliptical.
# There is no location and scale to separate, the support is a simplex, and
# the covariance is singular by construction.

test_that("the density is the formula and the parameters are flattened", {
  d <- dirichlet_distrib(4)
  expect_identical(d@params,
                   c("mean_alr1", "mean_alr2", "mean_alr3", "phi"))
  expect_identical(d@n_dim, 4L)

  th <- as.list(stats::setNames(c(0.4, -0.3, 0.8, 12), d@params))
  mu <- mv_location(d, th)
  a <- 12 * mu
  set.seed(1)
  y <- distrib_rng(d, 20, th)

  expect_equal(rowSums(y), rep(1, 20))
  expect_equal(distrib_pdf(d, y, th, log = TRUE),
               lgamma(12) - sum(lgamma(a)) + as.numeric(log(y) %*% (a - 1)))

  # a point off the simplex has no mass
  expect_identical(distrib_pdf(d, matrix(c(0.5, 0.5, 0.5, 0.5), 1), th), 0)
})


test_that("the moments are the ones the family promises", {
  d <- dirichlet_distrib(4)
  th <- as.list(stats::setNames(c(0.4, -0.3, 0.8, 12), d@params))
  mu <- mv_location(d, th)

  expect_equal(sum(mu), 1)
  # Cov = (diag(mu) - mu mu') / (phi + 1), written out here
  expect_equal(mv_sigma(d, th),
               (diag(mu) - tcrossprod(mu)) / 13, tolerance = 1e-12)
  # singular by construction: the coordinates sum to a constant
  expect_lt(abs(det(mv_sigma(d, th))), 1e-18)

  set.seed(2)
  big <- distrib_rng(d, 2e5, th)
  expect_equal(colMeans(big), mu, tolerance = 0.01)
  # absolute, because the off-diagonal entries are two orders of magnitude
  # smaller than the diagonal and a relative tolerance would be set by them
  S <- stats::cov(big) * (nrow(big) - 1) / nrow(big)
  expect_lt(max(abs(S - mv_sigma(d, th))), 1e-4)
})


test_that("the analytical derivatives match one Richardson differentiation", {
  skip_if_not_installed("numDeriv")
  d <- dirichlet_distrib(4)
  th <- as.list(stats::setNames(c(0.4, -0.3, 0.8, 12), d@params))
  q0 <- unlist(th)
  set.seed(1)
  y <- distrib_rng(d, 30, th)
  llv <- function(q) {
    sum(distrib_pdf(d, y, as.list(stats::setNames(q, d@params)), log = TRUE))
  }

  g <- distrib_gradient(d, y, th)
  expect_equal(vapply(g, sum, numeric(1)), numDeriv::grad(llv, q0),
               tolerance = 1e-6, ignore_attr = TRUE)

  H <- numDeriv::hessian(llv, q0)
  h <- distrib_hessian(d, y, th)
  nm <- hess_names(d@params)
  pr <- hess_pairs(d@params)
  pos <- stats::setNames(seq_along(d@params), d@params)
  for (m in seq_along(nm)) {
    expect_equal(sum(h[[m]]), H[pos[[pr[[m]][1]]], pos[[pr[[m]][2]]]],
                 tolerance = 1e-5, label = nm[m])
  }
})


test_that("the expected information is closed form", {
  # Two zero sums make it so, both from differentiating sum(mu) = 1: the
  # columns of A and every second-derivative vector of the simplex. Checked
  # against a Monte Carlo mean of the observed Hessian, which shares none of
  # that reasoning.
  d <- dirichlet_distrib(4)
  th <- as.list(stats::setNames(c(0.4, -0.3, 0.8, 12), d@params))
  set.seed(2)
  y <- distrib_rng(d, 2e5, th)
  eh <- distrib_expected_hessian(d, y[1, , drop = FALSE], th)
  hb <- distrib_hessian(d, y, th)
  # absolute against the size of the largest component: the mixed mean-phi
  # entries are two orders of magnitude smaller than the rest, and a relative
  # tolerance on them would be a tolerance on the Monte Carlo error alone
  scl <- max(abs(vapply(eh, function(z) z[1], numeric(1))))
  for (nm in names(eh)) {
    expect_lt(abs(eh[[nm]][1] - mean(hb[[nm]])), 5e-3 * scl, label = nm)
  }

  # the pure-phi component is free of the data, the family being an
  # exponential family in log y
  h1 <- distrib_hessian(d, y[1:5, ], th)
  expect_equal(length(unique(h1$phi_phi)), 1L)
})


test_that("the marginals are Beta with the same concentration", {
  d <- dirichlet_distrib(4)
  th <- as.list(stats::setNames(c(0.4, -0.3, 0.8, 12), d@params))
  mu <- mv_location(d, th)

  for (j in 1:4) {
    m <- mv_marginal(d, th, which = j)
    expect_true(S7::S7_inherits(m$distrib, BetaDistrib))
    expect_equal(m$theta$mu, mu[j])
    expect_equal(m$theta$phi, 12)
    # against dbeta with the shapes written out
    x <- c(0.05, 0.3, 0.7, 0.95)
    expect_equal(distrib_pdf(m$distrib, x, m$theta),
                 stats::dbeta(x, 12 * mu[j], 12 - 12 * mu[j]),
                 label = as.character(j))
  }

  # several coordinates are again Dirichlet only after the rest is collapsed,
  # which is a different object, so the request is refused
  expect_error(mv_marginal(d, th, which = 1:2), "one coordinate at a time")
})


test_that("a fit recovers the mean and the concentration", {
  d <- dirichlet_distrib(4)
  th <- as.list(stats::setNames(c(0.4, -0.3, 0.8, 12), d@params))
  set.seed(5)
  y <- distrib_rng(d, 3000, th)
  f <- fit_distrib(d, y)
  expect_true(f@converged)
  expect_equal(unname(coef(f)), unname(unlist(th)), tolerance = 0.1)
})


test_that("the validator passes, and still catches a broken density", {
  # The normalisation cannot be checked against a gaussian proposal here: the
  # simplex has no volume in R^p. The family supplies a uniform one instead,
  # and the injection confirms the substitution has not made the check vacuous.
  d <- dirichlet_distrib(4)
  th <- as.list(stats::setNames(c(0.4, -0.3, 0.8, 12), d@params))
  set.seed(1)
  res <- check_distrib(d, th, verbose = FALSE)
  expect_true(all(res$status == "OK"),
              info = paste(res$check[res$status != "OK"], collapse = ", "))
  # a discrete family's checks are not run on a continuous one, and the
  # response derivatives are not run on a family the base class refuses them for
  expect_false("response derivatives vs finite differences" %in% res$check)
  expect_false(has_mv_support(d))
  expect_false(has_mv_grad_y(d))

  good <- S7::method(distrib_pdf, DirichletDistrib)
  on.exit(S7::method(distrib_pdf, DirichletDistrib) <- good, add = TRUE)
  suppressMessages(
    S7::method(distrib_pdf, DirichletDistrib) <- function(distrib, y, theta, log = FALSE) {
      out <- good(distrib, y, theta, log = TRUE) + log(1.05)
      if (log) out else exp(out)
    }
  )
  set.seed(1)
  bad <- check_distrib(d, th, verbose = FALSE)
  expect_equal(bad$status[bad$check == "density integrates to 1"], "FAIL")
})


test_that("the proposal is uniform on the simplex", {
  d <- dirichlet_distrib(4)
  th <- as.list(stats::setNames(c(0.4, -0.3, 0.8, 12), d@params))
  set.seed(3)
  pr <- mv_reference_draw(d, th, 500)
  expect_equal(rowSums(pr$y), rep(1, 500))
  # constant, and equal to the density of Dirichlet(1, ..., 1)
  expect_equal(pr$logd, rep(lgamma(4), 500))
  # The base class's gaussian proposal does not fail loudly on the singular
  # covariance -- chol() accepts a matrix whose smallest eigenvalue is zero to
  # rounding -- it silently draws points off the simplex, and the estimate of
  # an integral that is one comes back around 1e-8.
  set.seed(3)
  base <- S7::method(mv_reference_draw, multivariate_distrib)(d, th, 2000)
  lf <- distrib_pdf(d, base$y, th, log = TRUE)
  expect_lt(mean(exp(lf - base$logd)), 1e-4)
})


test_that("the density is zero off the simplex, without complaining", {
  d <- dirichlet_distrib(4)
  th <- as.list(stats::setNames(c(0.4, -0.3, 0.8, 12), d@params))
  off <- rbind(c(0.25, 0.25, 0.25, 0.25), c(-0.1, 0.4, 0.4, 0.3),
               c(0.5, 0.5, 0.5, 0.5))
  expect_silent(f <- distrib_pdf(d, off, th))
  expect_identical(f, c(distrib_pdf(d, off[1, , drop = FALSE], th), 0, 0))
})


test_that("the constructor validates its dimension and its parameter", {
  expect_error(dirichlet_distrib(1), "at least 2")
  expect_error(dirichlet_distrib(3, mean = "nonsense"), "parameters7 parameter")
  expect_error(dirichlet_distrib(3, mean = parameters7::simplex(4)),
               "coordinates")
})
