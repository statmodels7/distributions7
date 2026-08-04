# Folding at zero. The map y -> |y| is two to one, so this is not a change of
# variable and cannot be a transformer(): the density adds the two preimages
# instead of carrying one through a Jacobian.

test_that("the density and the distribution function add the two preimages", {
  d <- folded(gaussian_distrib())
  th <- list(mu = 0.5, sigma = 1)
  x <- c(0, 0.2, 0.8, 2, 3.5)

  expect_equal(distrib_pdf(d, x, th),
               stats::dnorm(x, 0.5, 1) + stats::dnorm(-x, 0.5, 1))
  expect_equal(distrib_cdf(d, x, th),
               stats::pnorm(x, 0.5, 1) - stats::pnorm(-x, 0.5, 1))

  # nothing below zero, and the whole mass above it
  expect_identical(distrib_pdf(d, c(-1, -0.1), th), c(0, 0))
  expect_equal(stats::integrate(function(z) distrib_pdf(d, z, th), 0, Inf)$value,
               1, tolerance = 1e-8)

  # the support is the folded one and the parameters are the parent's
  expect_identical(d@bounds, c(0, Inf))
  expect_identical(d@params, gaussian_distrib()@params)
  expect_identical(d@n_params, gaussian_distrib()@n_params)
})


test_that("drawing is the absolute value of the parent's draw", {
  d <- folded(gaussian_distrib())
  th <- list(mu = 0.5, sigma = 1)
  set.seed(21)
  a <- distrib_rng(d, 50, th)
  set.seed(21)
  b <- abs(distrib_rng(gaussian_distrib(), 50, th))
  expect_identical(a, b)
  expect_true(all(a >= 0))
})


test_that("every order matches one Richardson differentiation of the order below", {
  skip_if_not_installed("numDeriv")
  d <- folded(gaussian_distrib())
  th <- list(mu = 0.5, sigma = 1)
  x <- c(0.2, 0.8, 2, 3.5)
  p0 <- c(0.5, 1)
  at <- function(p) list(mu = p[1], sigma = p[2])

  # order 1 against the log-density
  g <- distrib_gradient(d, x, th)
  expect_equal(vapply(g, sum, numeric(1)),
    numDeriv::grad(function(p) sum(distrib_pdf(d, x, at(p), log = TRUE)), p0),
    tolerance = 1e-7, ignore_attr = TRUE)

  # order 2 against the analytic gradient
  h <- distrib_hessian(d, x, th)
  J <- numDeriv::jacobian(function(p) {
    vapply(distrib_gradient(d, x, at(p)), sum, numeric(1))
  }, p0)
  expect_equal(sum(h$mu_mu), J[1, 1], tolerance = 1e-7)
  expect_equal(sum(h$mu_sigma), J[1, 2], tolerance = 1e-7)
  expect_equal(sum(h$sigma_sigma), J[2, 2], tolerance = 1e-7)

  # orders 3 and 4 against the analytic order below, never nested
  h_nm <- names(distrib_hessian(d, x, th))
  d3 <- distrib_deriv3(d, x, th)
  J3 <- numDeriv::jacobian(function(p) {
    vapply(distrib_hessian(d, x, at(p)), sum, numeric(1))
  }, p0)
  for (nm in names(d3)) {
    part <- strsplit(nm, "_")[[1]]
    k <- match(paste(part[1:2], collapse = "_"), h_nm)
    if (is.na(k)) k <- match(paste(rev(part[1:2]), collapse = "_"), h_nm)
    j <- match(part[3], c("mu", "sigma"))
    expect_equal(sum(d3[[nm]]), J3[k, j], tolerance = 1e-6, label = nm)
  }

  d4 <- distrib_deriv4(d, x, th)
  J4 <- numDeriv::jacobian(function(p) {
    vapply(distrib_deriv3(d, x, at(p)), sum, numeric(1))
  }, p0)
  for (nm in names(d4)) {
    part <- strsplit(nm, "_")[[1]]
    k <- match(paste(part[1:3], collapse = "_"), names(d3))
    j <- match(part[4], c("mu", "sigma"))
    if (!is.na(k)) {
      expect_equal(sum(d4[[nm]]), J4[k, j], tolerance = 1e-5, label = nm)
    }
  }
})


test_that("the score is the mixture score its derivation says it is", {
  # w s(x) + (1 - w) s(-x), the two preimages weighted by which one the point
  # came from. Written out here rather than taken from the package.
  d <- folded(gaussian_distrib())
  th <- list(mu = 0.5, sigma = 1)
  x <- c(0.2, 0.8, 2)

  fp <- stats::dnorm(x, 0.5, 1)
  fm <- stats::dnorm(-x, 0.5, 1)
  w <- fp / (fp + fm)
  sp <- distrib_gradient(gaussian_distrib(), x, th)
  sm <- distrib_gradient(gaussian_distrib(), -x, th)

  g <- distrib_gradient(d, x, th)
  expect_equal(g$mu, w * sp$mu + (1 - w) * sm$mu)
  expect_equal(g$sigma, w * sp$sigma + (1 - w) * sm$sigma)
})


test_that("the response derivatives carry the reflection's sign", {
  skip_if_not_installed("numDeriv")
  d <- folded(gaussian_distrib())
  th <- list(mu = 0.5, sigma = 1)
  x <- c(0.2, 0.8, 2, 3.5)

  expect_equal(distrib_grad_y(d, x, th),
    numDeriv::grad(function(z) sum(distrib_pdf(d, z, th, log = TRUE)), x),
    tolerance = 1e-8)
  expect_equal(distrib_hess_y(d, x, th),
    diag(numDeriv::jacobian(function(z) distrib_grad_y(d, z, th), x)),
    tolerance = 1e-8)
})


test_that("the half-normal is a folded gaussian with its location held at zero", {
  hn <- fixed(folded(gaussian_distrib()), mu = 0)
  expect_identical(hn@params, "sigma")

  s <- 2
  x <- c(0.2, 1, 3)
  expect_equal(distrib_pdf(hn, x, list(sigma = s)),
               sqrt(2 / pi) / s * exp(-x^2 / (2 * s^2)))
  # the mean of a half-normal is sigma sqrt(2/pi)
  expect_equal(mean(hn, list(sigma = s)), s * sqrt(2 / pi), tolerance = 1e-6)
  expect_equal(variance(hn, list(sigma = s)), s^2 * (1 - 2 / pi),
               tolerance = 1e-6)
})


test_that("a symmetric parent folds to twice its upper half", {
  # With mu = 0 the two preimages carry the same density, so w is one half
  # everywhere and the folded density is exactly doubled.
  d <- folded(gaussian_distrib())
  th <- list(mu = 0, sigma = 1.5)
  x <- c(0.3, 1, 2.5)
  expect_equal(distrib_pdf(d, x, th), 2 * stats::dnorm(x, 0, 1.5))

  # and the score in mu vanishes, the two preimages cancelling
  expect_equal(distrib_gradient(d, x, th)$mu, rep(0, length(x)))
})


test_that("the validator passes on the fold and on the half-normal", {
  for (d in list(folded(gaussian_distrib()),
                 fixed(folded(gaussian_distrib()), mu = 0))) {
    set.seed(4)
    res <- check_distrib(d, verbose = FALSE)
    expect_true(all(res$status == "OK"),
      label = sprintf("%s: %s", d@distrib_name,
                      paste(res$check[res$status != "OK"], collapse = ", ")))
  }
})


test_that("folded() refuses what it would silently mishandle", {
  # a parent already on the non-negative half line folds to itself, so the
  # call is a mistake rather than a no-op
  expect_error(folded(gamma_distrib()), "leaves alone")
  expect_error(folded(folded(gaussian_distrib())), "leaves alone")

  # an atom would be counted twice at zero or moved onto its reflection
  expect_error(folded(zero_adjusted(gaussian_distrib())), "atom")

  # and folding a lattice is not what the sum above computes
  expect_error(folded(poisson_distrib()), "continuous")
})


test_that("the sign of a symmetric parent's location is not identified", {
  # f(-x; mu) = f(x; -mu) for a parent symmetric about its location, so the
  # two terms of L merely swap: the likelihood is EXACTLY even in mu. This is
  # a property of the model, not of the arithmetic, so the comparison is
  # against zero rather than against a tolerance.
  x <- c(0.2, 0.9, 2.3, 4)
  for (parent in list(gaussian_distrib(), cauchy_distrib(),
                      logistic_distrib(), laplace_distrib())) {
    d <- folded(parent)
    nm <- d@params
    tp <- stats::setNames(as.list(c(1.2, 1.5)), nm)
    tm <- stats::setNames(as.list(c(-1.2, 1.5)), nm)
    expect_identical(distrib_pdf(d, x, tp, log = TRUE),
                     distrib_pdf(d, x, tm, log = TRUE),
                     label = parent@distrib_name)
  }

  # and it is the parent's symmetry that does it, not the folding: a skew
  # parent has no such invariance
  ds <- folded(skewnormal_distrib())
  expect_gt(max(abs(
    distrib_pdf(ds, x, list(mu = 1.2, sigma = 1.5, alpha = 3), log = TRUE) -
    distrib_pdf(ds, x, list(mu = -1.2, sigma = 1.5, alpha = 3), log = TRUE))), 1)
})


test_that("a fit through the fold recovers what generated the data", {
  set.seed(31)
  d <- folded(gaussian_distrib())
  true <- list(mu = 1.2, sigma = 2)
  y <- distrib_rng(d, 3000, true)

  # The location's sign is not identified, so the two starts reach the two
  # maxima and the estimable quantity is its magnitude. Their maximised
  # log-likelihoods must agree, which is what makes them the same maximum.
  fp <- fit_distrib(d, y, start = list(mu = 3, sigma = 1))
  fm <- fit_distrib(d, y, start = list(mu = -3, sigma = 1))
  expect_true(fp@converged, info = fit_report(fp, d, y))
  expect_true(fm@converged, info = fit_report(fm, d, y))

  expect_equal(as.numeric(logLik(fp)), as.numeric(logLik(fm)), tolerance = 1e-8)
  expect_equal(unname(coef(fp)[1]), -unname(coef(fm)[1]), tolerance = 1e-4)
  expect_equal(unname(coef(fp)[2]), unname(coef(fm)[2]), tolerance = 1e-6)

  expect_equal(abs(unname(coef(fp)[1])), 1.2, tolerance = 0.15)
  expect_equal(unname(coef(fp)[2]), 2, tolerance = 0.15)

  # holding the location at zero removes the question entirely
  set.seed(32)
  hn <- fixed(folded(gaussian_distrib()), mu = 0)
  yh <- distrib_rng(hn, 3000, list(sigma = 2))
  fh <- fit_distrib(hn, yh)
  expect_true(fh@converged, info = fit_report(fh, hn, yh))
  expect_equal(unname(coef(fh)), 2, tolerance = 0.1)
})
