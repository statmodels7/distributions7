# Components evaluated far outside the region a fit normally visits. The
# ordinary tests keep to sensible parameters, which is where a formula is
# easiest to check and where these defects are invisible: a derivative written
# as a ratio whose denominator carries a power of the scale, or a mass written
# as a difference of two log-gammas, is exact in the middle and silently wrong
# at the edge. An optimizer reaches the edge whenever a start is poor, and a
# component that returns zero there is read as stationarity.

test_that("the gaussian score survives a scale no power of which is representable", {
  d <- gaussian1_distrib()
  set.seed(31)
  y <- distrib_rng(d, 200L, list(mu = 30, sigma = 7))
  n <- length(y)

  # On the link scale the score in the scale is exactly sum(z^2 - 1)/n, which
  # tends to -1 as the scale grows. Written as (res^2 - sigma^2)/sigma^3 the
  # denominator overflows first and the whole component comes back as zero.
  for (sg in c(1e10, 1e102, 8e102, 1e150, 1e200, 1e300)) {
    z <- (y - 14) / sg
    got <- sum(distrib_gradient(d, y, list(mu = 14, sigma = sg),
                                scale = "link")$sigma) / n
    expect_equal(got, sum(z^2 - 1) / n, tolerance = 1e-12)
    expect_true(is.finite(got))
  }

  # and the second derivative, whose old form carried sigma^4 and so failed
  # from 1.2e77 upwards
  for (sg in c(1e50, 1e100, 1e150)) {
    z <- (y - 14) / sg
    h <- distrib_hessian(d, y, list(mu = 14, sigma = sg), scale = "link")
    expect_equal(sum(h$sigma_sigma) / n, sum(-2 * z^2) / n, tolerance = 1e-10)
    expect_true(all(is.finite(unlist(h))))
  }

  # On the parameter scale the kernel is right wherever the value itself is
  # representable, which for the second derivative is the whole range.
  for (sg in c(1e200, 1e300)) {
    z <- (y - 14) / sg
    h <- distrib_hessian(d, y, list(mu = 14, sigma = sg))
    expect_equal(h$sigma_sigma, (1 - 3 * z^2) / sg^2)
    expect_true(all(is.finite(unlist(h))))
  }

  # Above 1.3e154 the LINK-SCALE second order is not reachable, and the reason
  # is the chain rule rather than the kernel: it forms (h')^2, which overflows,
  # against a component of order 1/sigma^2, which underflows. What must not
  # happen is a plausible number -- a zero there would be read as a flat
  # likelihood -- so the contract is that the result is not finite.
  hb <- distrib_hessian(d, y, list(mu = 14, sigma = 1e200), scale = "link")
  expect_true(all(is.nan(hb$sigma_sigma)))
})

test_that("the rewritten gaussian components agree with the forms they replace", {
  # The old expressions are transcribed here rather than called, so this pins
  # the new bodies against the algebra and not against themselves. They are
  # compared only where the old ones were reliable.
  d <- gaussian1_distrib()
  set.seed(4)
  y <- rnorm(50, 3, 2)
  for (s in c(0.01, 0.5, 7, 1e3, 1e30, 1e60)) {
    th <- list(mu = 3, sigma = s)
    r <- y - 3
    expect_equal(distrib_gradient(d, y, th)$sigma, (r^2 - s^2) / s^3)
    expect_equal(distrib_hessian(d, y, th)$sigma_sigma, (s^2 - 3 * r^2) / s^4)
    expect_equal(distrib_hessian(d, y, th)$mu_sigma, -2 * r / s^3)
    expect_equal(distrib_deriv3(d, y, th)$sigma_sigma_sigma,
                 -2 * (s^2 - 6 * r^2) / s^5)
    expect_equal(distrib_deriv4(d, y, th)$sigma_sigma_sigma_sigma,
                 6 * (s^2 - 10 * r^2) / s^6)
  }

  d2 <- gaussian2_distrib()
  for (v in c(0.01, 4, 1e6, 1e40)) {
    th <- list(mu = 3, sigma2 = v)
    r <- y - 3
    expect_equal(distrib_gradient(d2, y, th)$sigma2, (r^2 / v - 1) / (2 * v))
    expect_equal(distrib_hessian(d2, y, th)$sigma2_sigma2,
                 1 / (2 * v^2) - r^2 / v^3)
    expect_equal(distrib_deriv3(d2, y, th)$sigma2_sigma2_sigma2,
                 -1 / v^3 + 3 * r^2 / v^4)
    expect_equal(distrib_deriv4(d2, y, th)$sigma2_sigma2_sigma2_sigma2,
                 3 / v^4 - 12 * r^2 / v^5)
  }
})

test_that("the beta-binomial mass stays a probability at any shapes", {
  bb <- betabinom2_distrib(size = 10L)
  for (k in c(1, 1e3, 1e6, 1e12, 1e15, 1e18, 1e23)) {
    th <- list(alpha = 27.38 * k, beta = 52.68 * k)
    expect_equal(sum(distrib_pdf(bb, 0:10, th)), 1, tolerance = 1e-9)
    expect_true(all(distrib_pdf(bb, 0:10, th, log = TRUE) <= 0))
  }

  # As the shapes grow at a fixed ratio the family becomes the binomial at the
  # limiting proportion, which is an independent expression of the same number.
  p <- 27.38 / (27.38 + 52.68)
  expect_equal(distrib_pdf(bb, 0:10, list(alpha = 27.38e20, beta = 52.68e20),
                           log = TRUE),
               distrib_pdf(binomial_distrib(size = 10L), 0:10, list(mu = p),
                           log = TRUE),
               tolerance = 1e-10)

  # the same through the mean-dispersion chart, where huge shapes are a tiny
  # dispersion and the kernel is the compiled one
  b1 <- betabinom1_distrib(size = 10L)
  for (sg in c(1e-2, 1e-8, 1e-16, 1e-22)) {
    expect_equal(sum(distrib_pdf(b1, 0:10, list(mu = 0.35, sigma = sg))), 1,
                 tolerance = 1e-9)
  }
  expect_equal(distrib_pdf(b1, 0:10, list(mu = 0.35, sigma = 1e-20),
                           log = TRUE),
               distrib_pdf(binomial_distrib(size = 10L), 0:10,
                           list(mu = 0.35), log = TRUE),
               tolerance = 1e-10)

  # and the two charts agree with each other at ordinary values, so the route
  # taken at large shapes has not been bought at the cost of the ordinary one
  s <- 1 / 80.06
  expect_equal(distrib_pdf(b1, 0:10, list(mu = 0.35, sigma = s), log = TRUE),
               distrib_pdf(bb, 0:10, list(alpha = 0.35 / s, beta = 0.65 / s),
                           log = TRUE),
               tolerance = 1e-10)
})

test_that("a beta-binomial fit lands on the optimum rather than on huge shapes", {
  # Binomial data ask the beta-binomial for shapes running to infinity at a
  # fixed ratio. One start in six used to arrive there and report a
  # log-likelihood of zero, which beat every real fit.
  bb <- betabinom2_distrib(size = 10L)
  set.seed(31)
  y <- rbinom(200L, 10L, 0.35)
  f <- suppressWarnings(fit_distrib(bb, y))
  expect_true(as.numeric(logLik(f)) < -300)
  expect_true(all(unlist(coef(f)) < 1e6))

  set.seed(20260810L)
  for (s in distrib_start(bb, y, 6L)) {
    fs <- tryCatch(suppressWarnings(fit_distrib(bb, y, start = s)),
                   error = function(e) NULL)
    if (!is.null(fs)) expect_true(as.numeric(logLik(fs)) <= 0)
  }
})

test_that("fit_distrib rejects a discrete run whose log-likelihood is positive", {
  # The mass of a discrete family is a probability, so its logarithm cannot be
  # positive. The defect is injected on a SUBCLASS: registering a method
  # mutates the generic in place, so a broken mass on a shipped family would
  # stay broken for the rest of the session.
  Impossible <- S7::new_class("ImpossibleMass", parent = discrete_distrib,
                              package = NULL)
  S7::method(distrib_pdf, Impossible) <- function(distrib, y, theta,
                                                  log = FALSE) {
    # every mass slightly above one, as the beta-binomial's was at shapes
    # around 1e15 before its log-mass was rewritten
    ld <- rep(0.01, length(y))
    if (log) ld else exp(ld)
  }
  imp <- Impossible(
    distrib_name = "impossible", dimension = "univariate",
    bounds = c(0, Inf), params = "mu",
    params_interpretation = c(mu = "mean"), n_params = 1,
    params_bounds = list(mu = c(0, Inf)),
    link_params = list(mu = linkfunctions7::log_link()),
    params_smooth = c(mu = TRUE)
  )
  expect_error(suppressWarnings(fit_distrib(imp, c(1, 2, 3, 2, 1))),
               "positive log-likelihood")

  # and the guard has not blunted the ordinary case: a discrete family whose
  # mass is a probability still fits, and a log-likelihood of exactly zero is
  # left alone, being what a degenerate but legitimate fit reports
  set.seed(7)
  p <- poisson_distrib()
  fp <- suppressWarnings(fit_distrib(p, rpois(100L, 4)))
  expect_true(isTRUE(fp@converged))
  expect_true(as.numeric(logLik(fp)) < 0)
})
