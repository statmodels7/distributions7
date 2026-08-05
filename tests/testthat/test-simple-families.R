# The exponential, geometric and chi-squared, all written in their mean.
#
# Two of the three coincide with a wrapper over a larger family, and that
# comparison is the point rather than a redundancy: two independent
# implementations of one object agree exactly or one of them is wrong, so
# there is no tolerance to choose. It also exercises fixed() against something
# that was derived rather than delegated.

test_that("the three families are what their densities say", {
  mu <- 2.5
  y <- c(0.4, 1, 2.5, 7)

  expect_equal(distrib_pdf(exponential_distrib(), y, list(mu = mu)),
               stats::dexp(y, rate = 1 / mu))
  k <- 0:6
  expect_equal(distrib_pdf(geometric_distrib(), k, list(mu = mu)),
               stats::dgeom(k, prob = 1 / (1 + mu)))
  expect_equal(distrib_pdf(chisq_distrib(), y, list(mu = 4)),
               stats::dchisq(y, df = 4))

  # the parametrization is the mean, so the first moment is the parameter
  for (d in list(exponential_distrib(), geometric_distrib(), chisq_distrib())) {
    expect_equal(mean(d, list(mu = 3)), 3, tolerance = 1e-6,
                 label = d@distrib_name)
  }
  expect_equal(variance(exponential_distrib(), list(mu = 3)), 9,
               tolerance = 1e-6)
  expect_equal(variance(geometric_distrib(), list(mu = 3)), 12,
               tolerance = 1e-6)
  expect_equal(variance(chisq_distrib(), list(mu = 3)), 6, tolerance = 1e-6)
})


test_that("every order matches one Richardson differentiation of the order below", {
  skip_if_not_installed("numDeriv")
  cases <- list(
    list(d = exponential_distrib(), y = stats::rexp(60, 1 / 2.5), mu = 2.5),
    list(d = geometric_distrib(), y = stats::rgeom(60, 1 / 3.5), mu = 2.5),
    list(d = chisq_distrib(), y = stats::rchisq(60, 4), mu = 4)
  )
  set.seed(1)
  for (cs in cases) {
    d <- cs$d
    at <- function(order, m) {
      th <- list(mu = m)
      sum(switch(order,
        distrib_pdf(d, cs$y, th, log = TRUE),
        distrib_gradient(d, cs$y, th)[[1]],
        distrib_hessian(d, cs$y, th)[[1]],
        distrib_deriv3(d, cs$y, th)[[1]],
        distrib_deriv4(d, cs$y, th)[[1]]))
    }
    for (k in 1:4) {
      # never nested: the reference differentiates the ANALYTIC order below
      ref <- numDeriv::grad(function(m) at(k, m), cs$mu)
      expect_equal(at(k + 1L, cs$mu), ref, tolerance = 1e-7,
                   label = sprintf("%s order %d", d@distrib_name, k))
    }
  }
})


test_that("the expected orders are the observed ones averaged", {
  set.seed(7)
  for (cs in list(list(d = exponential_distrib(), mu = 2.5),
                  list(d = geometric_distrib(), mu = 2.5))) {
    d <- cs$d
    th <- list(mu = cs$mu)
    ysim <- distrib_rng(d, 4e5, th)
    got <- c(distrib_expected_hessian(d, cs$mu, th)[[1]][1],
             distrib_deriv3(d, cs$mu, th, expected = TRUE)[[1]][1],
             distrib_deriv4(d, cs$mu, th, expected = TRUE)[[1]][1])
    mc <- c(mean(distrib_hessian(d, ysim, th)[[1]]),
            mean(distrib_deriv3(d, ysim, th)[[1]]),
            mean(distrib_deriv4(d, ysim, th)[[1]]))
    # Monte Carlo over 4e5 draws, so the agreement is to its own noise
    expect_equal(got, mc, tolerance = 5e-3, label = d@distrib_name)
  }

  # and the first Bartlett identity, which the expected first order is
  for (d in list(exponential_distrib(), geometric_distrib(), chisq_distrib())) {
    th <- list(mu = 2)
    set.seed(8)
    ysim <- distrib_rng(d, 2e5, th)
    expect_equal(mean(distrib_gradient(d, ysim, th)[[1]]), 0,
                 tolerance = 5e-2, label = d@distrib_name)
  }
})


test_that("the chi-squared observed and expected derivatives coincide exactly", {
  # The family is a one-parameter exponential family in log y, so from the
  # second order on the derivatives do not involve the response: there is
  # nothing to average, and the agreement is exact rather than close.
  d <- chisq_distrib()
  set.seed(5)
  y <- stats::rchisq(30, 6)
  th <- list(mu = 6)

  expect_identical(distrib_hessian(d, y, th)[[1]],
                   distrib_expected_hessian(d, y, th)[[1]])
  expect_identical(distrib_deriv3(d, y, th)[[1]],
                   distrib_deriv3(d, y, th, expected = TRUE)[[1]])
  expect_identical(distrib_deriv4(d, y, th)[[1]],
                   distrib_deriv4(d, y, th, expected = TRUE)[[1]])

  # every entry is the same number, the response having left the expression
  expect_length(unique(distrib_hessian(d, y, th)[[1]]), 1L)

  # On the LINK scale they part company, by exactly the term the order-two
  # chain rule adds: h''(eta) times the score, which has mean zero but is not
  # zero for a sample. This is why Fisher scoring and Newton differ here.
  lk <- d@link_params$mu
  eta <- linkfunctions7::linkfun(lk, th$mu)
  h2 <- linkfunctions7::d2linkinv(lk, eta)
  gap <- distrib_hessian(d, y, th, scale = "link")[[1]] -
    distrib_expected_hessian(d, y, th, scale = "link")[[1]]
  expect_equal(gap, h2 * distrib_gradient(d, y, th)[[1]], tolerance = 1e-12)
  expect_gt(max(abs(gap)), 1)
})


test_that("the exponential is the unit-shape Weibull, to machine precision", {
  # Two implementations of one object: exponential_distrib() derives its own
  # kernels, fixed() splices a value into the Weibull's. Nothing here needs a
  # tolerance argument.
  a <- exponential_distrib()
  b <- fixed(weibull1_distrib(), sigma = 1)
  set.seed(2)
  y <- stats::rexp(40, 1 / 2.5)
  th <- list(mu = 2.5)

  expect_equal(distrib_pdf(a, y, th, log = TRUE), distrib_pdf(b, y, th, log = TRUE))
  expect_equal(distrib_cdf(a, y, th), distrib_cdf(b, y, th))
  for (sc in c("parameter", "link")) {
    expect_equal(distrib_gradient(a, y, th, scale = sc),
                 distrib_gradient(b, y, th, scale = sc))
    expect_equal(distrib_hessian(a, y, th, scale = sc),
                 distrib_hessian(b, y, th, scale = sc))
    expect_equal(distrib_expected_hessian(a, y, th, scale = sc),
                 distrib_expected_hessian(b, y, th, scale = sc))
  }
  for (ex in c(FALSE, TRUE)) {
    expect_equal(distrib_deriv3(a, y, th, expected = ex),
                 distrib_deriv3(b, y, th, expected = ex))
    expect_equal(distrib_deriv4(a, y, th, expected = ex),
                 distrib_deriv4(b, y, th, expected = ex))
  }
})


test_that("the geometric is the negative binomial at theta one, to machine precision", {
  a <- geometric_distrib()
  b <- fixed(negbin2_distrib(), theta = 1)
  set.seed(3)
  y <- stats::rgeom(40, 1 / 3.5)
  th <- list(mu = 2.5)

  expect_equal(distrib_pdf(a, y, th, log = TRUE), distrib_pdf(b, y, th, log = TRUE))
  expect_equal(distrib_cdf(a, y, th), distrib_cdf(b, y, th))
  for (sc in c("parameter", "link")) {
    expect_equal(distrib_gradient(a, y, th, scale = sc),
                 distrib_gradient(b, y, th, scale = sc))
    expect_equal(distrib_hessian(a, y, th, scale = sc),
                 distrib_hessian(b, y, th, scale = sc))
    expect_equal(distrib_expected_hessian(a, y, th, scale = sc),
                 distrib_expected_hessian(b, y, th, scale = sc))
  }
  for (ex in c(FALSE, TRUE)) {
    expect_equal(distrib_deriv3(a, y, th, expected = ex),
                 distrib_deriv3(b, y, th, expected = ex))
    expect_equal(distrib_deriv4(a, y, th, expected = ex),
                 distrib_deriv4(b, y, th, expected = ex))
  }
})


test_that("neither is a Gamma with a parameter held fixed", {
  # The package writes the Gamma in (mu, sigma2), whose shape is mu^2/sigma2,
  # so unit shape is the RELATION sigma2 = mu^2 and not a value sigma2 can be
  # held at: a fixed() agrees with the exponential at one mean and nowhere
  # else. The same argument rules the chi-squared out, where it is
  # sigma2 = 2*mu.
  gf <- fixed(gamma2_distrib(), sigma2 = 2.5^2)
  y <- c(0.5, 1, 3)
  expect_equal(distrib_pdf(gf, y, list(mu = 2.5)),
               stats::dexp(y, rate = 1 / 2.5))
  expect_gt(max(abs(distrib_pdf(gf, y, list(mu = 4)) -
                    stats::dexp(y, rate = 1 / 4))), 0.1)

  cf <- fixed(gamma2_distrib(), sigma2 = 2 * 4)
  expect_equal(distrib_pdf(cf, y, list(mu = 4)), stats::dchisq(y, df = 4))
  expect_gt(max(abs(distrib_pdf(cf, y, list(mu = 7)) -
                    stats::dchisq(y, df = 7))), 0.01)
})


test_that("the validator passes and the fits reproduce the closed-form MLEs", {
  for (d in list(exponential_distrib(), geometric_distrib(), chisq_distrib())) {
    set.seed(11)
    res <- check_distrib(d, verbose = FALSE)
    expect_true(all(res$status == "OK"),
      label = sprintf("%s: %s", d@distrib_name,
                      paste(res$check[res$status != "OK"], collapse = ", ")))
  }

  # for both one-parameter mean families the MLE of the mean IS the sample
  # mean, so the fit is checked against arithmetic rather than against itself
  set.seed(3)
  ye <- stats::rexp(2000, 1 / 2.5)
  expect_equal(unname(coef(fit_distrib(exponential_distrib(), ye))), mean(ye),
               tolerance = 1e-6)
  yg <- stats::rgeom(2000, 1 / 3.5)
  expect_equal(unname(coef(fit_distrib(geometric_distrib(), yg))), mean(yg),
               tolerance = 1e-6)

  yc <- stats::rchisq(2000, 5)
  fc <- fit_distrib(chisq_distrib(), yc)
  expect_true(fc@converged)
  expect_equal(unname(coef(fc)), 5, tolerance = 0.2)
})
