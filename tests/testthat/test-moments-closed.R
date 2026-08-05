# Closed-form moments. Every method registered on a family replaces a
# quadrature (or, for a discrete family, an exact sum), so the reference here
# is that same numerical route: it shares no code with the formulas.
#
# Kurtosis is the EXCESS throughout, so a gaussian is 0 and a Laplace is 3.

num_moments <- function(d, th) {
  m1 <- moment(d, th, p = 1, central = FALSE)
  c2 <- moment(d, th, p = 2, central = TRUE, mu = m1)
  c3 <- moment(d, th, p = 3, central = TRUE, mu = m1)
  c4 <- moment(d, th, p = 4, central = TRUE, mu = m1)
  c(mean = m1, var = c2, skew = c3 / c2^1.5, kurt = c4 / c2^2 - 3)
}

moment_cases <- function() {
  list(
    list("gaussian1",   gaussian1_distrib(),           list(mu = 2, sigma = 3)),
    list("logistic",    logistic_distrib(),            list(mu = 1, sigma = 0.7)),
    list("student_t1",  student_t1_distrib(),          list(mu = 1, sigma = 2, nu = 9)),
    list("gamma2",      gamma2_distrib(),              list(mu = 3, sigma2 = 2)),
    list("exponential", exponential_distrib(),         list(mu = 2.5)),
    list("chisq",       chisq_distrib(),               list(mu = 7)),
    list("lognormal1",  lognormal1_distrib(),          list(mu = 0.3, sigma2 = 0.4)),
    list("invgauss1",   invgauss1_distrib(),           list(mu = 2, phi = 0.3)),
    list("beta1",       beta1_distrib(),               list(mu = 0.3, phi = 9)),
    list("gpd",         gpd_distrib(),                 list(sigma = 2, xi = 0.15)),
    list("gengamma1",   gengamma1_distrib(),           list(a = 2, d = 3, p = 1.5)),
    list("poisson",     poisson_distrib(),             list(mu = 4)),
    list("bernoulli",   bernoulli_distrib(),           list(mu = 0.3)),
    list("binomial",    binomial_distrib(size = 10),   list(mu = 0.3)),
    list("geometric",   geometric_distrib(),           list(mu = 2)),
    list("negbin1",     negbin1_distrib(),             list(mu = 5, theta = 1.5)),
    list("negbin2",     negbin2_distrib(),             list(mu = 5, theta = 3)),
    list("betabinom1",  betabinom1_distrib(size = 12), list(mu = 0.3, sigma = 0.5))
  )
}


test_that("every closed form agrees with the numerical route", {
  for (cs in moment_cases()) {
    d <- cs[[2]]
    th <- cs[[3]]
    a <- c(mean = mean(d, th), var = variance(d, th),
           skew = skewness(d, th), kurt = kurtosis(d, th))
    e <- num_moments(d, th)
    expect_equal(unname(a), unname(e), tolerance = 1e-5, label = cs[[1]])
  }
})


test_that("a discrete family agrees to machine precision", {
  # For a discrete family the numerical route is an exact sum rather than a
  # quadrature, so the comparison has no truncation error to hide behind.
  for (cs in moment_cases()) {
    d <- cs[[2]]
    if (!S7::S7_inherits(d, discrete_distrib)) next
    th <- cs[[3]]
    a <- c(mean(d, th), variance(d, th), skewness(d, th), kurtosis(d, th))
    expect_lt(max(abs(a - unname(num_moments(d, th)))), 1e-12,
              label = cs[[1]])
  }
})


test_that("the values that are exact come back exact", {
  # These are the ones a quadrature could only approach: the motivating case
  # was the gaussian excess kurtosis, which the fallback returns as -4.8e-07.
  expect_identical(kurtosis(gaussian1_distrib(), list(mu = 0, sigma = 2)), 0)
  expect_identical(skewness(gaussian1_distrib(), list(mu = 0, sigma = 2)), 0)
  expect_identical(skewness(logistic_distrib(), list(mu = 3, sigma = 2)), 0)
  expect_identical(kurtosis(logistic_distrib(), list(mu = 3, sigma = 2)), 6 / 5)
  expect_identical(skewness(exponential_distrib(), list(mu = 3)), 2)
  expect_identical(kurtosis(exponential_distrib(), list(mu = 3)), 6)
  expect_identical(variance(chisq_distrib(), list(mu = 7)), 14)
  expect_identical(variance(gaussian1_distrib(), list(mu = 0, sigma = 3)), 9)
})


test_that("a moment that does not exist is NaN, not an error", {
  # The Cauchy has no moment at all. The numerical route does not merely lose
  # accuracy on it, it raises "maximum number of subdivisions reached", which
  # names the quadrature rather than the mathematics.
  d <- cauchy_distrib()
  th <- list(mu = 0, sigma = 1)
  expect_true(is.nan(mean(d, th)))
  expect_true(is.nan(variance(d, th)))
  expect_true(is.nan(skewness(d, th)))
  expect_true(is.nan(kurtosis(d, th)))
  expect_error(moment(d, th, p = 1, central = FALSE))
})


test_that("the Student t moments appear one threshold at a time", {
  d <- student_t1_distrib()
  th <- function(nu) list(mu = 0, sigma = 1, nu = nu)

  expect_true(is.nan(mean(d, th(0.5))))
  expect_equal(mean(d, th(1.5)), 0)

  expect_true(is.nan(variance(d, th(0.5))))
  expect_identical(variance(d, th(1.5)), Inf)
  expect_equal(variance(d, th(6)), 1.5)

  expect_true(is.nan(skewness(d, th(2.5))))
  expect_equal(skewness(d, th(3.5)), 0)

  expect_true(is.nan(kurtosis(d, th(1.5))))
  expect_identical(kurtosis(d, th(3.5)), Inf)
  expect_equal(kurtosis(d, th(6)), 3)
})


test_that("the generalized Pareto moments follow xi < 1/k", {
  d <- gpd_distrib()
  th <- function(xi) list(sigma = 1, xi = xi)

  expect_equal(mean(d, th(0.1)), 1 / 0.9)
  expect_identical(mean(d, th(1.2)), Inf)
  expect_true(is.finite(variance(d, th(0.4))))
  expect_identical(variance(d, th(0.6)), Inf)
  expect_true(is.finite(skewness(d, th(0.3))))
  expect_identical(skewness(d, th(0.4)), Inf)
  expect_true(is.finite(kurtosis(d, th(0.2))))
  expect_identical(kurtosis(d, th(0.3)), Inf)
})


test_that("the beta-binomial reduces to the families it contains", {
  # Two independent implementations of one object, so no tolerance has to be
  # chosen: at size 1 with a vanishing dispersion it is a Bernoulli, and at a
  # vanishing dispersion with any size it is a binomial.
  eps <- 1e-9
  b1 <- betabinom1_distrib(size = 1)
  be <- bernoulli_distrib()
  for (p in c(0.2, 0.5, 0.8)) {
    expect_equal(mean(b1, list(mu = p, sigma = eps)), mean(be, list(mu = p)),
                 tolerance = 1e-6, label = as.character(p))
    expect_equal(kurtosis(b1, list(mu = p, sigma = eps)),
                 kurtosis(be, list(mu = p)), tolerance = 1e-6,
                 label = as.character(p))
  }
  bb <- betabinom1_distrib(size = 12)
  bn <- binomial_distrib(size = 12)
  expect_equal(variance(bb, list(mu = 0.3, sigma = eps)),
               variance(bn, list(mu = 0.3)), tolerance = 1e-6)
  expect_equal(skewness(bb, list(mu = 0.3, sigma = eps)),
               skewness(bn, list(mu = 0.3)), tolerance = 1e-6)
})


test_that("NB1 and NB2 tend to the Poisson from their own directions", {
  po <- poisson_distrib()
  mu <- 4
  expect_equal(variance(negbin1_distrib(), list(mu = mu, theta = 1e-10)),
               variance(po, list(mu = mu)), tolerance = 1e-8)
  expect_equal(skewness(negbin1_distrib(), list(mu = mu, theta = 1e-10)),
               skewness(po, list(mu = mu)), tolerance = 1e-8)
  expect_equal(kurtosis(negbin2_distrib(), list(mu = mu, theta = 1e10)),
               kurtosis(po, list(mu = mu)), tolerance = 1e-8)
})


test_that("the moments recycle over vector parameters", {
  d <- gaussian1_distrib()
  th <- list(mu = 0, sigma = c(1, 2, 3))
  expect_identical(variance(d, th), c(1, 4, 9))
  expect_length(mean(d, th), 3L)
  expect_length(skewness(d, th), 3L)
  expect_length(kurtosis(d, th), 3L)

  g <- gamma2_distrib()
  thg <- list(mu = c(1, 2), sigma2 = c(0.5, 3))
  expect_identical(variance(g, thg), c(0.5, 3))
  expect_length(skewness(g, thg), 2L)

  # a constant moment still takes the length the parameters imply
  expect_length(skewness(cauchy_distrib(), list(mu = c(0, 1, 2), sigma = 1)), 3L)
})


test_that("a wrong closed form would be caught", {
  # The paired injection: the agreement above is only evidence if a corrupted
  # formula fails it.
  good <- S7::method(variance, Gamma2Distrib)
  on.exit(S7::method(variance, Gamma2Distrib) <- good, add = TRUE)
  suppressMessages(
    S7::method(variance, Gamma2Distrib) <- function(x, theta, ...) {
      good(x, theta, ...) * 1.05
    }
  )
  d <- gamma2_distrib()
  th <- list(mu = 3, sigma2 = 2)
  expect_false(isTRUE(all.equal(variance(d, th),
                                unname(num_moments(d, th)[["var"]]),
                                tolerance = 1e-5)))
})


test_that("a multivariate family refuses a scalar skewness by name", {
  # The refusal used to surface as "Can't find method for expectation(...)",
  # which names the quadrature the base class fell through to rather than the
  # reason. Mardia's skewness, Malkovich-Afifi's and the vector of marginal
  # skewnesses are three different quantities, so the bare name is refused.
  for (d in list(mvgaussian_distrib(2), mvstudent_t_distrib(2),
                 dirichlet_distrib(3), multinomial_distrib(3, size = 4))) {
    th <- generate_random_theta(d)
    expect_error(skewness(d, th), "no single skewness")
    expect_error(kurtosis(d, th), "no single kurtosis")
  }

  # mv_marginal() is not a way round it for an elliptical family: the marginal
  # of one coordinate is a one-dimensional MULTIVARIATE gaussian, which refuses
  # in turn. It is a way round it for the two families whose marginals change
  # family, and there the univariate answer is the one to have.
  d <- mvgaussian_distrib(2)
  th <- list(mu1 = 1, mu2 = -1, sigma_log_L1 = 0, sigma_log_L2 = 0,
             sigma_L2.1 = 0.3)
  m <- mv_marginal(d, th, which = 1L)
  expect_true(S7::S7_inherits(m$distrib, multivariate_distrib))
  expect_error(skewness(m$distrib, m$theta), "no single skewness")

  dd <- dirichlet_distrib(3)
  thd <- list(mean_alr1 = 0.2, mean_alr2 = -0.1, phi = 8)
  md <- mv_marginal(dd, thd, which = 1L)
  expect_true(S7::S7_inherits(md$distrib, Beta1Distrib))
  expect_equal(skewness(md$distrib, md$theta),
               skewness(beta1_distrib(), md$theta))
})


test_that("the von Mises mean is not its mean direction", {
  # mu is the mean DIRECTION. The ordinary expectation of Y as a number on
  # [-pi, pi) differs from it whenever mu is not zero, because the interval is
  # cut at +/- pi rather than at mu +/- pi, and that is what mean() returns.
  d <- vonmises1_distrib()
  expect_equal(mean(d, list(mu = 0, kappa = 2)), 0, tolerance = 1e-6)
  em <- mean(d, list(mu = 1.2, kappa = 2))
  expect_lt(em, 1.2)
  expect_equal(em, 1.079, tolerance = 1e-3)
})
