# Strategies for approximating expected derivatives when no closed form exists:
# "opg", "bartlett", "integrate" and "mc". The first two are the two sides of
# the second Bartlett identity taken literally -- "bartlett" evaluates the
# expectation, "opg" reads the integrand at the observation -- so they agree in
# expectation and not per observation. See test-expected-opg.R.

BareG <- S7::new_class("BareG", parent = continuous_distrib, package = NULL)
S7::method(distrib_pdf, BareG) <- function(distrib, y, theta, log = FALSE) {
  stats::dnorm(y, theta[[1]], theta[[2]], log = log)
}
S7::method(distrib_rng, BareG) <- function(distrib, n, theta) {
  stats::rnorm(n, theta[[1]], theta[[2]])
}

bare_gauss_rng <- function() {
  BareG(
    distrib_name = "bare gaussian", dimension = "univariate", bounds = c(-Inf, Inf),
    params = c("mu", "sigma"), params_interpretation = c(mu = "m", sigma = "s"),
    n_params = 2, params_bounds = list(mu = c(-Inf, Inf), sigma = c(0, Inf)),
    link_params = list(mu = linkfunctions7::identity_link(),
                       sigma = linkfunctions7::log_link())
  )
}

test_that("all strategies recover the analytical expected Hessian", {
  bg <- bare_gauss_rng()
  th <- list(mu = 1.5, sigma = 2)
  truth <- distrib_expected_hessian(gaussian1_distrib(), 0, th)

  # "opg" is deliberately not in this loop: it estimates this quantity rather
  # than evaluating it, so per observation it is a rank-one matrix read at the
  # data point. Its own claim is the test below.
  for (a in c("bartlett", "integrate")) {
    e <- distrib_expected_hessian(bg, 0, th, approx = a)
    for (k in names(truth)) {
      expect_equal(e[[k]][1], truth[[k]][1], tolerance = 1e-5,
                   label = paste("expected hessian", a, k))
    }
  }

  set.seed(11)
  e <- distrib_expected_hessian(bg, 0, th, approx = "mc", nsim = 1e5)
  for (k in names(truth)) {
    expect_equal(e[[k]][1], truth[[k]][1], tolerance = 0.02,
                 label = paste("expected hessian mc", k))
  }
})

test_that("'opg' is NOT an alias for 'bartlett', and estimates what it computes", {
  # It used to be one: `expected_derivative()` rewrote "opg" as "bartlett", on
  # the reading that at order 2 the outer product of gradients IS the second
  # Bartlett identity. The identity is the same and the computation is not --
  # one evaluates the expectation, the other reads the integrand at the
  # observation -- so the two differ per observation and agree in the mean.
  bg <- bare_gauss_rng()
  th <- list(mu = 1, sigma = 2)
  opg <- distrib_expected_hessian(bg, 0, th, approx = "opg")
  bar <- distrib_expected_hessian(bg, 0, th, approx = "bartlett")
  expect_false(isTRUE(all.equal(opg, bar)))

  # what opg returns, written out: -s_i s_j at the observation
  g <- distrib_gradient(bg, 0, th)
  expect_equal(opg$mu_mu, -g$mu * g$mu)
  expect_equal(opg$mu_sigma, -g$mu * g$sigma)

  # and its mean over a sample is the quantity bartlett evaluates
  set.seed(9)
  y <- distrib_rng(bg, 30000L, th)
  opg_n <- distrib_expected_hessian(bg, y, th, approx = "opg")
  for (k in names(bar)) {
    expect_equal(mean(opg_n[[k]]), bar[[k]][1], tolerance = 0.05, label = k)
  }
})

test_that("approx is ignored when the distribution has a closed form", {
  g <- gaussian1_distrib()
  th <- list(mu = 1.5, sigma = 2)
  base <- distrib_expected_hessian(g, 0, th)
  for (a in c("opg", "bartlett", "integrate", "mc")) {
    expect_equal(distrib_expected_hessian(g, 0, th, approx = a), base)
  }
})

BareLap <- S7::new_class("BareLap", parent = continuous_distrib, package = NULL)
S7::method(distrib_pdf, BareLap) <- function(distrib, y, theta, log = FALSE) {
  ld <- -log(2 * theta[[2]]) - abs(y - theta[[1]]) / theta[[2]]
  if (log) ld else exp(ld)
}
S7::method(distrib_gradient, BareLap) <- function(distrib, y, theta,
                                                  scale = c("parameter", "link"), ...) {
  r <- y - theta[[1]]
  list(mu = sign(r) / theta[[2]], b = (abs(r) / theta[[2]] - 1) / theta[[2]])
}

test_that("only the score-based forms give the information for a non-regular model", {
  bl <- BareLap(
    distrib_name = "bare laplace", dimension = "univariate", bounds = c(-Inf, Inf),
    params = c("mu", "b"), params_interpretation = c(mu = "loc", b = "scale"),
    n_params = 2, params_bounds = list(mu = c(-Inf, Inf), b = c(0, Inf)),
    link_params = list(mu = linkfunctions7::identity_link(),
                       b = linkfunctions7::log_link()),
    params_smooth = c(mu = FALSE, b = TRUE)
  )
  th <- list(mu = 1, b = 2)

  # Bartlett = -E[score^2] = -1/b^2, the true Fisher information
  bart <- distrib_expected_hessian(bl, 0, th, approx = "bartlett")
  expect_equal(bart$mu_mu[1], -1 / th$b^2, tolerance = 1e-6)
  expect_equal(bart$b_b[1], -1 / th$b^2, tolerance = 1e-6)

  # quadrature of the observed Hessian misses the kink and returns 0 for mu
  integ <- distrib_expected_hessian(bl, 0, th, approx = "integrate")
  expect_equal(integ$mu_mu[1], 0, tolerance = 1e-6)

  # opg, the default, is EXACT here and needs no expectation at all: the
  # squared score of a Laplace is 1/b^2 almost surely rather than in
  # expectation, so the integrand is constant. y = 0 is not the kink, mu is 1.
  opg <- distrib_expected_hessian(bl, 0, th, approx = "opg")
  expect_equal(opg$mu_mu[1], -1 / th$b^2)

  # AT the kink it is zero instead, sign(0) being 0, and that is the one place
  # the per-observation reading parts company with the expectation. The set has
  # probability zero, so nothing a fit aggregates is affected -- but the claim
  # above must not be written without this one beside it.
  at_kink <- distrib_expected_hessian(bl, th$mu, th, approx = "opg")
  expect_equal(at_kink$mu_mu[1], 0)
})

test_that("higher-order strategies agree for a distribution with analytic derivatives", {
  st <- student_t1_distrib()
  th <- list(mu = 0, sigma = 1, nu = 8)
  r_int <- distrib_deriv3(st, 0, th, expected = TRUE, approx = "integrate")
  r_bar <- distrib_deriv3(st, 0, th, expected = TRUE, approx = "bartlett")
  for (k in names(r_int)) {
    expect_equal(r_bar[[k]][1], r_int[[k]][1], tolerance = 1e-6, label = paste("d3", k))
  }
})

test_that("the order-3 Bartlett identity is exact for the Poisson", {
  p <- poisson_distrib()
  th <- list(mu = 4)
  analytic <- distrib_deriv3(p, 0, th, expected = TRUE)$mu_mu_mu
  bart <- distrib_deriv3(p, 0, th, expected = TRUE, approx = "bartlett")$mu_mu_mu
  expect_equal(bart, analytic, tolerance = 1e-8)
})

test_that("integrate fails informatively on a purely numerical high-order derivative", {
  bg <- bare_gauss_rng()
  th <- list(mu = 1.5, sigma = 2)
  expect_error(
    distrib_deriv3(bg, 0, th, expected = TRUE, approx = "integrate"),
    "bartlett"
  )
  # and the suggested alternative works
  ref <- distrib_deriv3(gaussian1_distrib(), 0, th, expected = TRUE)
  alt <- distrib_deriv3(bg, 0, th, expected = TRUE, approx = "bartlett")
  for (k in names(ref)) {
    expect_equal(alt[[k]][1], ref[[k]][1], tolerance = 1e-5, label = paste("d3", k))
  }
})

test_that("the partition enumeration comes from numericals7", {
  # The local copy is gone -- it was one of the two independent copies whose
  # existence motivated numericals7 -- so what remains to check here is that
  # the Bartlett machinery reaches the shared one.
  expect_identical(lengths(lapply(1:4, numericals7::set_partitions)),
                   c(1L, 2L, 5L, 15L))
})
