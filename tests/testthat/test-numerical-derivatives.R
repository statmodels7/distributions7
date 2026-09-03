# Fallback numerical derivatives: a distribution that implements ONLY
# distrib_pdf must still provide gradient, observed Hessian and expected
# Hessian through the default methods on the base `distrib` class.

# A "bare" Gaussian: same math as gaussian1_distrib(), but with no analytical
# derivative methods registered, so every derivative call hits the fallbacks.
BareGauss <- S7::new_class("BareGauss", parent = continuous_distrib, package = NULL)
S7::method(distrib_pdf, BareGauss) <- function(distrib, y, theta, log = FALSE) {
  stats::dnorm(y, mean = theta[[1]], sd = theta[[2]], log = log)
}
S7::method(distrib_quantile, BareGauss) <- function(distrib, p, theta, ...) {
  stats::qnorm(p, mean = theta[[1]], sd = theta[[2]])
}
S7::method(distrib_rng, BareGauss) <- function(distrib, n, theta) {
  stats::rnorm(n, mean = theta[[1]], sd = theta[[2]])
}

bare_gauss <- function() {
  BareGauss(
    distrib_name = "bare gaussian",
    dimension = "univariate",
    bounds = c(-Inf, Inf),
    params = c("mu", "sigma"),
    params_interpretation = c(mu = "mean", sigma = "standard deviation"),
    n_params = 2,
    params_bounds = list(mu = c(-Inf, Inf), sigma = c(0, Inf)),
    link_params = list(
      mu = linkfunctions7::identity_link(),
      sigma = linkfunctions7::log_link()
    )
  )
}

test_that("fallback gradient matches the analytical one", {
  set.seed(61)
  bg <- bare_gauss()
  ref <- gaussian1_distrib()
  th <- list(mu = 1.5, sigma = 2.0)
  y <- distrib_rng(ref, 30, th)

  g_num <- distrib_gradient(bg, y, th)
  g_ana <- distrib_gradient(ref, y, th)

  expect_named(g_num, c("mu", "sigma"))
  expect_equal(g_num$mu, g_ana$mu, tolerance = 1e-7)
  expect_equal(g_num$sigma, g_ana$sigma, tolerance = 1e-7)
})

test_that("fallback hessian matches the analytical one, in hess_names order", {
  set.seed(62)
  bg <- bare_gauss()
  ref <- gaussian1_distrib()
  th <- list(mu = 1.5, sigma = 2.0)
  y <- distrib_rng(ref, 30, th)

  h_num <- distrib_hessian(bg, y, th)
  h_ana <- distrib_hessian(ref, y, th)

  expect_named(h_num, hess_names(c("mu", "sigma")))
  for (nm in names(h_num)) {
    expect_equal(h_num[[nm]], h_ana[[nm]], tolerance = 1e-6, label = paste("hess", nm))
  }
})

test_that("fallback expected hessian matches the analytical one", {
  # approx = "bartlett" is asked for explicitly: it evaluates the expectation,
  # which is the quantity being checked against the closed form. "opg" (the
  # default since 0.44.0) reads the score at each observation instead of
  # averaging it, so it would vary across y = c(0, 1, 2) where the closed form
  # is constant -- see test-expected-opg.R for what it agrees with and where.
  bg <- bare_gauss()
  ref <- gaussian1_distrib()
  th <- list(mu = 1.5, sigma = 2.0)

  eh_num <- distrib_expected_hessian(bg, c(0, 1, 2), th, approx = "bartlett")
  eh_ana <- distrib_expected_hessian(ref, c(0, 1, 2), th)

  expect_named(eh_num, hess_names(c("mu", "sigma")))
  expect_length(eh_num$mu_mu, 3)
  for (nm in names(eh_num)) {
    expect_equal(eh_num[[nm]], eh_ana[[nm]], tolerance = 1e-5, label = paste("expected hess", nm))
  }
})

test_that("fallback derivatives are vectorized over theta and validated like the analytic ones", {
  bg <- bare_gauss()
  ref <- gaussian1_distrib()
  y <- c(-1, 0, 3)
  th_vec <- list(mu = c(0, 1, 2), sigma = c(1, 2, 3))

  expect_equal(
    distrib_gradient(bg, y, th_vec)$mu,
    distrib_gradient(ref, y, th_vec)$mu,
    tolerance = 1e-7
  )

  # Same argument validation path as analytic methods (generic-level checks)
  expect_error(distrib_gradient(bg, rnorm(10), list(mu = 1:3, sigma = 1)), "dimension mismatch")
  expect_error(distrib_gradient(bg, 0, list(mu = 0)), "Missing parameter")
})

test_that("finite-difference steps respect parameter domain boundaries", {
  bg <- bare_gauss()
  # sigma extremely close to its lower bound 0: steps must shrink, not cross it
  g <- distrib_gradient(bg, 0.1, list(mu = 0, sigma = 1e-4))
  expect_true(all(is.finite(unlist(g))))

  # exactly on the boundary -> rejected up front by the domain check
  expect_error(
    numerical_gradient(bg, 0.1, list(mu = 0, sigma = 0)),
    "outside its domain"
  )
  # the finite-difference step helper keeps its own guard for internal use
  expect_error(fd_steps(0, c(0, Inf), 1e-5), "boundary")
})

test_that("numerical_gradient/numerical_hessian are exported and usable directly", {
  d <- gaussian1_distrib()
  th <- list(mu = 1, sigma = 2)
  y <- c(-1, 0.5, 4)

  expect_equal(numerical_gradient(d, y, th)$mu, (y - 1) / 4, tolerance = 1e-7)
  expect_equal(
    numerical_hessian(d, y, th)$mu_mu,
    rep(-1 / 4, 3),
    tolerance = 1e-6
  )
})
