# The fallbacks of orders 2 to 4 for a family without its own Hessian: one
# tensor-product stencil on the score or on the log-density. The reference is the
# shipped gamma family's analytical kernels, which these user classes share a
# density with.

GammaPdfOnly <- S7::new_class("GammaPdfOnly", parent = continuous_distrib)
S7::method(distrib_pdf, GammaPdfOnly) <- function(distrib, y, theta, log = FALSE, ...) {
  distrib_pdf(gamma1_distrib(), y, theta, log = log)
}
GammaWithScore <- S7::new_class("GammaWithScore", parent = GammaPdfOnly)
S7::method(distrib_gradient, GammaWithScore) <- function(distrib, y, theta,
                                                         scale = c("parameter", "link"), ...) {
  distrib_gradient(gamma1_distrib(), y, theta)
}
gamma_like <- function(cls) {
  g <- gamma1_distrib()
  cls(distrib_name = "gamma copy", dimension = "univariate", bounds = g@bounds,
      params = g@params, params_interpretation = g@params_interpretation,
      n_params = g@n_params, params_bounds = g@params_bounds,
      link_params = g@link_params)
}

y <- c(0.3, 1.1, 2.4, 5.0)
th <- list(mu = c(1.5, 1.5, 2.5, 2.5), phi = c(0.4, 0.4, 0.8, 0.8))
rel <- function(a, b) max(abs(unlist(a) - unlist(b))) / max(abs(unlist(b)))

test_that("analytic_order() reads what the family implements", {
  expect_identical(analytic_order(gamma_like(GammaPdfOnly), 2L), 0L)
  expect_identical(analytic_order(gamma_like(GammaWithScore), 2L), 1L)
  expect_identical(analytic_order(gamma1_distrib(), 2L), 2L)
})

test_that("a density alone gives orders 3 and 4 from one stencil each", {
  d <- gamma_like(GammaPdfOnly)
  ref3 <- distrib_deriv3(gamma1_distrib(), y, th)
  ref4 <- distrib_deriv4(gamma1_distrib(), y, th)
  got3 <- distrib_deriv3(d, y, th)
  got4 <- distrib_deriv4(d, y, th)
  expect_identical(names(got3), names(ref3))
  expect_identical(names(got4), names(ref4))
  expect_lt(rel(got3, ref3), 1e-4)
  expect_lt(rel(got4, ref4), 1e-3)
})

test_that("a score gives the Hessian and orders 3 and 4 from the score", {
  d <- gamma_like(GammaWithScore)
  got2 <- distrib_hessian(d, y, th)
  expect_identical(names(got2), hess_names(d@params))
  ref2 <- distrib_hessian(gamma1_distrib(), y, th)[names(got2)]
  expect_lt(rel(got2, ref2), 1e-8)
  expect_lt(rel(distrib_deriv3(d, y, th), distrib_deriv3(gamma1_distrib(), y, th)), 1e-6)
  expect_lt(rel(distrib_deriv4(d, y, th), distrib_deriv4(gamma1_distrib(), y, th)), 1e-4)
})

test_that("the stencil matches exact derivatives where nesting missed them", {
  # A generalized Poisson with a probability function alone. Differencing its
  # numerical Hessian missed the fourth order by factors of 3 to 71; the exact
  # reference is symbolic differentiation of the log probability.
  GP <- S7::new_class("GP", parent = discrete_distrib)
  S7::method(distrib_pdf, GP) <- function(distrib, y, theta, log = FALSE, ...) {
    mu <- theta[[1]]; phi <- theta[[2]]
    s <- sqrt(phi); a <- mu / s + (1 - 1 / s) * y
    ld <- log(mu / s) + (y - 1) * log(a) - a - lgamma(y + 1)
    if (log) ld else exp(ld)
  }
  d <- GP(distrib_name = "gp", dimension = "univariate", bounds = c(0, Inf),
          params = c("mu", "phi"), params_interpretation = c(mu = "mean", phi = "dispersion"),
          n_params = 2, params_bounds = list(mu = c(0, Inf), phi = c(0, Inf)),
          link_params = list(mu = log_link(), phi = log_link()))
  yy <- c(10, 20, 26, 40, 70)
  tt <- list(mu = rep(28, 5), phi = rep(6, 5))
  e <- quote(log(mu / sqrt(phi)) + (y - 1) * log(mu / sqrt(phi) + (1 - 1 / sqrt(phi)) * y) -
               (mu / sqrt(phi) + (1 - 1 / sqrt(phi)) * y))
  exact <- function(nm) {
    ex <- e
    for (v in strsplit(nm, "_")[[1]]) ex <- D(ex, v)
    eval(ex, list(mu = 28, phi = 6, y = yy))
  }
  g3 <- distrib_deriv3(d, yy, tt)
  g4 <- distrib_deriv4(d, yy, tt)
  for (nm in names(g3)) expect_lt(rel(g3[[nm]], exact(nm)), 1e-4, label = nm)
  for (nm in names(g4)) expect_lt(rel(g4[[nm]], exact(nm)), 1e-3, label = nm)
})


test_that("skip and an explicit step reach the stencil route", {
  d <- gamma_like(GammaPdfOnly)
  s <- numerical_deriv3(d, y, th, skip = "mu_mu_mu")
  expect_null(s$mu_mu_mu)
  expect_identical(names(s), deriv_names(d@params, 3))
  a <- numerical_deriv3(d, y, th, h_rel = 1e-3)
  b <- numerical_deriv3(d, y, th)
  expect_false(identical(a, b))
  expect_lt(rel(a, distrib_deriv3(gamma1_distrib(), y, th)), 1e-3)
})

test_that("a shipped family keeps the route through its own Hessian", {
  d <- gamma1_distrib()
  h <- .Machine$double.eps^(1 / 3)
  expect_identical(numerical_deriv3(d, y, th), numerical_deriv3(d, y, th, h_rel = h))
})
