# Third and fourth derivatives with respect to the response.
#
# For a family whose response enters only as y - mu there is nothing to
# derive: d^k l / dy^k = (-1)^k d^k l / dmu^k. The identity is checked at
# orders one and two, where both sides already exist and it must hold
# exactly, and the new orders are then checked against one stencil on the
# log-density.

LOCATION_FAMILIES <- function() {
  list(
    list(gaussian1_distrib(), list(mu = 0.3, sigma = 1.4)),
    list(gaussian2_distrib(), list(mu = 0.3, sigma2 = 2)),
    list(gaussian3_distrib(), list(mu = 0.3, tau = 0.5)),
    list(cauchy_distrib(), list(mu = 0.2, sigma = 1.1)),
    list(logistic_distrib(), list(mu = 0.2, sigma = 1.1)),
    list(laplace_distrib(), list(mu = 0.2, sigma = 1.1)),
    list(laplace2_distrib(), list(mu = 0.2, lambda = 0.8)),
    list(enet_distrib(), list(mu = 0.2, lambda = 1.5, alpha = 0.6)),
    list(pseudohuber_distrib(), list(mu = 0.2, sigma = 1.1, nu = 3)),
    list(student_t1_distrib(), list(mu = 0.2, sigma = 1.1, nu = 6)),
    list(skewnormal1_distrib(), list(mu = 0.5, sigma = 1.5, alpha = 2)),
    list(skewnormal2_distrib(), list(mu = 0.5, sigma = 1.5, gamma1 = 0.4)),
    list(skewt_distrib(), list(mu = 0.5, sigma = 1.5, alpha = 2, nu = 6)),
    list(gumbel_distrib(), list(mu = 0.5, sigma = 1.5))
  )
}

test_that("the identity holds exactly where both sides already exist", {
  # dl/dy = -dl/dmu and d2l/dy2 = d2l/dmu2. These two orders are written
  # independently of each other in every family, so exact agreement is what
  # says the identity may be used at the orders above, where only one side
  # is written.
  y <- c(-1.7, 0.9, 2.4)
  for (cs in LOCATION_FAMILIES()) {
    d <- cs[[1]]; th <- cs[[2]]
    expect_equal(distrib_grad_y(d, y, th),
                 -distrib_gradient(d, y, th)$mu,
                 tolerance = 1e-12, label = d@distrib_name)
    expect_equal(distrib_hess_y(d, y, th),
                 rep_len(distrib_hessian(d, y, th)$mu_mu, length(y)),
                 tolerance = 1e-12, label = d@distrib_name)
  }
})

test_that("orders three and four agree with one stencil on the log-density", {
  y <- c(-1.7, 0.9, 2.4)
  for (cs in LOCATION_FAMILIES()) {
    d <- cs[[1]]; th <- cs[[2]]
    for (ord in 3:4) {
      got <- if (ord == 3L) distrib_deriv3_y(d, y, th) else
                            distrib_deriv4_y(d, y, th)
      ref <- numerical_deriv_y(d, y, th, ord)
      # the stencil's own accuracy, which for a fourth difference is bounded
      # by the eps/h^4 rounding it carries
      tol <- if (ord == 3L) 1e-4 else 5e-3
      expect_lt(max(abs(got - ref)) / max(abs(got), 1), tol,
                label = paste(d@distrib_name, ord))
    }
  }
})

dy_families <- function() {
  list(
    list("gamma1", gamma1_distrib(), list(mu = 2, phi = 0.6), c(0.4, 1.5, 4)),
    list("gamma2", gamma2_distrib(), list(mu = 2, sigma2 = 1.1), c(0.4, 1.5, 4)),
    list("chisq", chisq_distrib(), list(mu = 5), c(1, 3, 8)),
    list("exponential", exponential_distrib(), list(mu = 2), c(0.3, 1, 4)),
    list("beta1", beta1_distrib(), list(mu = 0.4, phi = 6), c(0.2, 0.5, 0.8)),
    list("beta2", beta2_distrib(), list(alpha = 2, beta = 3), c(0.2, 0.5, 0.8)),
    list("weibull1", weibull1_distrib(), list(mu = 2, sigma = 1.5), c(0.5, 2, 4)),
    list("weibull3", weibull3_distrib(), list(mean = 2, sigma = 1.5), c(0.5, 2, 4)),
    list("gengamma1", gengamma1_distrib(), list(a = 2, d = 3, p = 1.4), c(0.6, 2, 4)),
    list("gengamma2", gengamma2_distrib(), list(mean = 2, d = 3, p = 1.4), c(0.6, 2, 4)),
    list("invgauss1", invgauss1_distrib(), list(mu = 2, phi = 0.5), c(0.5, 2, 5)),
    list("invgauss2", invgauss2_distrib(), list(mu = 2, lambda = 3), c(0.5, 2, 5)),
    list("lognormal1", lognormal1_distrib(), list(mu = 0.3, sigma2 = 0.7), c(0.4, 1.2, 4)),
    list("lognormal2", lognormal2_distrib(), list(mean = 1.5, var = 1.2), c(0.4, 1.2, 4)),
    list("gpd", gpd_distrib(), list(sigma = 1.4, xi = 0.3), c(0.3, 1.5, 5)),
    list("vonmises1", vonmises1_distrib(), list(mu = 0.4, kappa = 2), c(-2, 0, 1.5)),
    list("vonmises2", vonmises2_distrib(), list(mu = 0.4, rho = 0.6), c(-2, 0, 1.5)),
    list("student_t2", student_t2_distrib(), list(mu = 0.2, sigma = 1.3, nu = 6),
         c(-2, 0.5, 3))
  )
}

test_that("a family whose response is not a pure location is closed too", {
  # every one of these already carried a closed second response derivative, so
  # the third is the same elementary function once more; the reference is one
  # differentiation of that analytic second derivative, never of the third,
  # which would be a difference of a difference
  skip_if_not_installed("numDeriv")
  for (cs in dy_families()) {
    nm <- cs[[1]]; d <- cs[[2]]; th <- cs[[3]]; y <- cs[[4]]
    r3 <- vapply(seq_along(y), function(i)
      numDeriv::grad(function(v) distrib_hess_y(d, v, th), y[i]), numeric(1))
    expect_lt(max(abs(distrib_deriv3_y(d, y, th) - r3)) / max(1, max(abs(r3))),
              1e-8, label = paste(nm, "order 3"))
    # order 4 against one second-order stencil on the same analytic quantity
    h <- 1e-3 * pmax(1, abs(y))
    r4 <- (distrib_hess_y(d, y + h, th) - 2 * distrib_hess_y(d, y, th) +
             distrib_hess_y(d, y - h, th)) / h^2
    expect_lt(max(abs(distrib_deriv4_y(d, y, th) - r4)) / max(1, max(abs(r4))),
              1e-3, label = paste(nm, "order 4"))
  }
})

test_that("no continuous family is left on the stencil", {
  # the audit this batch was built from, run as a test so that a family added
  # later without a closed rule is caught
  exported <- grep("_distrib$", getNamespaceExports("distributions7"), value = TRUE)
  skip_these <- c("check_distrib", "fit_distrib", "continuous_distrib",
                  "discrete_distrib", "multivariate_distrib")
  open <- character()
  for (ctor in setdiff(exported, skip_these)) {
    d <- tryCatch(get(ctor, asNamespace("distributions7"))(),
                  error = function(e) NULL)
    if (is.null(d) || !S7::S7_inherits(d, distributions7:::continuous_distrib)) next
    m <- S7::method(distrib_deriv3_y, S7::S7_class(d))
    if (identical(attr(attr(m, "signature")[[1L]], "name"), "continuous_distrib")) {
      open <- c(open, ctor)
    }
  }
  expect_identical(open, character(),
                   info = paste("still on the stencil:", paste(open, collapse = ", ")))
})

test_that("the generalized Pareto stays finite as the shape goes to zero", {
  # the coefficient is written as xi^k + xi^(k-1) rather than (1 + 1/xi) xi^k,
  # which is the same number and does not divide by a vanishing shape; at xi = 0
  # the family is exponential and every order above the first is exactly zero
  d <- gpd_distrib()
  y <- c(0.3, 1.5, 5)
  for (xi in c(1e-12, 1e-9, 1e-6)) {
    th <- list(sigma = 1.4, xi = xi)
    expect_true(all(is.finite(distrib_deriv3_y(d, y, th))))
    expect_lt(max(abs(distrib_deriv3_y(d, y, th))), 1e-5)
    expect_lt(max(abs(distrib_deriv4_y(d, y, th))), 1e-5)
  }
  expect_equal(distrib_deriv3_y(exponential_distrib(), y, list(mu = 2)),
               rep(0, 3), tolerance = 1e-15)
})

test_that("the step stays inside the support", {
  # the stencil reaches two steps either side, so a point near a boundary
  # must not be differentiated outside it
  d <- gamma1_distrib()
  th <- list(mu = 2, phi = 0.6)
  expect_true(all(is.finite(distrib_deriv4_y(d, c(1e-6, 1e-3, 0.01), th))))
  b <- beta1_distrib()
  expect_true(all(is.finite(distrib_deriv4_y(b, c(1e-6, 1 - 1e-6),
                                             list(mu = 0.4, phi = 6)))))
})

test_that("a discrete family has no response derivative at these orders", {
  # as at the orders below: a difference across the lattice is not a
  # derivative, and returning one would be a fiction
  expect_error(distrib_deriv3_y(poisson_distrib(), 2, list(mu = 3)))
  expect_error(distrib_deriv4_y(poisson_distrib(), 2, list(mu = 3)))
})
