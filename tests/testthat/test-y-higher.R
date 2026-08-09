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

test_that("a family on a half line takes the fallback and it is finite", {
  # the location and the response are unrelated directions there, so the
  # identity does not apply and the stencil is what answers
  for (cs in list(list(gamma1_distrib(), list(mu = 2, phi = 0.6), c(0.7, 1.5, 3)),
                  list(weibull1_distrib(), list(mu = 2, sigma = 1.5), c(0.7, 1.5, 3)),
                  list(beta1_distrib(), list(mu = 0.4, phi = 6), c(0.2, 0.5, 0.8)))) {
    d <- cs[[1]]; th <- cs[[2]]; y <- cs[[3]]
    m <- S7::method(distrib_deriv3_y, S7::S7_class(d))
    expect_identical(attr(attr(m, "signature")[[1L]], "name"),
                     "continuous_distrib")
    expect_true(all(is.finite(distrib_deriv3_y(d, y, th))))
    expect_true(all(is.finite(distrib_deriv4_y(d, y, th))))
  }
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
