# The von Mises: the first family here whose support is a circle, so the two
# ends of the interval are the same point and the density need not vanish at
# either.

test_that("the density is the formula and integrates to one", {
  d <- vonmises1_distrib()
  th <- list(mu = 0.5, kappa = 2)
  y <- c(-2, -0.3, 0.5, 1.8, 3)

  expect_equal(distrib_pdf(d, y, th),
               exp(2 * cos(y - 0.5)) / (2 * pi * besselI(2, 0)))
  expect_equal(stats::integrate(function(z) distrib_pdf(d, z, th),
                                -pi, pi)$value, 1, tolerance = 1e-10)

  # the density does not vanish at the ends: they are the same point
  expect_gt(distrib_pdf(d, -pi, th), 0)
  expect_identical(d@bounds, c(-pi, pi))
})


test_that("a large concentration does not overflow the Bessel constant", {
  # I_0(kappa) itself is infinite past about 700; the scaled form with the
  # exponent added back is finite for any concentration.
  d <- vonmises1_distrib()
  expect_true(is.infinite(besselI(900, 0)))
  expect_true(is.finite(distrib_pdf(d, 0.5, list(mu = 0.5, kappa = 900))))
  expect_true(is.finite(distrib_pdf(d, 0.5, list(mu = 0.5, kappa = 5000),
                                    log = TRUE)))
})


test_that("the analytical derivatives match one Richardson differentiation", {
  skip_if_not_installed("numDeriv")
  d <- vonmises1_distrib()
  th <- list(mu = 0.5, kappa = 2)
  p0 <- c(0.5, 2)
  at <- function(p) list(mu = p[1], kappa = p[2])
  set.seed(1)
  y <- distrib_rng(d, 60, th)

  g <- distrib_gradient(d, y, th)
  expect_equal(vapply(g, sum, numeric(1)),
    numDeriv::grad(function(p) sum(distrib_pdf(d, y, at(p), log = TRUE)), p0),
    tolerance = 1e-7, ignore_attr = TRUE)

  J <- numDeriv::jacobian(function(p) {
    vapply(distrib_gradient(d, y, at(p)), sum, numeric(1))
  }, p0)
  h <- distrib_hessian(d, y, th)
  expect_equal(sum(h$mu_mu), J[1, 1], tolerance = 1e-7)
  expect_equal(sum(h$mu_kappa), J[1, 2], tolerance = 1e-7)
  expect_equal(sum(h$kappa_kappa), J[2, 2], tolerance = 1e-7)
})


test_that("the direction and the concentration are orthogonal", {
  # E[sin(Y - mu)] = 0 by symmetry, so the expected information is diagonal
  # and Fisher scoring updates the two independently.
  d <- vonmises1_distrib()
  for (th in list(list(mu = 0.5, kappa = 2), list(mu = -1, kappa = 0.3),
                  list(mu = 2, kappa = 12))) {
    eh <- distrib_expected_hessian(d, 0, th)
    expect_identical(eh$mu_kappa[1], 0)

    # and the diagonal against quadrature, which shares no code with it
    for (nm in c("mu_mu", "kappa_kappa")) {
      q <- stats::integrate(function(z) {
        distrib_hessian(d, z, th)[[nm]] * distrib_pdf(d, z, th)
      }, -pi, pi, rel.tol = 1e-11)$value
      expect_equal(eh[[nm]][1], q, tolerance = 1e-9, label = nm)
    }

    # A'(kappa) is the variance of cos(Y - mu), hence positive
    expect_gt(numericals7::bessel_i_ratio_derivs(th$kappa)$d1, 0)
  }
})


test_that("the generator reproduces the direction and the resultant length", {
  d <- vonmises1_distrib()
  th <- list(mu = 0.5, kappa = 2)
  set.seed(2)
  y <- distrib_rng(d, 2e5, th)

  expect_true(all(y >= -pi & y < pi + 1e-9))
  expect_equal(atan2(mean(sin(y)), mean(cos(y))), th$mu, tolerance = 0.02)
  expect_equal(sqrt(mean(cos(y))^2 + mean(sin(y))^2), numericals7::bessel_i_ratio(th$kappa),
               tolerance = 0.01)
})


test_that("the validator passes and a fit recovers the parameters", {
  d <- vonmises1_distrib()
  set.seed(3)
  res <- check_distrib(d, verbose = FALSE)
  expect_true(all(res$status == "OK"),
    label = paste(res$check[res$status != "OK"], collapse = ", "))

  set.seed(5)
  th <- list(mu = 0.5, kappa = 2)
  y <- distrib_rng(d, 3000, th)
  f <- fit_distrib(d, y)
  expect_true(f@converged, info = fit_report(f, d, th))
  expect_equal(unname(coef(f)), c(0.5, 2), tolerance = 0.15)
})
