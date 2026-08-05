# The two second parametrisations that are families in their own right rather
# than reparametrize() wrappers, each for a reason that is mathematical and not
# a matter of taste: the von Mises because A has no elementary inverse, and the
# skew normal because its map carries a sign and because what makes the centred
# parametrisation worth having is a cancellation a generic chain rule would
# compute as a difference of large numbers.


test_that("the von Mises inverse is an inverse", {
  rho <- c(0.05, 0.2, 0.5, 0.7, 0.9, 0.99)
  k <- vm_kappa_of_rho(rho)$kappa
  expect_equal(vm_A(k), rho, tolerance = 1e-10)
  expect_true(all(diff(k) > 0))
  # outside (0, 1) there is no concentration to find
  expect_true(is.na(vm_kappa_of_rho(0)$kappa))
  expect_true(is.na(vm_kappa_of_rho(1)$kappa))
})


test_that("the derivatives of A come from the recurrences", {
  skip_if_not_installed("numDeriv")
  # Each order against ONE Richardson pass on the analytic order below, which
  # is the rule everywhere in this package.
  for (k in c(0.5, 2, 10)) {
    a <- vm_A_derivs(k)
    expect_equal(a$d1, numDeriv::grad(vm_A, k), tolerance = 1e-8)
    expect_equal(a$d2, numDeriv::grad(function(z) vm_A_derivs(z)$d1, k),
                 tolerance = 1e-7, label = paste("A''", k))
    expect_equal(a$d3, numDeriv::grad(function(z) vm_A_derivs(z)$d2, k),
                 tolerance = 1e-6, label = paste("A'''", k))
    expect_equal(a$d4, numDeriv::grad(function(z) vm_A_derivs(z)$d3, k),
                 tolerance = 1e-6, label = paste("A''''", k))
  }
  # A' is the variance of cos(Y - mu) and so is strictly positive, which is
  # what keeps the inverse function rule from dividing by zero
  expect_true(all(vm_A_derivs(c(0.01, 1, 100))$d1 > 0))
})


test_that("the derivatives of the von Mises inverse follow the inverse rule", {
  skip_if_not_installed("numDeriv")
  for (r in c(0.2, 0.5, 0.8)) {
    k <- vm_kappa_of_rho(r)
    expect_equal(k$d1, numDeriv::grad(function(z) vm_kappa_of_rho(z)$kappa, r),
                 tolerance = 1e-6, label = paste("k'", r))
    expect_equal(k$d2, numDeriv::grad(function(z) vm_kappa_of_rho(z)$d1, r),
                 tolerance = 1e-6, label = paste("k''", r))
    expect_equal(k$d3, numDeriv::grad(function(z) vm_kappa_of_rho(z)$d2, r),
                 tolerance = 1e-5, label = paste("k'''", r))
    expect_equal(k$d4, numDeriv::grad(function(z) vm_kappa_of_rho(z)$d3, r),
                 tolerance = 1e-4, label = paste("k''''", r))
    # the first derivative is the reciprocal of A', by construction
    expect_equal(k$d1, 1 / vm_A_derivs(k$kappa)$d1)
  }
})


test_that("vonmises2 is the von Mises at the matching concentration", {
  d <- vonmises2_distrib()
  th <- list(mu = 0.5, rho = 0.7)
  k <- vm_kappa_of_rho(0.7)$kappa
  set.seed(2)
  y <- distrib_rng(d, 200, th)
  expect_identical(distrib_pdf(d, y, th),
                   distrib_pdf(vonmises1_distrib(), y,
                               list(mu = 0.5, kappa = k)))

  set.seed(1)
  res <- check_distrib(d, th, verbose = FALSE)
  expect_identical(nrow(res), 13L)
  expect_true(all(res$status == "OK"),
              info = paste(res$check[res$status != "OK"], collapse = ", "))

  # the expected information in rho is the inverse of the one in kappa, which
  # is what a one-to-one change of a single parameter must give
  eh <- distrib_expected_hessian(d, 0, th)
  expect_equal(eh[["rho_rho"]][1], -1 / vm_A_derivs(k)$d1, tolerance = 1e-10)
  expect_identical(eh[["mu_rho"]][1], 0)

  set.seed(9)
  yb <- distrib_rng(d, 2000, th)
  f <- fit_distrib(d, yb)
  expect_true(f@converged)
  expect_equal(unname(coef(f)), c(0.5, 0.7), tolerance = 0.1)
})


test_that("the skew normal's centred parameters ARE its moments", {
  d <- skewnormal2_distrib()
  for (g1 in c(-0.9, -0.3, 0.3, 0.9)) {
    th <- list(mu = 2, sigma = 3, gamma1 = g1)
    expect_equal(mean(d, th), 2, tolerance = 1e-10, label = as.character(g1))
    expect_equal(sqrt(variance(d, th)), 3, tolerance = 1e-10,
                 label = as.character(g1))
    expect_equal(skewness(d, th), g1, tolerance = 1e-10,
                 label = as.character(g1))
  }
  # the ceiling, which is also the reason the skew t exists
  expect_equal(sn_max_skew(), 0.9952717, tolerance = 1e-6)
  expect_identical(unname(skewnormal2_distrib()@params_bounds$gamma1),
                   c(-sn_max_skew(), sn_max_skew()))
})


test_that("skewnormal2 passes the validator on both signs", {
  d <- skewnormal2_distrib()
  for (g1 in c(-0.7, -0.2, 0.2, 0.7)) {
    th <- list(mu = 0, sigma = 1, gamma1 = g1)
    set.seed(1)
    res <- check_distrib(d, th, verbose = FALSE)
    expect_identical(nrow(res), 13L, label = as.character(g1))
    expect_true(all(res$status == "OK"),
                info = paste(g1, paste(res$check[res$status != "OK"],
                                       collapse = ", ")))
  }
})


test_that("the centred information is non-singular where the direct one is not", {
  # This is the property the parametrisation exists for, and the comparison is
  # the whole point: in the direct parametrisation the score for alpha is
  # proportional to the score for the location at alpha = 0, so the information
  # loses a rank. In the centred one it does not.
  info <- function(d, theta) {
    eh <- distrib_expected_hessian(d, 0, theta)
    pr <- hess_pairs(d@params)
    I <- matrix(0, d@n_params, d@n_params)
    for (nm in names(eh)) {
      k <- match(pr[[nm]], d@params)
      I[k[1], k[2]] <- I[k[2], k[1]] <- -eh[[nm]][1]
    }
    I
  }

  d2 <- skewnormal2_distrib()
  dets <- vapply(c(0.5, 0.1, 0.01, 1e-3), function(g) {
    det(info(d2, list(mu = 0, sigma = 1, gamma1 = g)))
  }, numeric(1))
  # it settles rather than collapsing
  expect_true(all(dets > 0.2))
  expect_lt(abs(dets[4] - dets[3]) / dets[3], 0.05)

  d1 <- skewnormal1_distrib()
  dets1 <- vapply(c(1, 0.1, 0.01), function(a) {
    det(info(d1, list(mu = 0, sigma = 1, alpha = a)))
  }, numeric(1))
  # and the direct one falls like the fourth power of alpha
  expect_true(all(diff(dets1) < 0))
  expect_lt(dets1[3], 1e-8)
  expect_gt(dets1[1] / dets1[3], 1e6)
})


test_that("the score does not follow the map that diverges", {
  # The cancellation, stated as a number: the Jacobian grows without bound as
  # the skewness goes to zero while the information in gamma1 does not move.
  d <- skewnormal2_distrib()
  dalpha <- function(g, h = 1e-6) {
    (sn2_theta(list(mu = 0, sigma = 1, gamma1 = g + h))$alpha -
       sn2_theta(list(mu = 0, sigma = 1, gamma1 = g - h))$alpha) / (2 * h)
  }
  expect_gt(dalpha(1e-4) / dalpha(0.5), 50)

  i_gamma <- function(g) {
    -distrib_expected_hessian(d, 0, list(mu = 0, sigma = 1,
                                         gamma1 = g))[["gamma1_gamma1"]][1]
  }
  expect_equal(i_gamma(0.01), i_gamma(0.05), tolerance = 0.05)
  expect_true(is.finite(i_gamma(1e-3)))
})


test_that("skewnormal2 agrees with the direct family at matching parameters", {
  d2 <- skewnormal2_distrib()
  d1 <- skewnormal1_distrib()
  th2 <- list(mu = 1, sigma = 2, gamma1 = -0.6)
  th1 <- sn2_theta(th2)
  set.seed(3)
  y <- distrib_rng(d2, 60, th2)
  expect_identical(distrib_pdf(d2, y, th2), distrib_pdf(d1, y, th1))
  expect_identical(distrib_cdf(d2, y, th2), distrib_cdf(d1, y, th1))
  # the response derivatives are the parent's, the coordinates having changed
  # and the response not
  expect_identical(distrib_grad_y(d2, y, th2), distrib_grad_y(d1, y, th1))
})


test_that("skewnormal2 fits from both sides", {
  d <- skewnormal2_distrib()
  for (g1 in c(-0.6, 0.6)) {
    th <- list(mu = 2, sigma = 3, gamma1 = g1)
    set.seed(11)
    y <- distrib_rng(d, 3000, th)
    f <- fit_distrib(d, y)
    expect_true(f@converged, label = as.character(g1))
    expect_equal(unname(coef(f)), unname(unlist(th)), tolerance = 0.15,
                 label = as.character(g1))
  }
})


test_that("the skew normal gradient agrees with Richardson", {
  skip_if_not_installed("numDeriv")
  d <- skewnormal2_distrib()
  th <- list(mu = 0, sigma = 1, gamma1 = 0.5)
  set.seed(2)
  y <- distrib_rng(d, 60, th)
  q0 <- unlist(th)
  llv <- function(q) {
    sum(distrib_pdf(d, y, as.list(stats::setNames(q, d@params)), log = TRUE))
  }
  g <- vapply(distrib_gradient(d, y, th), sum, numeric(1))
  expect_equal(unname(g), numDeriv::grad(llv, q0), tolerance = 1e-6)
})
