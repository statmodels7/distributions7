# The closed mixed response-parameter derivatives.
#
# Each is checked against Richardson extrapolation of the analytic response
# derivative in the parameter -- a route that shares no code with the closed
# forms, which differentiate the response derivative by hand. Richardson
# rather than a plain central difference, because several of these families
# reach the closed form through quantities that are themselves differences.

cross_ref <- function(d, y, th, nm) {
  j <- match(nm, d@params)
  vapply(seq_along(y), function(i) {
    numDeriv::grad(function(z) {
      t2 <- th
      t2[[j]] <- z
      distrib_grad_y(d, y[i], t2)
    }, th[[j]])
  }, numeric(1))
}

expect_cross <- function(d, th, y, tol, only = NULL) {
  got <- distrib_cross_y(d, y, th)
  for (nm in if (is.null(only)) names(got) else only) {
    ref <- cross_ref(d, y, th, nm)
    # the scale is taken over the whole component and not point by point: the
    # generalized Pareto's shape component passes through an exact zero at
    # y = sigma, and there no pointwise relative measure exists -- the
    # reference is its own truncation error. That zero is asserted separately
    # below, as the structural fact it is.
    den <- max(abs(ref), abs(got[[nm]]), 1e-8)
    expect_lt(max(abs(got[[nm]] - ref)) / den, tol)
  }
}


test_that("the location-scale identity gives the mixed block", {
  skip_if_not_installed("numDeriv")
  y <- c(-1.3, 0.4, 2.2)
  expect_cross(logistic_distrib(), list(mu = 0.3, sigma = 1.4), y, 1e-7)
  expect_cross(cauchy_distrib(), list(mu = 0.3, sigma = 1.4), y, 1e-7)
  expect_cross(gumbel_distrib(), list(mu = 0.3, sigma = 1.4), y, 1e-7)
  # the Laplace away from its kink, where the quantity does not exist
  expect_cross(laplace_distrib(), list(mu = 0.3, sigma = 1.4), y, 1e-9)
})


test_that("the location-scale block closes for the shape families too", {
  skip_if_not_installed("numDeriv")
  y <- c(-1.3, 0.4, 2.2)
  expect_cross(pseudohuber_distrib(), list(mu = 0.2, sigma = 1.1, nu = 2), y, 1e-7)
  expect_cross(skewnormal1_distrib(),
               list(mu = 0.2, sigma = 1.2, alpha = 1.5), y, 1e-7)
  expect_cross(skewt_distrib(),
               list(mu = 0.2, sigma = 1.2, alpha = 1.0, nu = 6), y, 1e-6)
})


test_that("the shape component of two of them is closed rather than differenced", {
  # The values agree with Richardson either way, so a tolerance cannot tell a
  # closed form from the difference it replaced. What separates them is
  # whether the body reaches numerical_cross_y at all.
  differences <- function(cls) {
    m <- S7::method(distrib_cross_y, cls)
    any(grepl("numerical_cross_y", deparse(body(m)), fixed = TRUE))
  }
  expect_false(differences(distributions7:::SkewNormal1Distrib))
  expect_false(differences(distributions7:::PseudoHuberDistrib))
  # and the skew t is not claimed to be closed: its nu direction carries the
  # derivative of a Student t distribution function in its degrees of freedom
  expect_true(differences(distributions7:::SkewTDistrib))
})


test_that("the closed shape components hold far into the tail", {
  skip_if_not_installed("numDeriv")
  # where phi and Phi both underflow, so that the ratio is carried on the log
  # scale, and where D is dominated by r^2/sigma^2
  d <- skewnormal1_distrib()
  th <- list(mu = 1, sigma = 2, alpha = 3)
  y <- th$mu + th$sigma * c(-30, -12, 12)
  expect_cross(d, th, y, 1e-6, only = "alpha")
  expect_true(all(is.finite(distrib_cross_y(d, y, th)$alpha)))

  d2 <- pseudohuber_distrib()
  th2 <- list(mu = 0.5, sigma = 1.3, nu = 2.5)
  y2 <- th2$mu + th2$sigma * c(-500, -20, 20, 500)
  expect_cross(d2, th2, y2, 1e-6, only = "nu")
})


test_that("the scale families close both parameters", {
  skip_if_not_installed("numDeriv")
  y <- c(0.3, 1.2, 3.1)
  expect_cross(exponential_distrib(), list(mu = 1.6), y, 1e-8)
  expect_cross(weibull1_distrib(), list(mu = 1.4, sigma = 1.9), y, 1e-8)
  expect_cross(gpd_distrib(), list(sigma = 1.2, xi = 0.35), y, 1e-7)
  expect_cross(gpd_distrib(), list(sigma = 2.0, xi = -0.25), y, 1e-7)
  expect_cross(gengamma1_distrib(), list(a = 1.3, d = 2.1, p = 1.6), y, 1e-8)
})


test_that("the generalized Pareto's shape component vanishes exactly at y = sigma", {
  # d^2 l / dy dxi = -1/(sigma t) + (xi+1) y /(sigma t)^2 is zero when
  # (xi+1) y = sigma + xi y, that is at y = sigma whatever the shape. A
  # structural zero, so it is asserted as one rather than to a tolerance.
  for (xi in c(-0.4, 0.01, 0.35, 2)) {
    got <- distrib_cross_y(gpd_distrib(), 1.2, list(sigma = 1.2, xi = xi))
    expect_lt(abs(got$xi), 1e-15)
  }
})


test_that("the response derivatives newly written match the log-density", {
  skip_if_not_installed("numDeriv")
  cases <- list(
    list(gpd_distrib(), list(sigma = 1.2, xi = 0.35), c(0.3, 1.2, 3.1)),
    list(gengamma1_distrib(), list(a = 1.3, d = 2.1, p = 1.6), c(0.3, 1.2, 3.1)),
    list(invgauss2_distrib(), list(mu = 1.3, lambda = 2.2), c(0.3, 1.2, 3.1)),
    list(vonmises2_distrib(), list(mu = 0.3, rho = 0.6), c(-2, 0.5, 2.4))
  )
  for (cs in cases) {
    d <- cs[[1]]; th <- cs[[2]]; y <- cs[[3]]
    ref <- vapply(y, function(v) {
      numDeriv::grad(function(z) distrib_pdf(d, z, th, log = TRUE), v)
    }, numeric(1))
    expect_equal(distrib_grad_y(d, y, th), ref, tolerance = 1e-7)
    refh <- vapply(y, function(v) {
      numDeriv::grad(function(z) distrib_grad_y(d, z, th), v)
    }, numeric(1))
    expect_equal(distrib_hess_y(d, y, th), refh, tolerance = 1e-7)
  }
})


test_that("a reparametrized family chains the parent's mixed block", {
  skip_if_not_installed("numDeriv")
  y <- c(0.3, 1.2, 3.1)
  expect_cross(invgauss2_distrib(), list(mu = 1.3, lambda = 2.2), y, 1e-7)
  expect_cross(lognormal2_distrib(), list(mean = 1.4, var = 0.7), y, 1e-6)
  expect_cross(weibull3_distrib(), list(mean = 1.5, sigma = 1.8), y, 1e-8)
  expect_cross(student_t2_distrib(), list(mu = 0.2, sigma = 1.2, nu = 7),
               c(-1.3, 0.4, 2.2), 1e-7)
  expect_cross(vonmises2_distrib(), list(mu = 0.3, rho = 0.6),
               c(-2, 0.5, 2.4), 1e-7)
})


test_that("the families whose score in y is short close it directly", {
  skip_if_not_installed("numDeriv")
  y <- c(-1.3, 0.4, 2.2)
  yp <- c(0.3, 1.2, 3.1)
  yu <- c(0.15, 0.5, 0.82)
  expect_cross(vonmises1_distrib(), list(mu = 0.3, kappa = 1.7),
               c(-2, 0.5, 2.4), 1e-7)
  expect_cross(gengamma1_distrib(), list(a = 1.3, d = 2.1, p = 1.6), yp, 1e-8)
  expect_cross(gaussian2_distrib(), list(mu = 0.5, sigma2 = 1.7), y, 1e-8)
  expect_cross(gaussian3_distrib(), list(mu = 0.5, tau = 0.8), y, 1e-8)
  expect_cross(lognormal1_distrib(), list(mu = 0.2, sigma2 = 0.6), yp, 1e-7)
  expect_cross(invgauss1_distrib(), list(mu = 1.3, phi = 0.5), yp, 1e-8)
  expect_cross(beta1_distrib(), list(mu = 0.4, phi = 3), yu, 1e-8)
  expect_cross(beta2_distrib(), list(alpha = 2, beta = 3), yu, 1e-8)
  expect_cross(chisq_distrib(), list(mu = 4), yp, 1e-9)
  expect_cross(gamma1_distrib(), list(mu = 1.5, phi = 0.6), yp, 1e-8)
  expect_cross(gamma2_distrib(), list(mu = 1.5, sigma2 = 0.9), yp, 1e-7)
  expect_cross(skewnormal2_distrib(),
               list(mu = 0.3, sigma = 1.2, gamma1 = 0.4), y, 1e-7)
})


test_that("every continuous family has a closed mixed block", {
  # The audit this batch was built from, run as a test so that a family added
  # later without one is caught. A method is the base fallback when the class
  # it was registered on is continuous_distrib itself.
  exported <- grep("_distrib$", getNamespaceExports("distributions7"), value = TRUE)
  skip_these <- c("mvgaussian1_distrib", "mvstudent_t1_distrib",
                  "dirichlet_distrib", "multinomial_distrib",
                  "check_distrib", "continuous_distrib", "discrete_distrib",
                  "fit_distrib", "multivariate_distrib")
  open <- character()
  for (ctor in setdiff(exported, skip_these)) {
    d <- tryCatch(get(ctor, asNamespace("distributions7"))(),
                  error = function(e) NULL)
    if (is.null(d) || !S7::S7_inherits(d, distributions7:::continuous_distrib)) next
    m <- S7::method(distrib_cross_y, S7::S7_class(d))
    owner <- attr(attr(m, "signature")[[1L]], "name")
    if (identical(owner, "continuous_distrib")) open <- c(open, ctor)
  }
  expect_identical(open, character(),
                   info = paste("still on the fallback:", paste(open, collapse = ", ")))
})
