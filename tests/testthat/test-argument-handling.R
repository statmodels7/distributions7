# Argument handling shared by every generic: empty input, parameter recycling,
# and parameter values that carry names.

test_that("an empty y gives empty derivatives rather than an error", {
  # distrib_pdf has always returned numeric(0) here; the derivative generics used
  # to reject it with "'y' must have length 1 or 1, not 0".
  d <- gaussian1_distrib()
  th <- list(mu = 1.5, sigma = 2)

  expect_length(distrib_pdf(d, numeric(0), th), 0)
  for (fn in list(distrib_gradient, distrib_hessian, distrib_deriv3, distrib_deriv4)) {
    g <- fn(d, numeric(0), th)
    expect_type(g, "list")
    expect_true(all(lengths(g) == 0))
  }
})

test_that("rng recycles per-observation parameters up to n", {
  # rnorm(4, c(0, 10)) draws two from each mean; the package generators follow
  # the same convention.
  set.seed(21)
  d <- gaussian1_distrib()
  y <- distrib_rng(d, 4000, list(mu = rep(c(0, 50), each = 2000), sigma = 1))
  expect_length(y, 4000)
  expect_equal(mean(y[1:2000]), 0, tolerance = 0.15)
  expect_equal(mean(y[2001:4000]), 50, tolerance = 0.15)

  # a parameter shorter than n is recycled, not rejected
  set.seed(22)
  ph <- pseudohuber_distrib()
  expect_length(
    distrib_rng(ph, 6, list(mu = c(0, 10), sigma = 1, nu = 2)),
    6
  )

  # the zero-adjusted wrapper used to subset a short parameter with a longer
  # logical index, which silently produced NA parameters
  set.seed(23)
  za <- zero_adjusted(gamma2_distrib())
  r <- distrib_rng(za, 8, list(mu = c(2, 5), sigma2 = 1, za = 0.3))
  expect_length(r, 8)
  expect_true(all(is.finite(r)))
})

test_that("parameters carrying names are handled like bare numbers", {
  # Regression: theta elements built through a link function inherit the name of
  # the parameter. transpose_params keyed rows by those inner names, so a
  # three-parameter theta collapsed to a single column and every consumer of
  # expectation() -- the numerical cdf, the Monte Carlo strategy, Fisher scoring
  # -- saw a theta with only its first parameter.
  d <- pseudohuber_distrib()
  bare <- list(mu = 0.5, sigma = 1.4, nu = 2.5)
  named <- list(mu = c(mu = 0.5), sigma = c(sigma = 1.4), nu = c(nu = 2.5))

  expect_equal(
    unname(transpose_params(named)[[1]]),
    unname(transpose_params(bare)[[1]])
  )
  expect_named(transpose_params(named)[[1]], c("mu", "sigma", "nu"))

  expect_equal(distrib_pdf(d, 1, named), distrib_pdf(d, 1, bare))
  expect_equal(
    distrib_expected_hessian(d, 0, named)$mu_mu,
    distrib_expected_hessian(d, 0, bare)$mu_mu
  )
})

test_that("fit_distrib survives a starting value where the quadrature fails", {
  # Fisher scoring on a distribution whose expected Hessian is approximated
  # numerically can hit "the integral is probably divergent" at an awkward
  # parameter value. That must count as a failed start, not abort the fit.
  skip_on_cran()
  d <- pseudohuber_distrib()
  th <- list(mu = 0.5, sigma = 1.4, nu = 2.5)
  set.seed(24)
  y <- distrib_rng(d, 300, th)

  for (s in 1:4) {
    set.seed(s)
    f <- expect_no_error(fit_distrib(d, y))
    expect_type(stats::coef(f), "double")
    expect_named(stats::coef(f), c("mu", "sigma", "nu"))
  }
})
