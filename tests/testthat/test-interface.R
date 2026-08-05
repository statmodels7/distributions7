# Interface contract: theta normalization, argument validation, printing,
# and class validators.

test_that("named theta is aligned to the distribution's parameter order", {
  d <- gaussian1_distrib()

  # Same value regardless of the order of a *named* theta
  expect_equal(
    distrib_pdf(d, 0, list(sigma = 2, mu = 1)),
    distrib_pdf(d, 0, list(mu = 1, sigma = 2))
  )
  expect_equal(
    distrib_gradient(d, 0.5, list(sigma = 2, mu = 1)),
    distrib_gradient(d, 0.5, list(mu = 1, sigma = 2))
  )

  # Named numeric vectors work too
  expect_equal(
    distrib_pdf(d, 0, c(sigma = 2, mu = 1)),
    distrib_pdf(d, 0, list(mu = 1, sigma = 2))
  )

  # Unnamed theta is interpreted positionally
  expect_equal(
    distrib_pdf(d, 0, list(1, 2)),
    distrib_pdf(d, 0, list(mu = 1, sigma = 2))
  )
})

test_that("malformed theta produces informative errors", {
  d <- gaussian1_distrib()

  expect_error(distrib_pdf(d, 0, list(mu = 0)), "Missing parameter")
  expect_error(distrib_pdf(d, 0, list(0)), "2 parameter")
  expect_error(
    distrib_gradient(d, rnorm(10), list(mu = 1:3, sigma = 1)),
    "dimension mismatch"
  )
  expect_error(
    distrib_hessian(d, rnorm(10), list(mu = 1:3, sigma = 1)),
    "dimension mismatch"
  )
})

test_that("derivatives recycle a scalar y over vectorized theta", {
  d <- gaussian1_distrib()
  g <- distrib_gradient(d, 0, list(mu = 1:5, sigma = 1))
  expect_length(g$mu, 5)
  expect_equal(g$mu, (0 - 1:5) / 1)
})

test_that("print method renders all distributions without error", {
  for (nm in names(all_distrib_cases())) {
    d <- all_distrib_cases()[[nm]]$d
    out <- capture.output(print(d))
    expect_true(any(grepl("^Distribution:", out)), label = nm)
    expect_true(any(grepl("Link:", out)), label = nm)
  }
})

test_that("distrib validator rejects inconsistent objects", {
  # n_params mismatch
  expect_error(
    distrib(
      distrib_name = "x", dimension = "univariate", bounds = c(0, 1),
      params = c("a", "b"), params_interpretation = c(a = "x", b = "y"),
      n_params = 5, params_bounds = list(a = c(0, 1), b = c(0, 1)),
      link_params = list(
        a = linkfunctions7::identity_link(),
        b = linkfunctions7::identity_link()
      )
    ),
    "n_params"
  )

  # missing params_bounds / link_params entries
  expect_error(
    distrib(
      distrib_name = "x", dimension = "univariate", bounds = c(0, 1),
      params = c("a", "b"), params_interpretation = c(a = "x", b = "y"),
      n_params = 2, params_bounds = list(), link_params = list()
    ),
    "params_bounds"
  )

  # non-link objects in link_params
  expect_error(gaussian1_distrib(link_mu = "not a link"), "link")

  # invalid dimension / bounds
  d <- gamma2_distrib()
  expect_error(d@dimension <- "bivariate", "univariate")
})

test_that("expectation() handles probability mass far from the origin", {
  d <- gaussian1_distrib()
  expect_equal(expectation(d, function(y, theta) y, list(mu = 200, sigma = 1)), 200, tolerance = 1e-6)
  expect_equal(expectation(d, function(y, theta) y, list(mu = -500, sigma = 1)), -500, tolerance = 1e-6)
  expect_equal(
    expectation(d, function(y, theta) (y - theta$mu)^2, list(mu = 50, sigma = 3)),
    9,
    tolerance = 1e-6
  )
})

test_that("expectation() works on heavy-tailed and bounded supports", {
  # Cauchy: E[1/(1+Y^2)] with mu=0, sigma=1 is exactly 1/2
  expect_equal(
    expectation(cauchy_distrib(), function(y, theta) 1 / (1 + y^2), list(mu = 0, sigma = 1)),
    0.5,
    tolerance = 1e-6
  )
  # Student t location
  expect_equal(
    expectation(student_t1_distrib(), function(y, theta) y, list(mu = 100, sigma = 1.3, nu = 6)),
    100,
    tolerance = 1e-5
  )
  # Beta on (0, 1)
  expect_equal(
    expectation(beta1_distrib(), function(y, theta) y, list(mu = 0.4, phi = 6)),
    0.4,
    tolerance = 1e-6
  )
  # Lognormal: E[Y] = exp(mu + sigma2/2)
  expect_equal(
    expectation(lognormal1_distrib(), function(y, theta) y, list(mu = 0.5, sigma2 = 1.3)),
    exp(0.5 + 1.3 / 2),
    tolerance = 1e-6
  )
})

test_that("expectation() is vectorized over theta", {
  expect_equal(
    expectation(poisson_distrib(), function(y, theta) y, list(mu = c(2, 4, 8))),
    c(2, 4, 8),
    tolerance = 1e-8
  )
  expect_equal(
    expectation(gaussian1_distrib(), function(y, theta) y, list(mu = c(0, 100), sigma = 2)),
    c(0, 100),
    tolerance = 1e-6
  )
})
