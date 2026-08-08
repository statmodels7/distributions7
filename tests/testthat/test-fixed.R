# fixed(): parameters held at known values.
#
# The wrapper computes nothing of its own -- every method splices the fixed
# values into theta and delegates -- so what the tests establish is that the
# splicing and the subsetting are exact: each answer equals the parent's at the
# full parameter vector, and each derivative equals the corresponding component
# of the parent's. Components are always extracted with [[ ]], never $, which
# partially matches.

test_that("the constructor validates its arguments", {
  g <- gaussian1_distrib()

  expect_error(fixed(g), "at least one named value")
  expect_error(fixed(g, nu = 1), "not a parameter")
  expect_error(fixed(g, 0.5), "must be named")
  expect_error(fixed(g, mu = 0, mu = 1), "more than once")
  expect_error(fixed(g, sigma = c(1, 2)), "single finite number")
  expect_error(fixed(g, sigma = Inf), "single finite number")
  expect_error(fixed(g, sigma = "a"), "single finite number")

  # The domains are open, so a value exactly on a bound is refused.
  expect_error(fixed(g, sigma = 0), "open domain")
  expect_error(fixed(beta1_distrib(), mu = 1), "open domain")
})

test_that("the free parameter set is the parent's minus the fixed ones", {
  d <- fixed(gaussian1_distrib(), mu = 0.5)

  expect_identical(d@params, "sigma")
  expect_identical(d@n_params, 1L)
  expect_named(d@fixed_params, "mu")
  expect_identical(d@params_bounds, list(sigma = c(0, Inf)))
  expect_named(d@link_params, "sigma")
  expect_true(S7::S7_inherits(d, FixedContinuousDistrib))
  expect_true(S7::S7_inherits(fixed(poisson_distrib(), mu = 2), FixedDiscreteDistrib))
})

test_that("fixing a fixed distribution collapses into one wrapper", {
  g <- gaussian1_distrib()
  d <- fixed(fixed(g, mu = 0), sigma = 1)

  expect_true(S7::S7_inherits(d@parent_distrib, Gaussian1Distrib))
  expect_identical(names(d@fixed_params), c("mu", "sigma"))
  expect_identical(d@n_params, 0L)

  # A parameter that is no longer free cannot be fixed again.
  expect_error(fixed(fixed(g, mu = 0), mu = 1), "not a free parameter")
})

test_that("pdf, cdf, quantile and rng equal the parent's at the full theta", {
  g <- gaussian1_distrib()
  d <- fixed(g, mu = 0.5)
  th <- list(sigma = 2)
  full <- list(mu = 0.5, sigma = 2)
  y <- c(-1.5, 0.2, 3)

  expect_equal(distrib_pdf(d, y, th), distrib_pdf(g, y, full))
  expect_equal(
    distrib_pdf(d, y, th, log = TRUE),
    distrib_pdf(g, y, full, log = TRUE)
  )
  expect_equal(distrib_cdf(d, y, th), distrib_cdf(g, y, full))
  expect_equal(
    distrib_quantile(d, c(0.1, 0.5, 0.9), th),
    distrib_quantile(g, c(0.1, 0.5, 0.9), full)
  )

  set.seed(11)
  r_fixed <- distrib_rng(d, 20, th)
  set.seed(11)
  r_parent <- distrib_rng(g, 20, full)
  expect_identical(r_fixed, r_parent)
})

test_that("theta stays vectorized in the free parameters", {
  d <- fixed(gaussian1_distrib(), mu = 0)
  y <- c(-1, 0, 2)
  sig <- c(1, 2, 3)

  expect_equal(
    distrib_pdf(d, y, list(sigma = sig)),
    dnorm(y, 0, sig)
  )
})

test_that("the derivatives are the parent's components among the free parameters", {
  g <- gaussian1_distrib()
  d <- fixed(g, mu = 0.5)
  th <- list(sigma = 2)
  full <- list(mu = 0.5, sigma = 2)
  y <- c(-1.5, 0.2, 3)

  gr <- distrib_gradient(d, y, th)
  expect_named(gr, "sigma")
  expect_equal(gr[["sigma"]], distrib_gradient(g, y, full)[["sigma"]])

  he <- distrib_hessian(d, y, th)
  expect_named(he, "sigma_sigma")
  expect_equal(he[["sigma_sigma"]], distrib_hessian(g, y, full)[["sigma_sigma"]])

  d3 <- distrib_deriv3(d, y, th)
  expect_named(d3, "sigma_sigma_sigma")
  expect_equal(
    d3[["sigma_sigma_sigma"]],
    distrib_deriv3(g, y, full)[["sigma_sigma_sigma"]]
  )

  d4 <- distrib_deriv4(d, y, th)
  expect_named(d4, "sigma_sigma_sigma_sigma")
  expect_equal(
    d4[["sigma_sigma_sigma_sigma"]],
    distrib_deriv4(g, y, full)[["sigma_sigma_sigma_sigma"]]
  )

  eh <- distrib_expected_hessian(d, y, th)
  expect_equal(
    eh[["sigma_sigma"]],
    distrib_expected_hessian(g, y, full)[["sigma_sigma"]]
  )

  expect_equal(distrib_grad_y(d, y, th), distrib_grad_y(g, y, full))
  expect_equal(distrib_hess_y(d, y, th), distrib_hess_y(g, y, full))

  gc <- distrib_grad_cdf(d, y, th)
  expect_equal(gc[["sigma"]], distrib_grad_cdf(g, y, full)[["sigma"]])
  hc <- distrib_hess_cdf(d, y, th)
  expect_equal(hc[["sigma_sigma"]], distrib_hess_cdf(g, y, full)[["sigma_sigma"]])
})

test_that("subsetting by name survives a three-parameter parent", {
  # Fixing the middle parameter leaves components whose names are built from
  # non-adjacent parents; each must still line up with the parent's component
  # of the same name.
  p <- pseudohuber_distrib()
  expect_identical(p@params, c("mu", "sigma", "nu"))

  d <- fixed(p, sigma = 1.5)
  th <- list(mu = 0.3, nu = 2)
  full <- list(mu = 0.3, sigma = 1.5, nu = 2)
  y <- c(-1, 0.5, 2)

  gr <- distrib_gradient(d, y, th)
  grp <- distrib_gradient(p, y, full)
  expect_named(gr, c("mu", "nu"))
  expect_equal(gr[["mu"]], grp[["mu"]])
  expect_equal(gr[["nu"]], grp[["nu"]])

  he <- distrib_hessian(d, y, th)
  hep <- distrib_hessian(p, y, full)
  expect_named(he, c("mu_mu", "nu_nu", "mu_nu"))
  expect_equal(he[["mu_nu"]], hep[["mu_nu"]])

  d3 <- distrib_deriv3(d, y, th)
  d3p <- distrib_deriv3(p, y, full)
  expect_equal(d3[["mu_mu_nu"]], d3p[["mu_mu_nu"]])
})

test_that("the link scale is the parent's, restricted to the free parameters", {
  g <- gaussian1_distrib()
  d <- fixed(g, mu = 0.5)
  th <- list(sigma = 2)
  full <- list(mu = 0.5, sigma = 2)
  y <- c(-1, 0.2, 3)

  # The free parameter keeps its link, so the chain rule contributes the same
  # factor on both routes and the components must agree exactly.
  gl <- distrib_gradient(d, y, th, scale = "link")
  glp <- distrib_gradient(g, y, full, scale = "link")
  expect_equal(gl[["sigma"]], glp[["sigma"]])

  hl <- distrib_hessian(d, y, th, scale = "link")
  hlp <- distrib_hessian(g, y, full, scale = "link")
  expect_equal(hl[["sigma_sigma"]], hlp[["sigma_sigma"]])
})

test_that("moments delegate to the parent's closed forms", {
  d <- fixed(gaussian1_distrib(), mu = 0.5)
  th <- list(sigma = 2)

  expect_equal(mean(d, th), 0.5)
  expect_equal(variance(d, th), 4)
  expect_equal(std_dev(d, th), 2)
  expect_equal(skewness(d, th), 0)
})

test_that("params_smooth travels with the free parameters", {
  dl <- fixed(laplace_distrib(), sigma = 1)
  expect_false(param_smoothness(dl)[["mu"]])

  dm <- fixed(laplace_distrib(), mu = 0)
  expect_true(param_smoothness(dm)[["sigma"]])
})

test_that("fixing every parameter leaves a fully known distribution", {
  d0 <- fixed(gaussian1_distrib(), mu = 1, sigma = 2)

  expect_identical(d0@n_params, 0L)
  expect_identical(d0@params, character(0))
  expect_equal(distrib_pdf(d0, c(0, 1), list()), dnorm(c(0, 1), 1, 2))
  expect_equal(distrib_cdf(d0, 1, list()), 0.5)
  expect_output(print(d0), "none free")
})

test_that("a wrapper's own parameter can be fixed", {
  zip <- zero_inflated(poisson_distrib())
  fz <- fixed(zip, zi = 0.3)
  th <- list(mu = 2)
  full <- list(mu = 2, zi = 0.3)

  expect_identical(fz@params, "mu")
  expect_equal(distrib_pdf(fz, 0:5, th), distrib_pdf(zip, 0:5, full))
  expect_equal(
    distrib_gradient(fz, 0:5, th)[["mu"]],
    distrib_gradient(zip, 0:5, full)[["mu"]]
  )

  # The atoms of a mixed parent come through unchanged.
  za <- zero_adjusted(gamma2_distrib())
  fa <- fixed(za, za = 0.2)
  tha <- list(mu = 2, sigma2 = 1)
  fulla <- list(mu = 2, sigma2 = 1, za = 0.2)
  expect_equal(distrib_atoms(fa, tha), distrib_atoms(za, fulla))
})

test_that("fixed composes with truncation in both orders", {
  g <- gaussian1_distrib()
  full <- list(mu = 0, sigma = 1)
  y <- c(0.5, 1, 2)

  a <- truncated(fixed(g, mu = 0), lower = 0)
  b <- fixed(truncated(g, lower = 0), mu = 0)
  expect_equal(
    distrib_pdf(a, y, list(sigma = 1)),
    distrib_pdf(b, y, list(sigma = 1))
  )
  expect_equal(
    distrib_pdf(a, y, list(sigma = 1)),
    distrib_pdf(truncated(g, lower = 0), y, full)
  )
})

test_that("check_distrib passes on a fixed distribution", {
  set.seed(21)
  res <- check_distrib(fixed(gaussian1_distrib(), mu = 0.3),
    theta = list(sigma = 1.5), orders = 1:2, nsim = 2e4, verbose = FALSE
  )
  expect_true(all(res$status == "OK"))
})

test_that("fit_distrib estimates only the free parameters", {
  set.seed(31)
  y <- rnorm(500, mean = 0, sd = 1.7)
  d <- fixed(gaussian1_distrib(), mu = 0)

  fit <- fit_distrib(d, y)
  est <- coef(fit)

  expect_named(est, "sigma")
  # With mu known the MLE of sigma is sqrt(mean(y^2)), in closed form.
  expect_equal(unname(est[["sigma"]]), sqrt(mean(y^2)), tolerance = 1e-5)

  # And the score at the optimum vanishes.
  sc <- distrib_gradient(d, y, as.list(est))
  expect_lt(abs(sum(sc[["sigma"]])), 1e-4)
})

test_that("print shows the fixed values", {
  d <- fixed(gaussian1_distrib(), mu = 0.5)
  expect_output(print(d), "Fixed:")
  expect_output(print(d), "mu = 0.5")
  expect_output(print(d), "sigma")
})
