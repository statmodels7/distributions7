# Variable transformations: change of variables, delegation of derivatives,
# and support handling.

test_that("exp(gaussian) reproduces the lognormal distribution exactly", {
  d <- transformation(gaussian_distrib(), exp_transform())
  th <- list(mu = 0.4, sigma = 1.1)
  y <- c(0.3, 1, 2.5, 7)
  p <- c(0.1, 0.5, 0.9)

  expect_equal(distrib_pdf(d, y, th), dlnorm(y, 0.4, 1.1))
  expect_equal(distrib_pdf(d, y, th, log = TRUE), dlnorm(y, 0.4, 1.1, log = TRUE))
  expect_equal(distrib_cdf(d, y, th), plnorm(y, 0.4, 1.1))
  expect_equal(distrib_quantile(d, p, th), qlnorm(p, 0.4, 1.1))

  # Derivatives coincide with the equivalent native lognormal parameterization
  logn <- transformation(gaussian_distrib(), exp_transform())
  g_t <- distrib_gradient(logn, y, th)
  g_n <- fd_gradient_ref(logn, y, th)
  expect_equal(g_t$mu, g_n$mu, tolerance = 1e-5)
  expect_equal(g_t$sigma, g_n$sigma, tolerance = 1e-5)
})

test_that("all transformers produce correctly normalized densities", {
  thg <- list(mu = 0.4, sigma = 1.1)
  thga <- list(mu = 3, sigma2 = 2)
  cases <- list(
    log_gamma = transformation(gamma_distrib(), log_transform()),
    inv_gamma = transformation(gamma_distrib(), inverse_transform()),
    sqrt_gamma = transformation(gamma_distrib(), sqrt_transform()),
    power2_gamma = transformation(gamma_distrib(), power_transform(2)),
    bc_gamma = transformation(gamma_distrib(), bc_transform(0.5)),
    softplus_gamma = transformation(gamma_distrib(), softplus_transform()),
    asinh_gauss = transformation(gaussian_distrib(), asinh_transform()),
    yj_gauss = transformation(gaussian_distrib(), yj_transform(0.7)),
    expit_gauss = transformation(gaussian_distrib(), expit_transform()),
    affine_gauss = transformation(gaussian_distrib(), affine_transform(1, -2)),
    logit_beta = transformation(beta_distrib(), logit_transform())
  )
  thetas <- list(
    log_gamma = thga, inv_gamma = thga, sqrt_gamma = thga, power2_gamma = thga,
    bc_gamma = thga, softplus_gamma = thga,
    asinh_gauss = thg, yj_gauss = thg, expit_gauss = thg, affine_gauss = thg,
    logit_beta = list(mu = 0.4, phi = 6)
  )

  for (nm in names(cases)) {
    d <- cases[[nm]]
    b <- d@bounds
    expect_equal(
      stats::integrate(function(t) distrib_pdf(d, t, thetas[[nm]]), b[1], b[2])$value,
      1,
      tolerance = 1e-5,
      label = paste(nm, "integral")
    )
  }
})

test_that("decreasing transformations swap tails coherently", {
  d <- transformation(gamma_distrib(), inverse_transform())
  th <- list(mu = 3, sigma2 = 2)
  p <- c(0.1, 0.5, 0.9)
  q <- distrib_quantile(d, p, th)

  expect_true(all(diff(q) > 0)) # quantiles increasing in p
  expect_equal(distrib_cdf(d, q, th), p, tolerance = 1e-9)

  # Affine with negative scale
  d2 <- transformation(gaussian_distrib(), affine_transform(loc = 1, scale = -2))
  thg <- list(mu = 0, sigma = 1)
  expect_equal(distrib_pdf(d2, 0.5, thg), dnorm(0.5, 1, 2))
  expect_equal(distrib_cdf(d2, 1, thg), 0.5, tolerance = 1e-12)
})

test_that("transformed derivatives delegate correctly (finite differences)", {
  set.seed(51)
  d <- transformation(gamma_distrib(), log_transform())
  th <- list(mu = 3, sigma2 = 2)
  y <- distrib_rng(d, 25, th)

  a_grad <- distrib_gradient(d, y, th)
  n_grad <- fd_gradient_ref(d, y, th)
  for (p in names(th)) {
    expect_equal(a_grad[[p]], n_grad[[p]], tolerance = 1e-5, label = paste("grad", p))
  }

  a_hess <- distrib_hessian(d, y, th)
  n_hess <- fd_hessian_ref(d, y, th)
  for (p in names(n_hess)) {
    expect_equal(a_hess[[p]], n_hess[[p]], tolerance = 1e-4, label = paste("hess", p))
  }

  # Expected hessian is exactly the parent's (Jacobian is parameter-free)
  eh_t <- distrib_expected_hessian(d, y, th)
  eh_p <- distrib_expected_hessian(gamma_distrib(), y, th)
  for (p in names(eh_p)) {
    expect_equal(eh_t[[p]], eh_p[[p]], label = paste("expected hess", p))
  }
})

test_that("transformed moments are computed numerically", {
  # E[log X], X ~ gamma(shape a, rate r) = digamma(a) - log(r)
  d <- transformation(gamma_distrib(), log_transform())
  th <- list(mu = 3, sigma2 = 2)
  a <- 9 / 2
  r <- 3 / 2
  expect_equal(mean(d, th), digamma(a) - log(r), tolerance = 1e-6)

  # E[a + bX] = a + b E[X]
  d2 <- transformation(gaussian_distrib(), affine_transform(1, 2))
  expect_equal(mean(d2, list(mu = 3, sigma = 1)), 1 + 2 * 3, tolerance = 1e-6)
  expect_equal(variance(d2, list(mu = 3, sigma = 1)), 4, tolerance = 1e-6)
})

test_that("invalid supports and inputs are rejected", {
  expect_error(transformation(gaussian_distrib(), log_transform()), "not valid")
  expect_error(transformation(gaussian_distrib(), sqrt_transform()), "not valid")
  expect_error(transformation(poisson_distrib(), log_transform()), "continuous")
  expect_error(transformation(gaussian_distrib(), "not a transformer"), "transformer")
  expect_error(affine_transform(scale = 0), "zero")
  expect_error(softplus_transform(a = -1), "greater than 0")

  # Box-Cox lambda = 0 degenerates to log
  bc0 <- bc_transform(0)
  expect_equal(bc0@name, "box_cox_0")
  expect_equal(bc0@trans_fun(exp(2)), 2)
})

test_that("transformations compose", {
  # log(exp(gaussian)) round-trips to the original gaussian
  d <- transformation(transformation(gaussian_distrib(), exp_transform()), log_transform())
  th <- list(mu = 0.4, sigma = 1.1)
  y <- c(-1, 0, 2)
  expect_equal(distrib_pdf(d, y, th), dnorm(y, 0.4, 1.1), tolerance = 1e-12)
})
