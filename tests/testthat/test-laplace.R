# Laplace distribution and the handling of a parameter the log-likelihood is not
# differentiable in.

lap_cdf <- function(q, mu, b) ifelse(q < mu, 0.5 * exp((q - mu) / b), 1 - 0.5 * exp(-(q - mu) / b))

test_that("laplace pdf/cdf/quantile/rng/moments are correct", {
  d <- laplace_distrib()
  th <- list(mu = 1, sigma = 2)

  expect_equal(stats::integrate(function(t) distrib_pdf(d, t, th), -Inf, Inf)$value, 1, tolerance = 1e-6)
  expect_equal(distrib_pdf(d, 1, th), 1 / (2 * 2))
  q <- c(-3, -1, 0, 1, 4)
  expect_equal(distrib_cdf(d, q, th), lap_cdf(q, 1, 2))
  p <- c(0.05, 0.25, 0.5, 0.75, 0.95)
  expect_equal(distrib_quantile(d, 0.5, th), 1)
  expect_equal(distrib_cdf(d, distrib_quantile(d, p, th), th), p, tolerance = 1e-9)

  expect_equal(mean(d, th), 1)
  expect_equal(variance(d, th), 2 * 2^2)
  expect_equal(skewness(d, th), 0)
  expect_equal(kurtosis(d, th), 3)

  set.seed(1)
  y <- distrib_rng(d, 1e5, th)
  expect_equal(mean(y), 1, tolerance = 0.05)
  expect_equal(stats::var(y), 8, tolerance = 0.3)
})

test_that("laplace score matches finite differences away from the kink", {
  d <- laplace_distrib()
  th <- list(mu = 1, sigma = 2)
  y <- c(-2, 0.3, 3.7) # none equal to mu
  a <- distrib_gradient(d, y, th)
  fd <- function(j, h = 1e-6) {
    tp <- tm <- th; tp[[j]] <- th[[j]] + h; tm[[j]] <- th[[j]] - h
    (distrib_pdf(d, y, tp, log = TRUE) - distrib_pdf(d, y, tm, log = TRUE)) / (2 * h)
  }
  expect_equal(a$mu, fd(1), tolerance = 1e-5)
  expect_equal(a$sigma, fd(2), tolerance = 1e-5)
})

test_that("laplace observed Hessian is degenerate in mu but the expected Hessian is the Fisher information", {
  d <- laplace_distrib()
  th <- list(mu = 1, sigma = 2)
  y <- c(-2, 0.3, 3.7)

  h <- distrib_hessian(d, y, th)
  expect_equal(h$mu_mu, rep(0, 3))          # kink -> zero observed curvature in mu
  expect_equal(h$sigma_sigma, (2 - 2 * abs(y - 1)) / 2^3)
  expect_equal(h$mu_sigma, -sign(y - 1) / 2^2)

  # Expected Hessian = -Fisher information = c(mu_mu = -1/s^2, sigma_sigma = -1/s^2, mu_sigma = 0)
  eh <- distrib_expected_hessian(d, c(0, 1, 2), th)
  expect_equal(eh$mu_mu, rep(-1 / 4, 3))
  expect_equal(eh$sigma_sigma, rep(-1 / 4, 3))
  expect_equal(eh$mu_sigma, rep(0, 3))

  # confirm it disagrees with the (degenerate) expectation of the observed Hessian
  set.seed(2)
  ys <- distrib_rng(d, 2e5, th)
  expect_equal(mean(distrib_hessian(d, ys, th)$mu_mu), 0) # E[H mu_mu] = 0, NOT the Fisher info
})

test_that("params_smooth metadata is exposed and validated", {
  d <- laplace_distrib()
  expect_equal(param_smoothness(d), c(mu = FALSE, sigma = TRUE))
  # defaults to all TRUE when unset
  expect_true(all(param_smoothness(gaussian1_distrib())))
  expect_named(param_smoothness(gaussian1_distrib()), c("mu", "sigma"))

  # print marks the non-smooth parameter
  out <- paste(utils::capture.output(print(d)), collapse = "\n")
  expect_match(out, "non-smooth", fixed = TRUE)
})

test_that("OPG fallback yields the correct Fisher information for a bare distribution", {
  # a user-defined gaussian with only the density: expected Hessian via OPG must
  # match the analytical -I of the built-in gaussian
  BareGaussOPG <- S7::new_class("BareGaussOPG", parent = continuous_distrib, package = NULL)
  S7::method(distrib_pdf, BareGaussOPG) <- function(distrib, y, theta, log = FALSE) {
    stats::dnorm(y, theta[[1]], theta[[2]], log = log)
  }
  bg <- BareGaussOPG(
    distrib_name = "bareg", dimension = "univariate", bounds = c(-Inf, Inf),
    params = c("mu", "sigma"), params_interpretation = c(mu = "m", sigma = "s"), n_params = 2,
    params_bounds = list(mu = c(-Inf, Inf), sigma = c(0, Inf)),
    link_params = list(mu = linkfunctions7::identity_link(), sigma = linkfunctions7::log_link())
  )
  th <- list(mu = 1.5, sigma = 2)
  eh_opg <- distrib_expected_hessian(bg, 0, th)
  eh_ana <- distrib_expected_hessian(gaussian1_distrib(), 0, th)
  for (nm in names(eh_ana)) {
    expect_equal(eh_opg[[nm]][1], eh_ana[[nm]][1], tolerance = 1e-5, label = nm)
  }
})

test_that("response derivatives (grad_y / hess_y) work analytically and via fallback", {
  y <- c(-2, 0.3, 3.7)

  # gaussian: closed form
  g <- gaussian1_distrib(); thg <- list(mu = 0.5, sigma = 2)
  expect_equal(distrib_grad_y(g, y, thg), -(y - 0.5) / 4)
  expect_equal(distrib_hess_y(g, y, thg), rep(-1 / 4, 3))

  # gamma: numerical fallback vs the analytical response derivatives
  gam <- gamma2_distrib(); thga <- list(mu = 3, sigma2 = 2)
  a <- 9 / 2; rate <- 3 / 2; yg <- c(1, 3, 6)
  expect_equal(distrib_grad_y(gam, yg, thga), (a - 1) / yg - rate, tolerance = 1e-5)
  expect_equal(distrib_hess_y(gam, yg, thga), -(a - 1) / yg^2, tolerance = 1e-5)

  # laplace: closed form (kink-aware)
  d <- laplace_distrib(); th <- list(mu = 1, sigma = 2)
  expect_equal(distrib_grad_y(d, y, th), -sign(y - 1) / 2)
  expect_equal(distrib_hess_y(d, y, th), rep(0, 3))
})
