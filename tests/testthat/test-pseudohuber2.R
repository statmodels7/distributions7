# Pseudo-Huber in its standard-deviation parametrization: the density against
# the scale parametrization, the moments, and every derivative order against
# an independent reference.

# R(nu) = sqrt(nu) K2 / K1 with base R's Bessel functions
ph2_R <- function(nu) sqrt(nu) * besselK(sqrt(nu), 2) / besselK(sqrt(nu), 1)

# the log-density written from the formula, sharing no code with the package
ph2_ll <- function(p, y) {
  R <- ph2_R(p[3])
  D <- sqrt(p[3] + R * (y - p[1])^2 / p[2]^2)
  -D - log(2) - log(p[2]) - 0.25 * log(p[3]) +
    0.5 * log(besselK(sqrt(p[3]), 2)) - 1.5 * log(besselK(sqrt(p[3]), 1))
}

test_that("pseudohuber2 is pseudohuber at sigma1 = sigma / sqrt(R(nu))", {
  d <- pseudohuber2_distrib()
  d1 <- pseudohuber_distrib()
  y <- c(-4, -0.3, 0.4, 2.2, 9)
  for (th in list(list(mu = 0.4, sigma = 1.5, nu = 2.2),
                  list(mu = -1, sigma = 0.3, nu = 0.05),
                  list(mu = 2, sigma = 4, nu = 400))) {
    th1 <- list(mu = th$mu, sigma = th$sigma / sqrt(ph2_R(th$nu)), nu = th$nu)
    expect_equal(distrib_pdf(d, y, th, log = TRUE),
                 distrib_pdf(d1, y, th1, log = TRUE), tolerance = 1e-12)
    expect_equal(distrib_pdf(d, y, th, log = TRUE),
                 vapply(y, function(v) ph2_ll(unlist(th), v), 0),
                 tolerance = 1e-12)
  }
})

test_that("pseudohuber2 has variance sigma^2 and the family's kurtosis", {
  d <- pseudohuber2_distrib()
  d1 <- pseudohuber_distrib()
  for (th in list(list(mu = 0.4, sigma = 1.5, nu = 2.2),
                  list(mu = 0, sigma = 0.7, nu = 0.1))) {
    v <- stats::integrate(function(t) (t - th$mu)^2 * distrib_pdf(d, t, th),
                          -Inf, Inf, rel.tol = 1e-10)$value
    expect_equal(v, th$sigma^2, tolerance = 1e-8)
    expect_equal(variance(d, th), th$sigma^2)
    expect_equal(kurtosis(d, th), kurtosis(d1, th))
  }
  # the gaussian limit at fixed sigma
  yy <- c(-1, 0.5, 2)
  expect_equal(distrib_pdf(d, yy, list(mu = 0.4, sigma = 1.5, nu = 1e12)),
               stats::dnorm(yy, 0.4, 1.5), tolerance = 1e-6)
})

test_that("pseudohuber2 orders one and two match numDeriv on the formula", {
  d <- pseudohuber2_distrib()
  y <- c(-2.3, -0.4, 0.4, 1.1, 3.7, 8)
  p0 <- c(0.4, 1.5, 2.2)
  th <- list(mu = p0[1], sigma = p0[2], nu = p0[3])
  g <- distrib_gradient(d, y, th)
  G <- t(sapply(y, function(v) numDeriv::grad(function(p) ph2_ll(p, v), p0)))
  expect_equal(unname(do.call(cbind, g)), G, tolerance = 1e-8)
  H <- distrib_hessian(d, y, th)
  expect_named(H, c("mu_mu", "sigma_sigma", "nu_nu", "mu_sigma", "mu_nu",
                    "sigma_nu"))
  Hn <- t(sapply(y, function(v) {
    h <- numDeriv::hessian(function(p) ph2_ll(p, v), p0)
    c(h[1, 1], h[2, 2], h[3, 3], h[1, 2], h[1, 3], h[2, 3])
  }))
  expect_equal(unname(do.call(cbind, H)), Hn, tolerance = 1e-7)
})

test_that("pseudohuber2 orders three to five are the derivative of the order below", {
  d <- pseudohuber2_distrib()
  y <- c(-2.3, 0.4, 3.7)
  p0 <- c(0.4, 1.5, 2.2)
  nm <- c("mu", "sigma", "nu")
  lower <- list(`3` = distrib_hessian, `4` = function(...) distrib_deriv3(...),
                `5` = function(...) distrib_deriv4(...))
  upper <- list(`3` = function(...) distrib_deriv3(...),
                `4` = function(...) distrib_deriv4(...),
                `5` = function(...) distrib_deriv5(...))
  for (k in c("3", "4", "5")) {
    hi <- upper[[k]](d, y, list(mu = p0[1], sigma = p0[2], nu = p0[3]))
    expect_length(hi, c(`3` = 10, `4` = 15, `5` = 21)[[k]])
    for (key in names(hi)) {
      parts <- strsplit(key, "_")[[1]]
      last <- match(parts[length(parts)], nm)
      lo_key <- paste(parts[-length(parts)], collapse = "_")
      num <- sapply(seq_along(y), function(i) numDeriv::grad(function(p) {
        lo <- lower[[k]](d, y[i], list(mu = p[1], sigma = p[2], nu = p[3]))
        if (is.null(lo[[lo_key]])) lo_key <- paste(rev(strsplit(lo_key, "_")[[1]]),
                                                   collapse = "_")
        lo[[lo_key]]
      }, p0)[last])
      expect_equal(hi[[key]], num, tolerance = 1e-7, label = paste("order", k, key))
    }
  }
  expect_identical(names(distrib_deriv5(d, y, list(mu = 0.4, sigma = 1.5, nu = 2.2))),
                   names(numerical_deriv5(d, y, list(mu = 0.4, sigma = 1.5, nu = 2.2))))
})

test_that("pseudohuber2 response derivatives are minus the location ones", {
  d <- pseudohuber2_distrib()
  th <- list(mu = 0.4, sigma = 1.5, nu = 2.2)
  y <- c(-1, 0.5, 3)
  expect_equal(distrib_grad_y(d, y, th), -distrib_gradient(d, y, th)$mu)
  expect_equal(distrib_hess_y(d, y, th), distrib_hessian(d, y, th)$mu_mu)
  cy <- distrib_cross_y(d, y, th)
  expect_equal(cy$nu, -distrib_hessian(d, y, th)$mu_nu)
})

test_that("pseudohuber2 cdf, quantile and generator agree", {
  d <- pseudohuber2_distrib()
  th <- list(mu = 0.4, sigma = 1.5, nu = 2.2)
  q <- distrib_quantile(d, c(0.05, 0.3, 0.5, 0.9), th)
  expect_equal(distrib_cdf(d, q, th), c(0.05, 0.3, 0.5, 0.9), tolerance = 1e-7)
  expect_equal(distrib_cdf(d, 1.7, th),
               stats::integrate(function(t) distrib_pdf(d, t, th), -Inf, 1.7)$value,
               tolerance = 1e-8)
})

test_that("pseudohuber2 starts from its moments and below gaussian kurtosis at nu = 10", {
  d <- pseudohuber2_distrib()
  set.seed(5)
  y <- distrib_rng(d, 400, list(mu = 1, sigma = 2, nu = 1))
  st <- distrib_start(d, y, 1L)[[1]]
  expect_equal(st$sigma, sd(y))
  # the kurtosis of the start is the sample's, held in [0.85, 2.9]
  z <- (y - mean(y)) / sd(y)
  expect_equal(kurtosis(d, st), min(max(mean(z^4) - 3, 0.85), 2.9),
               tolerance = 1e-6)
  s <- distrib_intercept_start(d, runif(500, 2, 6))
  expect_equal(s$nu, 10)
})

test_that("the pseudo-Huber generator is the normal variance mixture", {
  # the ratio-of-uniforms box against a numerical extremum of (x - 1) sqrt(h)
  for (om in c(1e-3, 0.3, 1, 5, 100)) {
    f <- function(x) (x - 1) * exp(-om * (x - 1)^2 / (4 * x))
    b <- pseudohuber_rou_box_cpp(om)
    expect_equal(b[1], optimize(f, c(1e-12, 1))$objective, tolerance = 1e-7)
    expect_equal(b[2], optimize(f, c(1, 10 / om + 10), maximum = TRUE)$objective,
                 tolerance = 1e-7)
  }
  # the empirical cdf at the deciles, within five binomial standard errors
  d <- pseudohuber_distrib()
  p <- c(0.05, 0.25, 0.5, 0.75, 0.95)
  n <- 2e5
  for (nu in c(1e-6, 1, 1e4)) {
    th <- list(mu = 0.4, sigma = 1.3, nu = nu)
    set.seed(3)
    x <- distrib_rng(d, n, th)
    Fe <- vapply(distrib_quantile(d, p, th), function(q) mean(x <= q), 0)
    expect_lt(max(abs(Fe - p) / sqrt(p * (1 - p) / n)), 5)
  }
  # pseudohuber2 draws with variance sigma^2
  set.seed(4)
  x2 <- distrib_rng(pseudohuber2_distrib(), n, list(mu = -1, sigma = 2, nu = 0.7))
  expect_equal(var(x2), 4, tolerance = 0.02)
})
