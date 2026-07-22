# Probability-function fallbacks: a distribution that implements ONLY
# distrib_pdf must still provide cdf, quantile and rng through the default
# methods registered on continuous_distrib / discrete_distrib.

# --- bare continuous distribution: only the density ---
BareGaussF <- S7::new_class("BareGaussF", parent = continuous_distrib, package = NULL)
S7::method(distrib_pdf, BareGaussF) <- function(distrib, y, theta, log = FALSE) {
  stats::dnorm(y, theta[[1]], theta[[2]], log = log)
}
bare_gauss_f <- function() {
  BareGaussF(
    distrib_name = "bare gaussian", dimension = "univariate", bounds = c(-Inf, Inf),
    params = c("mu", "sigma"), params_interpretation = c(mu = "mean", sigma = "sd"),
    n_params = 2, params_bounds = list(mu = c(-Inf, Inf), sigma = c(0, Inf)),
    link_params = list(mu = linkfunctions7::identity_link(), sigma = linkfunctions7::log_link())
  )
}

# --- bare discrete distribution: only the pmf ---
BarePoisF <- S7::new_class("BarePoisF", parent = discrete_distrib, package = NULL)
S7::method(distrib_pdf, BarePoisF) <- function(distrib, y, theta, log = FALSE) {
  stats::dpois(y, theta[[1]], log = log)
}
bare_pois_f <- function() {
  BarePoisF(
    distrib_name = "bare poisson", dimension = "univariate", bounds = c(0, Inf),
    params = "mu", params_interpretation = c(mu = "mean"),
    n_params = 1, params_bounds = list(mu = c(0, Inf)),
    link_params = list(mu = linkfunctions7::log_link())
  )
}

test_that("continuous CDF fallback matches the analytical CDF", {
  bg <- bare_gauss_f()
  th <- list(mu = 1.5, sigma = 2)
  q <- c(-4, -2, 0, 1.5, 4, 8)
  expect_equal(distrib_cdf(bg, q, th), pnorm(q, 1.5, 2), tolerance = 1e-7)
  expect_equal(distrib_cdf(bg, q, th, lower.tail = FALSE), pnorm(q, 1.5, 2, lower.tail = FALSE), tolerance = 1e-7)
  expect_equal(distrib_cdf(bg, 1.5, th, log.p = TRUE), log(0.5), tolerance = 1e-7)
})

test_that("continuous quantile fallback matches the analytical quantile", {
  bg <- bare_gauss_f()
  th <- list(mu = 1.5, sigma = 2)
  p <- c(0.001, 0.05, 0.25, 0.5, 0.75, 0.95, 0.999)
  expect_equal(distrib_quantile(bg, p, th), qnorm(p, 1.5, 2), tolerance = 1e-6)
  # roundtrip
  expect_equal(distrib_cdf(bg, distrib_quantile(bg, p, th), th), p, tolerance = 1e-7)
  # tails / log.p / lower.tail
  expect_equal(distrib_quantile(bg, log(0.9), th, log.p = TRUE), qnorm(0.9, 1.5, 2), tolerance = 1e-6)
  expect_equal(distrib_quantile(bg, 0.9, th, lower.tail = FALSE), qnorm(0.9, 1.5, 2, lower.tail = FALSE), tolerance = 1e-6)
})

test_that("continuous RNG fallback produces the right distribution", {
  set.seed(101)
  bg <- bare_gauss_f()
  th <- list(mu = 1.5, sigma = 2)
  y <- distrib_rng(bg, 20000, th)
  expect_equal(mean(y), 1.5, tolerance = 0.05)
  expect_equal(sd(y), 2, tolerance = 0.05)
})

test_that("continuous fallbacks handle vectorized theta and bounded supports", {
  bg <- bare_gauss_f()
  expect_equal(
    distrib_cdf(bg, c(0, 1), list(mu = c(0, 1), sigma = c(1, 2))),
    c(0.5, 0.5), tolerance = 1e-7
  )

  # A bounded-support density (Beta): quantile fallback must stay inside (0,1)
  BareBeta <- S7::new_class("BareBeta", parent = continuous_distrib, package = NULL)
  S7::method(distrib_pdf, BareBeta) <- function(distrib, y, theta, log = FALSE) {
    stats::dbeta(y, theta[[1]], theta[[2]], log = log)
  }
  bb <- BareBeta(
    distrib_name = "bare beta", dimension = "univariate", bounds = c(0, 1),
    params = c("a", "b"), params_interpretation = c(a = "shape1", b = "shape2"),
    n_params = 2, params_bounds = list(a = c(0, Inf), b = c(0, Inf)),
    link_params = list(a = linkfunctions7::log_link(), b = linkfunctions7::log_link())
  )
  th <- list(a = 2, b = 5)
  p <- c(0.1, 0.5, 0.9)
  expect_equal(distrib_quantile(bb, p, th), qbeta(p, 2, 5), tolerance = 1e-6)
  expect_equal(distrib_cdf(bb, c(0.2, 0.5), th), pbeta(c(0.2, 0.5), 2, 5), tolerance = 1e-7)
})

test_that("discrete CDF fallback matches the analytical CDF", {
  bp <- bare_pois_f()
  th <- list(mu = 4)
  q <- c(0, 1, 2, 4, 8, 15)
  expect_equal(distrib_cdf(bp, q, th), ppois(q, 4), tolerance = 1e-12)
  expect_equal(distrib_cdf(bp, 2.7, th), ppois(2, 4), tolerance = 1e-12) # floors q
  expect_equal(distrib_cdf(bp, q, th, lower.tail = FALSE), ppois(q, 4, lower.tail = FALSE), tolerance = 1e-12)
})

test_that("discrete quantile fallback matches the analytical quantile", {
  bp <- bare_pois_f()
  th <- list(mu = 4)
  p <- c(0.01, 0.1, 0.25, 0.5, 0.75, 0.9, 0.99)
  expect_equal(distrib_quantile(bp, p, th), qpois(p, 4))
  # large mean forces the table to grow
  expect_equal(distrib_quantile(bp, 0.5, list(mu = 5000)), qpois(0.5, 5000))
  # vectorized theta
  expect_equal(distrib_quantile(bp, c(0.5, 0.5), list(mu = c(4, 100))), qpois(c(0.5, 0.5), c(4, 100)))
})

test_that("discrete RNG fallback produces the right distribution", {
  set.seed(102)
  bp <- bare_pois_f()
  y <- distrib_rng(bp, 20000, list(mu = 4))
  expect_equal(mean(y), 4, tolerance = 0.05)
  expect_equal(var(y), 4, tolerance = 0.1)
})

test_that("a genuinely new distribution (Laplace) works from the pdf alone", {
  Laplace <- S7::new_class("Laplace", parent = continuous_distrib, package = NULL)
  S7::method(distrib_pdf, Laplace) <- function(distrib, y, theta, log = FALSE) {
    ld <- -log(2 * theta[[2]]) - abs(y - theta[[1]]) / theta[[2]]
    if (log) ld else exp(ld)
  }
  lp <- Laplace(
    distrib_name = "laplace", dimension = "univariate", bounds = c(-Inf, Inf),
    params = c("mu", "b"), params_interpretation = c(mu = "location", b = "scale"),
    n_params = 2, params_bounds = list(mu = c(-Inf, Inf), b = c(0, Inf)),
    link_params = list(mu = linkfunctions7::identity_link(), b = linkfunctions7::log_link())
  )
  th <- list(mu = 1, b = 2)
  lap_cdf <- function(q, m, b) ifelse(q < m, 0.5 * exp((q - m) / b), 1 - 0.5 * exp(-(q - m) / b))

  expect_equal(distrib_quantile(lp, 0.5, th), 1, tolerance = 1e-6)   # median = mu
  expect_equal(variance(lp, th), 8, tolerance = 1e-5)                # 2 b^2
  q <- c(-3, 0, 1, 4)
  expect_equal(distrib_cdf(lp, q, th), lap_cdf(q, 1, 2), tolerance = 1e-7)
  expect_equal(distrib_cdf(lp, distrib_quantile(lp, c(0.2, 0.8), th), th), c(0.2, 0.8), tolerance = 1e-6)
})

test_that("discrete fallback rejects an infinite lower bound", {
  BareInt <- S7::new_class("BareInt", parent = discrete_distrib, package = NULL)
  S7::method(distrib_pdf, BareInt) <- function(distrib, y, theta, log = FALSE) {
    d <- stats::dnorm(y, 0, 5); if (log) log(d) else d
  }
  bi <- BareInt(
    distrib_name = "bare integer", dimension = "univariate", bounds = c(-Inf, Inf),
    params = "s", params_interpretation = c(s = "scale"), n_params = 1,
    params_bounds = list(s = c(0, Inf)), link_params = list(s = linkfunctions7::log_link())
  )
  expect_error(distrib_cdf(bi, 0, list(s = 5)), "finite lower bound")
})
