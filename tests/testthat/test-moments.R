# Numerical moment machinery and analytical overrides.

test_that("numerical moments match known closed forms", {
  g <- gaussian_distrib()
  thg <- list(mu = 2, sigma = 3)
  expect_equal(mean(g, thg), 2, tolerance = 1e-6)
  expect_equal(variance(g, thg), 9, tolerance = 1e-6)
  expect_equal(skewness(g, thg), 0, tolerance = 1e-6)
  expect_equal(kurtosis(g, thg), 0, tolerance = 1e-5)

  # Gamma: shape a = mu^2/sigma2, skew = 2/sqrt(a), excess kurtosis = 6/a
  gam <- gamma_distrib()
  thga <- list(mu = 3, sigma2 = 2)
  a <- 9 / 2
  expect_equal(skewness(gam, thga), 2 / sqrt(a), tolerance = 1e-6)
  expect_equal(kurtosis(gam, thga), 6 / a, tolerance = 1e-5)

  # Poisson (discrete): all standardized cumulants known
  pois <- poisson_distrib()
  expect_equal(mean(pois, list(mu = 4)), 4, tolerance = 1e-8)
  expect_equal(variance(pois, list(mu = 4)), 4, tolerance = 1e-8)
  expect_equal(skewness(pois, list(mu = 4)), 1 / 2, tolerance = 1e-6)
  expect_equal(kurtosis(pois, list(mu = 4)), 1 / 4, tolerance = 1e-6)
})

test_that("moment() supports raw/central moments and vectorized theta", {
  g <- gaussian_distrib()
  expect_equal(moment(g, list(mu = 2, sigma = 3), p = 2), 4 + 9, tolerance = 1e-6)
  expect_equal(
    moment(g, list(mu = c(0, 1), sigma = c(1, 2)), p = 2, central = TRUE),
    c(1, 4),
    tolerance = 1e-6
  )
  expect_equal(
    variance(g, list(mu = c(0, 1), sigma = c(1, 2))),
    c(1, 4),
    tolerance = 1e-6
  )
})

test_that("negbin analytical moments agree with the numerical machinery", {
  d <- negbin_distrib()
  th <- list(mu = 4, theta = 1.7)

  expect_equal(mean(d, th), 4)
  expect_equal(variance(d, th), 4 + 16 / 1.7)

  m <- moment(d, th, p = 1)
  m2 <- moment(d, th, p = 2, central = TRUE, mu = m)
  m3 <- moment(d, th, p = 3, central = TRUE, mu = m)
  m4 <- moment(d, th, p = 4, central = TRUE, mu = m)

  expect_equal(variance(d, th), m2, tolerance = 1e-7)
  expect_equal(skewness(d, th), m3 / m2^1.5, tolerance = 1e-7)
  expect_equal(kurtosis(d, th), m4 / m2^2 - 3, tolerance = 1e-6)
})

test_that("pseudohuber analytical moments agree with the numerical machinery", {
  d <- pseudohuber_distrib()
  th <- list(mu = 0.5, sigma = 1.4, nu = 2.5)

  expect_equal(mean(d, th), 0.5, tolerance = 1e-7)
  expect_equal(variance(d, th), moment(d, th, p = 2, central = TRUE, mu = 0.5), tolerance = 1e-6)
  expect_equal(skewness(d, th), 0)
  m2 <- moment(d, th, p = 2, central = TRUE, mu = 0.5)
  m4 <- moment(d, th, p = 4, central = TRUE, mu = 0.5)
  expect_equal(kurtosis(d, th), m4 / m2^2 - 3, tolerance = 1e-5)
})

test_that("sample statistics methods work on numeric vectors", {
  x <- c(1, 2, 4, 8, 16)
  expect_equal(variance(x), stats::var(x))
  expect_equal(std_dev(x), stats::sd(x))

  m <- base::mean(x)
  s <- sqrt(base::mean((x - m)^2))
  expect_equal(skewness(x), base::mean((x - m)^3) / s^3)
  expect_equal(kurtosis(x), base::mean((x - m)^4) / s^4 - 3)
})

test_that("hess_names produces diagonal-first ordering", {
  expect_equal(hess_names("mu"), "mu_mu")
  expect_equal(hess_names(c("mu", "sigma")), c("mu_mu", "sigma_sigma", "mu_sigma"))
  expect_equal(
    hess_names(c("a", "b", "c")),
    c("a_a", "b_b", "c_c", "a_b", "a_c", "b_c")
  )
})
