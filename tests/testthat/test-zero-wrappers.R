# Zero-inflated and zero-adjusted wrappers: probability functions, moments,
# and derivatives (finite differences + expectation-based expected Hessian).

zi_za_cases <- function() {
  list(
    zip = list(d = zero_inflated(poisson_distrib()), theta = list(mu = 3, zi = 0.25)),
    zinb = list(d = zero_inflated(negbin_distrib()), theta = list(mu = 4, theta = 1.5, zi = 0.3)),
    zap = list(d = zero_adjusted(poisson_distrib()), theta = list(mu = 3, za = 0.4)),
    zagamma = list(d = zero_adjusted(gamma_distrib()), theta = list(mu = 3, sigma2 = 2, za = 0.3))
  )
}

test_that("wrapped pmfs/pdfs are normalized and place the right mass at zero", {
  # ZIP: P(0) = zi + (1-zi) dpois(0, mu)
  zip <- zero_inflated(poisson_distrib())
  th <- list(mu = 3, zi = 0.25)
  expect_equal(distrib_pdf(zip, 0, th), 0.25 + 0.75 * dpois(0, 3))
  expect_equal(distrib_pdf(zip, 2, th), 0.75 * dpois(2, 3))
  expect_equal(numerical_series(function(k) distrib_pdf(zip, k, th), 0, Inf), 1, tolerance = 1e-9)

  # ZAP (hurdle): P(0) = za exactly, positives renormalized
  zap <- zero_adjusted(poisson_distrib())
  th <- list(mu = 3, za = 0.4)
  expect_equal(distrib_pdf(zap, 0, th), 0.4)
  expect_equal(distrib_pdf(zap, 2, th), 0.6 * dpois(2, 3) / (1 - dpois(0, 3)))
  expect_equal(numerical_series(function(k) distrib_pdf(zap, k, th), 0, Inf), 1, tolerance = 1e-9)

  # ZA gamma: P(0) = za, continuous part scaled by (1-za)
  zag <- zero_adjusted(gamma_distrib())
  th <- list(mu = 3, sigma2 = 2, za = 0.3)
  expect_equal(distrib_pdf(zag, 0, th), 0.3)
  expect_equal(
    stats::integrate(function(t) distrib_pdf(zag, t, th), 1e-12, Inf)$value,
    0.7,
    tolerance = 1e-6
  )
})

test_that("wrapped cdf/quantile are mutually consistent", {
  for (nm in names(zi_za_cases())) {
    case <- zi_za_cases()[[nm]]
    d <- case$d
    th <- case$theta
    p <- c(0.05, 0.3, 0.7, 0.95)
    q <- distrib_quantile(d, p, th)

    if (S7::S7_inherits(d, discrete_distrib)) {
      expect_true(all(distrib_cdf(d, q, th) >= p - 1e-12), label = paste(nm, "F(q) >= p"))
      expect_true(all(distrib_cdf(d, q - 1, th) < p), label = paste(nm, "F(q-1) < p"))
    } else {
      pos <- q > 0
      expect_equal(distrib_cdf(d, q[pos], th), p[pos], tolerance = 1e-7, label = nm)
    }
  }
})

test_that("wrapped moments match closed forms", {
  # ZIP
  zip <- zero_inflated(poisson_distrib())
  th <- list(mu = 3, zi = 0.25)
  expect_equal(mean(zip, th), 0.75 * 3, tolerance = 1e-7)
  expect_equal(variance(zip, th), 0.75 * 3 + 0.25 * 0.75 * 9, tolerance = 1e-6)

  # ZAP (hurdle): E[Y] = (1-za) mu / (1 - exp(-mu))
  zap <- zero_adjusted(poisson_distrib())
  th <- list(mu = 3, za = 0.4)
  expect_equal(mean(zap, th), 0.6 * 3 / (1 - dpois(0, 3)), tolerance = 1e-7)

  # ZA gamma (GAMLSS 9.2): E[Y] = (1-p) mu, V = (1-p) s2 + p(1-p) mu^2
  zag <- zero_adjusted(gamma_distrib())
  th <- list(mu = 3, sigma2 = 2, za = 0.3)
  expect_equal(mean(zag, th), 0.7 * 3, tolerance = 1e-6)
  expect_equal(variance(zag, th), 0.7 * 2 + 0.3 * 0.7 * 9, tolerance = 1e-5)
})

test_that("wrapped gradients and hessians match finite differences", {
  set.seed(41)
  for (nm in names(zi_za_cases())) {
    case <- zi_za_cases()[[nm]]
    d <- case$d
    th <- case$theta
    y <- distrib_rng(d, 25, th)
    stopifnot(any(y == 0)) # make sure both branches are exercised

    a_grad <- distrib_gradient(d, y, th)
    n_grad <- fd_gradient_ref(d, y, th)
    for (p in names(th)) {
      expect_equal(a_grad[[p]], n_grad[[p]], tolerance = 1e-5,
        label = sprintf("%s grad %s", nm, p))
    }

    a_hess <- distrib_hessian(d, y, th)
    n_hess <- fd_hessian_ref(d, y, th)
    expect_setequal(names(a_hess), names(n_hess))
    for (p in names(n_hess)) {
      expect_equal(a_hess[[p]], n_hess[[p]], tolerance = 1e-4,
        label = sprintf("%s hess %s", nm, p))
    }
  }
})

test_that("wrapped expected hessians equal the expectation of the observed hessians", {
  for (nm in names(zi_za_cases())) {
    case <- zi_za_cases()[[nm]]
    d <- case$d
    th <- case$theta

    a_exp <- distrib_expected_hessian(d, 0, th)
    for (p in names(a_exp)) {
      e_num <- expectation(
        d,
        function(y, theta) distrib_hessian(d, y, theta)[[p]],
        th
      )
      expect_equal(a_exp[[p]][1], e_num, tolerance = 1e-5,
        label = sprintf("%s expected hess %s", nm, p),
        expected.label = "numerical expectation")
    }
  }
})

test_that("wrapper constructors validate their input", {
  expect_error(zero_inflated(gaussian_distrib()), "discrete")
  expect_error(zero_inflated(zero_inflated(poisson_distrib())), "zi")
  expect_error(zero_adjusted(zero_adjusted(poisson_distrib())), "za")
})
