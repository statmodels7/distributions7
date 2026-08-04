# Third/fourth-order derivatives: closed-form C++ kernels validated against the
# finite-difference fallback (observed) and against Monte Carlo (expected), plus
# the numerical fallback itself on distributions without analytic kernels.

test_that("deriv_names enumerates unique multi-indices, diagonal-first", {
  expect_equal(deriv_names("mu", 3), "mu_mu_mu")
  expect_equal(
    deriv_names(c("mu", "sigma"), 3),
    c("mu_mu_mu", "mu_mu_sigma", "mu_sigma_sigma", "sigma_sigma_sigma")
  )
  # counts = combinations with repetition
  expect_length(deriv_names(c("a", "b", "c"), 3), 10)
  expect_length(deriv_names(c("a", "b", "c"), 4), 15)
})

# distributions with a closed-form analytic higher-derivative kernel
analytic_hd_cases <- function() {
  list(
    gaussian  = list(d = gaussian_distrib(),  theta = list(mu = 1.5, sigma = 2)),
    cauchy    = list(d = cauchy_distrib(),    theta = list(mu = 0.5, sigma = 1.4)),
    invgauss  = list(d = invgauss_distrib(),  theta = list(mu = 2, phi = 0.7)),
    bernoulli = list(d = bernoulli_distrib(), theta = list(mu = 0.35)),
    binomial  = list(d = binomial_distrib(size = 10), theta = list(mu = 0.35)),
    gamma     = list(d = gamma_distrib(),     theta = list(mu = 3, sigma2 = 2)),
    beta      = list(d = beta_distrib(),      theta = list(mu = 0.4, phi = 6)),
    logistic  = list(d = logistic_distrib(),  theta = list(mu = 0.5, sigma = 1.4)),
    student_t = list(d = student_t_distrib(), theta = list(mu = 0.5, sigma = 1.3, nu = 6)),
    lognormal = list(d = lognormal_distrib(), theta = list(mu = 0.5, sigma2 = 1.3)),
    poisson   = list(d = poisson_distrib(),   theta = list(mu = 4)),
    negbin    = list(d = negbin_distrib(),    theta = list(mu = 4, theta = 1.7)),
    pseudohuber = list(d = pseudohuber_distrib(), theta = list(mu = 0.5, sigma = 1.4, nu = 2.5)),
    weibull    = list(d = weibull_distrib(),    theta = list(mu = 2, sigma = 1.6)),
    gumbel     = list(d = gumbel_distrib(),     theta = list(mu = 1, sigma = 2.5)),
    laplace    = list(d = laplace_distrib(),    theta = list(mu = 0.5, b = 1.3)),
    skewnormal = list(d = skewnormal_distrib(), theta = list(mu = 0.3, sigma = 1.4, alpha = 3))
  )
}

test_that("analytic observed 3rd/4th derivatives match the numerical fallback", {
  set.seed(1)
  for (nm in names(analytic_hd_cases())) {
    case <- analytic_hd_cases()[[nm]]
    d <- case$d
    th <- case$theta
    y <- distrib_rng(d, 40, th)

    a3 <- distrib_deriv3(d, y, th)
    n3 <- numerical_deriv3(d, y, th)
    expect_setequal(names(a3), deriv_names(d@params, 3))
    for (k in names(a3)) {
      expect_equal(a3[[k]], n3[[k]], tolerance = 1e-4,
        label = sprintf("%s d3 %s", nm, k), expected.label = "finite differences")
    }

    a4 <- distrib_deriv4(d, y, th)
    n4 <- numerical_deriv4(d, y, th)
    expect_setequal(names(a4), deriv_names(d@params, 4))
    for (k in names(a4)) {
      expect_equal(a4[[k]], n4[[k]], tolerance = 1e-3,
        label = sprintf("%s d4 %s", nm, k), expected.label = "finite differences")
    }
  }
})

test_that("closed-form expected derivatives match a Monte Carlo mean of the observed", {
  # only the distributions whose expected higher derivatives are closed-form
  cases <- analytic_hd_cases()
  # these have closed-form observed but no closed-form expected derivatives
  cases$student_t <- NULL
  cases$pseudohuber <- NULL
  cases$logistic <- NULL
  cases$skewnormal <- NULL
  set.seed(2)

  for (nm in names(cases)) {
    d <- cases[[nm]]$d
    th <- cases[[nm]]$theta
    ys <- distrib_rng(d, 5e5, th)

    for (ord in 3:4) {
      ea <- if (ord == 3) distrib_deriv3(d, 0, th, expected = TRUE) else distrib_deriv4(d, 0, th, expected = TRUE)
      oa <- if (ord == 3) distrib_deriv3(d, ys, th) else distrib_deriv4(d, ys, th)
      for (k in names(ea)) {
        mc <- mean(oa[[k]]); se <- stats::sd(oa[[k]]) / sqrt(length(ys))
        # allow a generous |z| (finite MC) plus an absolute floor
        expect_lt(abs(ea[[k]][1] - mc), 8 * se + 1e-6)
      }
    }
  }
})

test_that("expected higher derivatives are returned per observation", {
  g <- gaussian_distrib()
  e3 <- distrib_deriv3(g, c(0, 1, 2), list(mu = 0, sigma = 1), expected = TRUE)
  expect_length(e3$sigma_sigma_sigma, 3)
  expect_true(all(e3$sigma_sigma_sigma == e3$sigma_sigma_sigma[1]))
})

test_that("student_t expected derivatives (no closed form) fall back numerically", {
  st <- student_t_distrib()
  th <- list(mu = 0, sigma = 1, nu = 8)
  e3 <- distrib_deriv3(st, 0, th, expected = TRUE)
  expect_setequal(names(e3), deriv_names(st@params, 3))
  expect_true(all(is.finite(unlist(e3))))
  # symmetry zeros
  expect_equal(e3$mu_mu_mu, 0, tolerance = 1e-3)
  expect_equal(e3$mu_sigma_sigma, 0, tolerance = 1e-3)
})

test_that("closed-form response derivatives match the numerical fallback", {
  set.seed(7)
  cont <- list(
    gaussian    = list(d = gaussian_distrib(),    th = list(mu = 1.5, sigma = 2)),
    cauchy      = list(d = cauchy_distrib(),      th = list(mu = 0.5, sigma = 1.4)),
    logistic    = list(d = logistic_distrib(),    th = list(mu = 0.5, sigma = 1.4)),
    student_t   = list(d = student_t_distrib(),   th = list(mu = 0.5, sigma = 1.3, nu = 6)),
    pseudohuber = list(d = pseudohuber_distrib(), th = list(mu = 0.5, sigma = 1.4, nu = 2.5)),
    gamma       = list(d = gamma_distrib(),       th = list(mu = 3, sigma2 = 2)),
    invgauss    = list(d = invgauss_distrib(),    th = list(mu = 2, phi = 0.7)),
    lognormal   = list(d = lognormal_distrib(),   th = list(mu = 0.5, sigma2 = 1.3)),
    beta        = list(d = beta_distrib(),        th = list(mu = 0.4, phi = 6))
  )
  for (nm in names(cont)) {
    d <- cont[[nm]]$d; th <- cont[[nm]]$th
    y <- distrib_rng(d, 40, th)
    expect_equal(distrib_grad_y(d, y, th), numerical_grad_y(d, y, th),
      tolerance = 1e-4, label = paste(nm, "grad_y"))
    expect_equal(distrib_hess_y(d, y, th), numerical_hess_y(d, y, th),
      tolerance = 1e-3, label = paste(nm, "hess_y"))
  }
})

test_that("logistic higher derivatives stay finite at extreme residuals", {
  # The closed form is written in the sigmoid rather than in exp((mu-y)/sigma),
  # which the Wolfram output uses up to the fourth power and which overflows.
  lg <- logistic_distrib()
  th <- list(mu = 0, sigma = 1)
  y <- c(-4000, -800, 0, 800, 4000)
  d3 <- distrib_deriv3(lg, y, th)
  d4 <- distrib_deriv4(lg, y, th)
  expect_true(all(is.finite(unlist(d3))))
  expect_true(all(is.finite(unlist(d4))))
  # For |z| large, g1 -> -sign(z) and the other g's vanish exponentially
  expect_equal(d4$sigma_sigma_sigma_sigma, 6 - 24 * abs(y), tolerance = 1e-8)
})

test_that("the logistic expected higher derivatives match their exact values", {
  # Wolfram returns complex-valued, polylogarithmic and spuriously mu-dependent
  # expressions for these, so they are obtained through `expected_derivative()`.
  # These seven components are known exactly; the location-scale structure
  # forbids any dependence on mu.
  s <- 1.4
  th <- list(mu = 0.5, sigma = s)
  e3 <- distrib_deriv3(logistic_distrib(), 0, th, expected = TRUE)
  e4 <- distrib_deriv4(logistic_distrib(), 0, th, expected = TRUE)

  expect_equal(e3$mu_mu_mu[1], 0, tolerance = 1e-6)
  expect_equal(e3$mu_mu_sigma[1], 1 / (2 * s^3), tolerance = 1e-6)
  expect_equal(e3$mu_sigma_sigma[1], 0, tolerance = 1e-6)
  expect_equal(e3$sigma_sigma_sigma[1], (pi^2 + 2) / (2 * s^3), tolerance = 1e-6)
  expect_equal(e4$mu_mu_mu_mu[1], 1 / (15 * s^4), tolerance = 1e-6)
  expect_equal(e4$mu_mu_mu_sigma[1], 0, tolerance = 1e-6)
  expect_equal(e4$mu_sigma_sigma_sigma[1], 0, tolerance = 1e-6)
})

test_that("distributions without an analytic kernel still deliver 3rd/4th derivatives", {
  # a bare user-defined distribution: only a density, everything else numerical
  Bare <- S7::new_class("BareHD", parent = continuous_distrib, package = NULL)
  S7::method(distrib_pdf, Bare) <- function(distrib, y, theta, log = FALSE) {
    stats::dnorm(y, theta[[1]], theta[[2]], log = log)
  }
  bd <- Bare(
    distrib_name = "bare hd", dimension = "univariate", bounds = c(-Inf, Inf),
    params = c("mu", "sigma"), params_interpretation = c(mu = "loc", sigma = "scale"),
    n_params = 2, params_bounds = list(mu = c(-Inf, Inf), sigma = c(0, Inf)),
    link_params = list(
      mu = linkfunctions7::identity_link(),
      sigma = linkfunctions7::log_link()
    )
  )
  th <- list(mu = 0.5, sigma = 1.4)
  set.seed(3)
  y <- distrib_rng(bd, 20, th)
  d3 <- distrib_deriv3(bd, y, th)
  expect_setequal(names(d3), deriv_names(bd@params, 3))
  expect_true(all(is.finite(unlist(d3))))

  # a wrapper: zero-inflated poisson
  zip <- zero_inflated(poisson_distrib())
  d4 <- distrib_deriv4(zip, c(0, 1, 2, 3), list(mu = 3, zi = 0.2))
  expect_setequal(names(d4), deriv_names(zip@params, 4))
  expect_true(all(is.finite(unlist(d4))))
})
