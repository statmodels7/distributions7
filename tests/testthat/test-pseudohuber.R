# Pseudo-Huber distribution: normalization, numerical CDF/quantile/RNG.
# (Gradient/Hessian correctness is covered by test-derivatives via finite
# differences once the density is validated here.)

test_that("pseudohuber density is correctly normalized", {
  d <- pseudohuber_distrib()
  for (th in list(
    list(mu = 0.5, sigma = 1.4, nu = 2.5),
    list(mu = -3, sigma = 0.5, nu = 0.4),
    list(mu = 0, sigma = 2, nu = 30)
  )) {
    expect_equal(
      stats::integrate(function(t) distrib_pdf(d, t, th), -Inf, Inf)$value,
      1,
      tolerance = 1e-6,
      label = paste("integral, nu =", th$nu)
    )
  }
})

test_that("pseudohuber gradient and hessian match finite differences", {
  set.seed(31)
  d <- pseudohuber_distrib()
  th <- list(mu = 0.5, sigma = 1.4, nu = 2.5)
  y <- distrib_rng(d, 10, th)

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
})

test_that("pseudohuber cdf/quantile/rng are consistent", {
  d <- pseudohuber_distrib()
  th <- list(mu = 0.5, sigma = 1.4, nu = 2.5)

  # median is mu by symmetry
  expect_equal(distrib_quantile(d, 0.5, th), 0.5)

  p <- c(0.05, 0.3, 0.7, 0.99)
  q <- distrib_quantile(d, p, th)
  expect_equal(distrib_cdf(d, q, th), p, tolerance = 1e-7)

  # symmetry of quantiles around mu
  expect_equal(q[1] - 0.5, -(distrib_quantile(d, 0.95, th) - 0.5), tolerance = 1e-7)

  set.seed(5)
  y <- distrib_rng(d, 200, th)
  expect_length(y, 200)
  expect_true(all(is.finite(y)))
})

test_that("pseudohuber expected hessian has the symmetry zeros and negative-definite diagonal", {
  d <- pseudohuber_distrib()
  th <- list(mu = 0.5, sigma = 1.4, nu = 2.5)
  eh <- distrib_expected_hessian(d, c(0, 1), th)

  expect_named(eh, c("mu_mu", "sigma_sigma", "nu_nu", "mu_sigma", "mu_nu", "sigma_nu"))
  expect_equal(eh$mu_sigma, c(0, 0))
  expect_equal(eh$mu_nu, c(0, 0))
  expect_true(all(eh$mu_mu < 0))
  expect_true(all(eh$sigma_sigma < 0))
  expect_true(all(eh$nu_nu < 0))
})
