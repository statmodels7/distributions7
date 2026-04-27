test_that("gaussian_distrib standard functions match base R", {
  d <- gaussian_distrib()
  theta <- list(mu = 2, sigma = 1.5)
  y <- c(-1, 0, 2, 5)
  p <- c(0.1, 0.5, 0.9)
  
  # 1. PDF
  expect_equal(
    distrib_pdf(d, y, theta),
    stats::dnorm(y, mean = theta$mu, sd = theta$sigma)
  )
  
  # 2. CDF
  expect_equal(
    distrib_cdf(d, y, theta),
    stats::pnorm(y, mean = theta$mu, sd = theta$sigma)
  )
  
  # 3. Quantile
  expect_equal(
    distrib_quantile(d, p, theta),
    stats::qnorm(p, mean = theta$mu, sd = theta$sigma)
  )
  
  # 4. RNG
  set.seed(42)
  sims <- distrib_rng(d, 100, theta)
  expect_length(sims, 100)
  expect_true(is.numeric(sims))
})

test_that("gaussian_distrib derivatives are structurally correct and mathematically sound", {
  d <- gaussian_distrib()
  
  # Test exactly at the mean (y = mu)
  y_val <- 5
  theta <- list(mu = 5, sigma = 2)
  
  # Gradient
  grad <- distrib_gradient(d, y_val, theta)
  expect_named(grad, c("mu", "sigma"))
  # The gradient w.r.t mu evaluated at the mean should be exactly 0
  expect_equal(grad$mu, 0)
  
  # Hessian
  hess <- distrib_hessian(d, y_val, theta)
  expect_named(hess, c("mu_mu", "sigma_sigma", "mu_sigma"))
  expect_true(all(!is.na(unlist(hess))))
})