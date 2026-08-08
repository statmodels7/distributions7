# Laplace in location and rate: the same law as laplace_distrib at
# sigma = 1/lambda, so the two implementations check each other.

test_that("laplace2 pdf/cdf/quantile/rng/moments are correct", {
  d <- laplace2_distrib()
  th <- list(mu = 1, lambda = 0.5)

  expect_equal(stats::integrate(function(t) distrib_pdf(d, t, th), -Inf, Inf)$value, 1, tolerance = 1e-6)
  expect_equal(distrib_pdf(d, 1, th), 0.5 / 2)
  p <- c(0.05, 0.25, 0.5, 0.75, 0.95)
  expect_equal(distrib_quantile(d, 0.5, th), 1)
  expect_equal(distrib_cdf(d, distrib_quantile(d, p, th), th), p, tolerance = 1e-9)

  expect_equal(mean(d, th), 1)
  expect_equal(variance(d, th), 2 / 0.5^2)
  expect_equal(skewness(d, th), 0)
  expect_equal(kurtosis(d, th), 3)

  set.seed(1)
  y <- distrib_rng(d, 1e5, th)
  expect_equal(mean(y), 1, tolerance = 0.05)
  expect_equal(stats::var(y), 8, tolerance = 0.3)
})

test_that("laplace2 is laplace at sigma = 1/lambda, derivatives by the chain rule", {
  d2 <- laplace2_distrib()
  d1 <- laplace_distrib()
  lam <- 0.8
  th2 <- list(mu = 0.3, lambda = lam)
  th1 <- list(mu = 0.3, sigma = 1 / lam)
  y <- c(-2, 0.4, 1.7, 5)
  q <- c(-1, 0.3, 2)
  p <- c(0.1, 0.5, 0.9)

  expect_equal(distrib_pdf(d2, y, th2), distrib_pdf(d1, y, th1))
  expect_equal(distrib_cdf(d2, q, th2), distrib_cdf(d1, q, th1))
  expect_equal(distrib_quantile(d2, p, th2), distrib_quantile(d1, p, th1))
  expect_equal(distrib_grad_y(d2, y, th2), distrib_grad_y(d1, y, th1))
  expect_equal(distrib_hess_y(d2, y, th2), distrib_hess_y(d1, y, th1))

  # score and Hessian: sigma = 1/lambda, so dsigma/dlambda = -1/lambda^2 and
  # d2sigma/dlambda2 = 2/lambda^3
  J <- -1 / lam^2
  K <- 2 / lam^3
  g2 <- distrib_gradient(d2, y, th2)
  g1 <- distrib_gradient(d1, y, th1)
  expect_equal(g2$mu, g1$mu)
  expect_equal(g2$lambda, g1$sigma * J)

  h2 <- distrib_hessian(d2, y, th2)
  h1 <- distrib_hessian(d1, y, th1)
  expect_equal(h2$mu_mu, h1$mu_mu)
  expect_equal(h2$mu_lambda, h1$mu_sigma * J)
  expect_equal(h2$lambda_lambda, h1$sigma_sigma * J^2 + g1$sigma * K)

  # expected Hessian: the information transforms by the same congruence
  eh2 <- distrib_expected_hessian(d2, y, th2)
  eh1 <- distrib_expected_hessian(d1, y, th1)
  expect_equal(eh2$mu_mu, eh1$mu_mu)
  expect_equal(eh2$mu_lambda, eh1$mu_sigma * J)
  expect_equal(eh2$lambda_lambda, eh1$sigma_sigma * J^2)
})

test_that("laplace2 higher orders: only the pure-lambda components survive", {
  d <- laplace2_distrib()
  lam <- 1.6
  th <- list(mu = 0, lambda = lam)
  y <- c(-1.2, 0.7, 2.1)

  d3 <- distrib_deriv3(d, y, th)
  expect_equal(d3[["lambda_lambda_lambda"]], rep(2 / lam^3, 3))
  expect_equal(d3[["mu_mu_mu"]], rep(0, 3))
  expect_equal(d3[["mu_mu_lambda"]], rep(0, 3))
  expect_equal(d3[["mu_lambda_lambda"]], rep(0, 3))

  d4 <- distrib_deriv4(d, y, th)
  expect_equal(d4[["lambda_lambda_lambda_lambda"]], rep(-6 / lam^4, 3))
  expect_equal(d4[["mu_lambda_lambda_lambda"]], rep(0, 3))

  # data-free in lambda beyond the first order, so observed = expected there
  e3 <- distrib_deriv3(d, y, th, expected = TRUE)
  expect_equal(e3[["lambda_lambda_lambda"]], d3[["lambda_lambda_lambda"]])
  e4 <- distrib_deriv4(d, y, th, expected = TRUE)
  expect_equal(e4[["lambda_lambda_lambda_lambda"]], d4[["lambda_lambda_lambda_lambda"]])

  # the third derivative in lambda against a finite difference of the closed
  # Hessian (legal: one stencil on an analytic quantity)
  h <- 1e-5
  fd <- (distrib_hessian(d, y, list(mu = 0, lambda = lam + h))$lambda_lambda -
           distrib_hessian(d, y, list(mu = 0, lambda = lam - h))$lambda_lambda) / (2 * h)
  expect_equal(d3[["lambda_lambda_lambda"]], fd, tolerance = 1e-6)
})

test_that("laplace2 mixed response-parameter block matches a difference of grad_y", {
  d <- laplace2_distrib()
  th <- list(mu = 0.3, lambda = 1.4)
  y <- c(-2, 0.8, 3.1) # away from the kink
  cr <- distrib_cross_y(d, y, th)
  h <- 1e-6
  fd_lam <- (distrib_grad_y(d, y, list(mu = 0.3, lambda = 1.4 + h)) -
               distrib_grad_y(d, y, list(mu = 0.3, lambda = 1.4 - h))) / (2 * h)
  expect_equal(cr$lambda, fd_lam, tolerance = 1e-6)
  expect_equal(cr$mu, rep(0, 3))
})

test_that("laplace2 cdf derivatives match finite differences of the cdf", {
  d <- laplace2_distrib()
  th <- list(mu = 0.5, lambda = 1.2)
  q <- c(-1.5, 0.1, 2.4)
  h <- 1e-6

  fd1 <- function(j) {
    tp <- tm <- th; tp[[j]] <- th[[j]] + h; tm[[j]] <- th[[j]] - h
    (distrib_cdf(d, q, tp) - distrib_cdf(d, q, tm)) / (2 * h)
  }
  g <- distrib_grad_cdf(d, q, th, log = FALSE)
  expect_equal(g$mu, fd1(1), tolerance = 1e-7)
  expect_equal(g$lambda, fd1(2), tolerance = 1e-7)

  fd2 <- function(j, k) {
    tp <- tm <- th; tp[[k]] <- th[[k]] + h; tm[[k]] <- th[[k]] - h
    (distrib_grad_cdf(d, q, tp, log = FALSE)[[j]] -
       distrib_grad_cdf(d, q, tm, log = FALSE)[[j]]) / (2 * h)
  }
  hs <- distrib_hess_cdf(d, q, th, log = FALSE)
  expect_equal(hs$mu_mu, fd2("mu", 1), tolerance = 1e-5)
  expect_equal(hs$lambda_lambda, fd2("lambda", 2), tolerance = 1e-5)
  expect_equal(hs$mu_lambda, fd2("mu", 2), tolerance = 1e-5)
})

test_that("laplace2 passes the check_distrib battery", {
  set.seed(3)
  res <- check_distrib(laplace2_distrib(), theta = list(mu = 1, lambda = 0.5),
                       n = 40, nsim = 5e4, orders = 1:2, verbose = FALSE)
  expect_true(all(res$status == "OK"),
              label = paste(res$check[res$status != "OK"], collapse = "; "))
})

test_that("laplace2 records the kink and fits by Fisher scoring", {
  d <- laplace2_distrib()
  expect_equal(param_smoothness(d), c(mu = FALSE, lambda = TRUE))

  set.seed(4)
  y <- distrib_rng(d, 4000, list(mu = 2, lambda = 1.5))
  f <- fit_distrib(d, y)
  expect_true(f@converged)
  expect_equal(unname(coef(f)["mu"]), 2, tolerance = 0.1)
  expect_equal(unname(coef(f)["lambda"]), 1.5, tolerance = 0.1)
})
