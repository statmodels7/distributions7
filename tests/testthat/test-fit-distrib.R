# Maximum-likelihood estimation on the link scale, checked against closed-form
# MLEs where they exist.

test_that("gaussian fit reproduces the closed-form MLE and its standard errors", {
  set.seed(42)
  d <- gaussian_distrib()
  y <- distrib_rng(d, 2000, list(mu = 2, sigma = 3))
  f <- fit_distrib(d, y, start = list(mu = 0, sigma = 1))

  expect_true(f@converged)
  expect_equal(unname(coef(f)["mu"]), mean(y), tolerance = 1e-6)
  expect_equal(unname(coef(f)["sigma"]), sqrt(mean((y - mean(y))^2)), tolerance = 1e-6)

  # se(mu_hat) = sigma / sqrt(n)
  expect_equal(unname(f@se["mu"]), unname(coef(f)["sigma"]) / sqrt(length(y)), tolerance = 1e-4)
  expect_equal(as.numeric(logLik(f)), sum(distrib_pdf(d, y, as.list(coef(f)), log = TRUE)))
})

test_that("poisson fit reproduces the closed-form MLE and information", {
  set.seed(7)
  d <- poisson_distrib()
  y <- distrib_rng(d, 1000, list(mu = 4.5))
  f <- fit_distrib(d, y, start = list(mu = 1))

  expect_true(f@converged)
  expect_equal(unname(coef(f)), mean(y), tolerance = 1e-7)
  expect_equal(unname(f@se), sqrt(mean(y) / length(y)), tolerance = 1e-5)
})

test_that("bernoulli fit is exact and its interval respects the (0,1) domain", {
  set.seed(11)
  d <- bernoulli_distrib()
  y <- rbinom(60, 1, 0.92)
  f <- fit_distrib(d, y, start = list(mu = 0.5))

  expect_true(f@converged)
  expect_equal(unname(coef(f)), mean(y), tolerance = 1e-7)
  expect_true(all(f@ci > 0 & f@ci < 1))
  expect_lt(f@ci[1, "lower"], coef(f)[[1]])
  expect_gt(f@ci[1, "upper"], coef(f)[[1]])
})

test_that("laplace fit works despite the degenerate observed Hessian in mu", {
  set.seed(13)
  d <- laplace_distrib()
  y <- distrib_rng(d, 2000, list(mu = 1, b = 2))
  f <- fit_distrib(d, y, start = list(mu = 0, b = 1))

  expect_true(f@converged)
  # b_hat is the mean absolute deviation about the fitted location
  expect_equal(unname(coef(f)["b"]), mean(abs(y - coef(f)[["mu"]])), tolerance = 1e-4)
  # the location is at (essentially) the sample median
  expect_equal(unname(coef(f)["mu"]), median(y), tolerance = 0.05)
  expect_true(all(f@se > 0))
})

test_that("all optimisation methods agree", {
  set.seed(3)
  d <- gaussian_distrib()
  y <- distrib_rng(d, 800, list(mu = 1, sigma = 2))
  fits <- lapply(c("fisher", "newton", "bfgs"), function(m) {
    fit_distrib(d, y, start = list(mu = 0, sigma = 1), method = m)
  })
  for (f in fits) expect_true(f@converged)
  expect_equal(coef(fits[[1]]), coef(fits[[2]]), tolerance = 1e-5)
  expect_equal(coef(fits[[1]]), coef(fits[[3]]), tolerance = 1e-4)
})

test_that("extractors and the fitted object are coherent", {
  set.seed(5)
  d <- gamma_distrib()
  y <- distrib_rng(d, 800, list(mu = 3, sigma2 = 2))
  f <- fit_distrib(d, y, start = list(mu = 1, sigma2 = 1))

  expect_s3_class(logLik(f), "logLik")
  expect_equal(attr(logLik(f), "nobs"), length(y))
  expect_equal(dim(vcov(f)), c(2L, 2L))
  expect_equal(dim(vcov(f, scale = "link")), c(2L, 2L))
  expect_named(coef(f), d@params)
  expect_named(coef(f, scale = "link"), d@params)

  # delta method: V_theta = diag(h') V_eta diag(h')
  J <- vapply(seq_along(d@params), function(i) {
    linkfunctions7::linkinvderiv(d@link_params[[d@params[i]]], f@eta[i], order = 1)
  }, numeric(1))
  expect_equal(vcov(f), diag(J) %*% vcov(f, scale = "link") %*% diag(J),
               ignore_attr = TRUE, tolerance = 1e-10)

  # AIC / BIC definitions
  expect_equal(f@aic, -2 * f@loglik + 2 * length(coef(f)))
  expect_equal(f@bic, -2 * f@loglik + log(f@n) * length(coef(f)))
})

test_that("fitting works from the default random starting values", {
  set.seed(9)
  d <- gaussian_distrib()
  y <- distrib_rng(d, 500, list(mu = 2, sigma = 3))
  f <- fit_distrib(d, y)
  expect_true(f@converged)
  expect_equal(unname(coef(f)["mu"]), mean(y), tolerance = 1e-4)
})

test_that("the print method reports both scales", {
  set.seed(2)
  d <- gaussian_distrib()
  y <- distrib_rng(d, 200, list(mu = 0, sigma = 1))
  f <- fit_distrib(d, y, start = list(mu = 0, sigma = 1))
  out <- paste(utils::capture.output(print(f)), collapse = "\n")
  expect_match(out, "Maximum-likelihood fit")
  expect_match(out, "Parameter scale")
  expect_match(out, "Link scale")
  expect_match(out, "Fisher scoring")
})
