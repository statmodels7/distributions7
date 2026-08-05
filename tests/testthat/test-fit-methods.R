# Transformed distributions end to end, plus the simulate() and plot() methods
# of a fitted model.

test_that("log(lognormal) is the gaussian, all the way through", {
  # lognormal(mu, sigma2) has log Y ~ N(mu, sigma2), so the log transformation
  # must reproduce the gaussian exactly -- a sharp end-to-end check of the
  # change-of-variables machinery against a distribution implemented separately.
  lg <- transformation(lognormal1_distrib(), log_transform())
  g <- gaussian1_distrib()
  th <- list(mu = 0.7, sigma2 = 1.3)
  thg <- list(mu = 0.7, sigma = sqrt(1.3))

  z <- c(-2, -0.5, 0, 1.2, 3)
  p <- c(0.05, 0.25, 0.5, 0.9)

  expect_equal(distrib_pdf(lg, z, th), distrib_pdf(g, z, thg))
  expect_equal(distrib_pdf(lg, z, th, log = TRUE), distrib_pdf(g, z, thg, log = TRUE))
  expect_equal(distrib_cdf(lg, z, th), distrib_cdf(g, z, thg))
  expect_equal(distrib_quantile(lg, p, th), distrib_quantile(g, p, thg), tolerance = 1e-7)

  # derivatives are w.r.t. sigma2 rather than sigma, so compare against FD
  a <- distrib_gradient(lg, z, th)
  e <- numerical_gradient(lg, z, th)
  for (k in names(a)) expect_equal(a[[k]], e[[k]], tolerance = 1e-5, label = k)
})

test_that("fitting log(lognormal) and the gaussian to one sample agrees exactly", {
  set.seed(71)
  y <- stats::rnorm(1500, 0.7, sqrt(1.3))
  f_lg <- fit_distrib(transformation(lognormal1_distrib(), log_transform()), y)
  f_g <- fit_distrib(gaussian1_distrib(), y)

  expect_equal(coef(f_lg)[["mu"]], coef(f_g)[["mu"]], tolerance = 1e-6)
  expect_equal(sqrt(coef(f_lg)[["sigma2"]]), coef(f_g)[["sigma"]], tolerance = 1e-6)
  expect_equal(as.numeric(logLik(f_lg)), as.numeric(logLik(f_g)), tolerance = 1e-8)
})

test_that("transformed distributions pass their own validator and recover parameters", {
  cases <- list(
    log_gamma = list(d = transformation(gamma2_distrib(), log_transform()),
                     th = list(mu = 3, sigma2 = 2)),
    exp_gaussian = list(d = transformation(gaussian1_distrib(), exp_transform()),
                        th = list(mu = 0.4, sigma = 1.1)),
    logit_beta = list(d = transformation(beta1_distrib(), logit_transform()),
                      th = list(mu = 0.4, phi = 6))
  )
  for (nm in names(cases)) {
    d <- cases[[nm]]$d
    th <- cases[[nm]]$th
    set.seed(72)
    res <- check_distrib(d, th, n = 30, nsim = 2e4, orders = 1:2, verbose = FALSE)
    expect_true(all(res$status == "OK"),
      label = paste(nm, "-", paste(res$check[res$status != "OK"], collapse = "; ")))

    set.seed(73)
    f <- fit_distrib(d, distrib_rng(d, 2000, th))
    se <- sqrt(diag(vcov(f)))
    for (p in names(th)) {
      expect_lt(abs(coef(f)[[p]] - th[[p]]) / se[[p]], 4)
    }
  }
})

test_that("simulate() follows the stats::simulate contract", {
  set.seed(74)
  y <- stats::rnorm(200, 3, 2)
  fit <- fit_distrib(gaussian1_distrib(), y)

  s <- simulate(fit, 5)
  expect_s3_class(s, "data.frame")
  expect_named(s, paste0("sim_", 1:5))
  expect_equal(nrow(s), fit@n)

  # a supplied seed makes it reproducible ...
  a <- simulate(fit, 3, seed = 42)
  b <- simulate(fit, 3, seed = 42)
  expect_equal(a, b)
  # stats::simulate reports the seed with the RNG kind attached
  expect_equal(as.vector(attr(a, "seed")), 42)
  expect_named(attr(attr(a, "seed"), "kind"), NULL)

  # ... and leaves the caller's random stream untouched
  set.seed(99)
  before <- stats::runif(1)
  set.seed(99)
  invisible(simulate(fit, 2, seed = 7))
  expect_equal(stats::runif(1), before)

  expect_error(simulate(fit, 0), "positive integer")
})

test_that("simulate() draws from the fitted distribution", {
  set.seed(75)
  y <- stats::rgamma(500, shape = 4, rate = 2)
  fit <- fit_distrib(gamma2_distrib(), y)
  s <- simulate(fit, 1, seed = 3)[[1]]

  expect_length(s, 500)
  expect_true(all(s > 0))
  # the draws should look like the fitted distribution, not like anything else
  ks <- suppressWarnings(stats::ks.test(
    s, function(q) distrib_cdf(fit@distrib, q, as.list(coef(fit)))
  ))
  expect_gt(ks$p.value, 0.001)
})

test_that("simulate() works for a discrete and a transformed fit", {
  set.seed(76)
  fp <- fit_distrib(poisson_distrib(), stats::rpois(300, 4))
  sp <- simulate(fp, 2, seed = 1)
  expect_equal(dim(sp), c(300L, 2L))
  expect_true(all(sp[[1]] == floor(sp[[1]])))

  d <- transformation(gamma2_distrib(), log_transform())
  set.seed(77)
  ft <- fit_distrib(d, distrib_rng(d, 300, list(mu = 3, sigma2 = 2)))
  st <- simulate(ft, 2, seed = 1)
  expect_equal(dim(st), c(300L, 2L))
  expect_true(all(is.finite(as.matrix(st))))
})

test_that("plot() draws for continuous, discrete and transformed fits", {
  skip_if_not(capabilities("png"))
  tmp <- tempfile(fileext = ".png")
  grDevices::png(tmp)
  on.exit({
    grDevices::dev.off()
    unlink(tmp)
  }, add = TRUE)

  set.seed(78)
  f1 <- fit_distrib(gaussian1_distrib(), stats::rnorm(200, 3, 2))
  expect_identical(plot(f1), f1)

  f2 <- fit_distrib(poisson_distrib(), stats::rpois(200, 4))
  expect_identical(plot(f2), f2)

  d <- transformation(gamma2_distrib(), log_transform())
  f3 <- fit_distrib(d, distrib_rng(d, 200, list(mu = 3, sigma2 = 2)))
  expect_identical(plot(f3), f3)

  # a bounded support: the fit must not be evaluated outside it
  f4 <- fit_distrib(beta1_distrib(), distrib_rng(beta1_distrib(), 200, list(mu = 0.4, phi = 6)))
  expect_identical(plot(f4, legend = FALSE, rug = FALSE), f4)
})

test_that("the fit keeps the data it was estimated from", {
  set.seed(79)
  y <- stats::rnorm(50)
  f <- fit_distrib(gaussian1_distrib(), y)
  expect_equal(f@y, y)
  expect_equal(f@n, 50)
})
