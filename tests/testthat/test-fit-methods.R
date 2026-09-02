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

test_that("the information at the optimum does not depend on the SPELLING", {
  # `method` takes a string or an optimizer object, and the choice between
  # the observed and the expected information was made by comparing it with
  # the string "newton". An object is normalized to "custom", so newton()
  # took the expected branch -- which for a family that does not write the
  # expectation out is a quadrature, and on a family whose score grows
  # exponentially that quadrature does not converge. The fit took 0.15 s and
  # the standard errors never returned.
  skip_on_cran()
  Gomp <- S7::new_class("Gomp", parent = continuous_distrib)
  S7::method(distrib_pdf, Gomp) <- function(distrib, y, theta, log = FALSE) {
    eta <- theta[[1]]; b <- theta[[2]]
    ld <- log(eta) + log(b) + b * y - eta * expm1(b * y)
    ld[y < 0] <- -Inf
    if (log) ld else exp(ld)
  }
  gz <- Gomp(distrib_name = "gomp", dimension = "univariate",
             bounds = c(0, Inf), params = c("eta", "b"),
             params_interpretation = c(eta = "level", b = "rate"),
             n_params = 2,
             params_bounds = list(eta = c(0, Inf), b = c(0, Inf)),
             link_params = list(eta = linkfunctions7::log_link(),
                                b = linkfunctions7::log_link()))
  expect_false(has_exact_expected_hessian(gz))

  set.seed(2)
  yg <- distrib_rng(gz, 150, list(eta = 0.02, b = 0.9))

  t0 <- proc.time()[["elapsed"]]
  f_obj <- fit_distrib(gz, yg, method = optimizers7::newton())
  el <- proc.time()[["elapsed"]] - t0
  f_str <- fit_distrib(gz, yg, method = "newton")

  # the estimate is the same algorithm's, and so is the variance matrix
  expect_equal(unlist(coef(f_obj)), unlist(coef(f_str)), tolerance = 1e-6)
  expect_equal(vcov(f_obj), vcov(f_str), tolerance = 1e-6)
  # and it returns at all, which is the whole point: the divergent
  # quadrature ran for minutes and raised nothing, so no tryCatch caught it
  expect_lt(el, 60)
  expect_true(all(is.finite(vcov(f_obj))))

  # a family that DOES write the expectation out is unaffected: Fisher
  # scoring still reports the expected information there
  gg <- gamma1_distrib()
  expect_true(has_exact_expected_hessian(gg))
  set.seed(3)
  yy <- distrib_rng(gg, 400, list(mu = 3, phi = 2))
  fg <- fit_distrib(gg, yy)
  # on the LINK scale, which is where the information is assembled; vcov()
  # defaults to the parameter scale and is that matrix carried through the
  # delta method
  expect_equal(vcov(fg, scale = "link"),
               solve(-fit_hess_matrix(gg, yy,
                                      as.list(coef(fg, scale = "parameter")),
                                      expected = TRUE)),
               tolerance = 1e-6, ignore_attr = TRUE)
})

test_that("distrib_kernel agrees with the generic route it bypasses", {
  # The kernel resolves the family's methods and the link's once and applies
  # the chain rule for ONE component; the generic validates, aligns, and
  # assembles every component of its order. They must not differ.
  fams <- list(
    list(d = gaussian1_distrib(), th = list(mu = 0.3, sigma = 1.4)),
    list(d = gamma1_distrib(),    th = list(mu = 2, phi = 3)),
    list(d = poisson_distrib(),   th = list(mu = 2.5)),
    list(d = negbin2_distrib(),   th = list(mu = 3, theta = 2)),
    list(d = weibull1_distrib(),  th = list(mu = 2, sigma = 1.5)),
    list(d = gpd_distrib(),       th = list(sigma = 1.2, xi = 0.3))
  )
  for (f in fams) {
    set.seed(5)
    y <- distrib_rng(f$d, 40, f$th)
    for (p in f$d@params) {
      lk <- f$d@link_params[[p]]
      eta <- linkfunctions7::linkfun(lk, f$th[[p]])
      k <- distrib_kernel(f$d, p)
      # the entry for p on the way in is immaterial, being replaced
      th_in <- f$th
      th_in[[p]] <- NA_real_
      expect_equal(as.numeric(k$logdens(y, th_in, eta)),
                   as.numeric(rep_len(distrib_pdf(f$d, y, f$th, log = TRUE),
                                      length(y))), tolerance = 1e-10)
      expect_equal(as.numeric(k$score(y, th_in, eta)),
                   as.numeric(rep_len(distrib_gradient(
                     f$d, y, f$th, scale = "link")[[p]], length(y))),
                   tolerance = 1e-10)
      expect_equal(as.numeric(k$curvature(y, th_in, eta)),
                   as.numeric(rep_len(distrib_hessian(
                     f$d, y, f$th, scale = "link")[[paste0(p, "_", p)]],
                     length(y))), tolerance = 1e-10)
    }
  }

  # a vector eta, which a filter never uses and a caller may
  d <- gaussian1_distrib()
  set.seed(6)
  y <- stats::rnorm(30)
  ev <- stats::rnorm(30, 0, 0.3)
  k <- distrib_kernel(d, "sigma")
  expect_equal(k$score(y, list(mu = 0.2, sigma = NA), ev),
               distrib_gradient(d, y, list(mu = 0.2, sigma = exp(ev)),
                                scale = "link")[["sigma"]],
               ignore_attr = TRUE)

  # the inverse link is clamped strictly inside its bounds, which is what
  # linkinv()'s own generic body does and is a correctness property rather
  # than something the fast path may skip: an unclamped exp(-800) is zero,
  # and a gaussian with a scale of exactly zero is not a distribution
  expect_true(is.finite(k$logdens(0, list(mu = 0, sigma = NA), -800)))

  expect_error(distrib_kernel(d, "nope"), "not a parameter")
  expect_error(distrib_kernel(mvgaussian1_distrib(2), "mu1"), "univariate")
})

test_that("every moment estimate inverts its own family's moments", {
  # The check that matters: simulate a large sample from a KNOWN theta with
  # the family's own rng and require the moment estimate to return it. That
  # tests the inversion against the family, not against memory, and it is
  # what catches a variance function written from the wrong parametrization.
  skip_on_cran()
  cases <- list(
    list(gaussian1_distrib(),  list(mu = 30, sigma = 7)),
    list(gaussian2_distrib(),  list(mu = 30, sigma2 = 49)),
    list(gaussian3_distrib(),  list(mu = 30, tau = 1 / 49)),
    list(cauchy_distrib(),     list(mu = 30, sigma = 7)),
    list(laplace_distrib(),    list(mu = 30, sigma = 7)),
    list(laplace2_distrib(),   list(mu = 30, lambda = 1 / 7)),
    list(logistic_distrib(),   list(mu = 30, sigma = 7)),
    list(student_t1_distrib(), list(mu = 30, sigma = 7, nu = 8)),
    list(student_t2_distrib(), list(mu = 30, sigma = 7, nu = 8)),
    list(gamma1_distrib(),     list(mu = 30, phi = 0.5)),
    list(gamma2_distrib(),     list(mu = 30, sigma2 = 400)),
    list(exponential_distrib(), list(mu = 30)),
    list(chisq_distrib(),      list(mu = 9)),
    list(poisson_distrib(),    list(mu = 40)),
    list(geometric_distrib(),  list(mu = 12)),
    list(negbin1_distrib(),    list(mu = 40, theta = 3)),
    list(negbin2_distrib(),    list(mu = 40, theta = 3)),
    list(beta1_distrib(),      list(mu = 0.35, phi = 9)),
    list(beta2_distrib(),      list(alpha = 3, beta = 5)),
    list(bernoulli_distrib(),  list(mu = 0.35)),
    list(weibull1_distrib(),   list(mu = 300, sigma = 2.5)),
    list(lognormal1_distrib(), list(mu = 4, sigma2 = 0.3)),
    list(lognormal2_distrib(), list(mean = 60, var = 900)),
    list(invgauss1_distrib(),  list(mu = 30, phi = 0.02)),
    list(invgauss2_distrib(),  list(mu = 30, lambda = 400)),
    list(gumbel_distrib(),     list(mu = 30, sigma = 7)),
    list(gpd_distrib(),        list(sigma = 20, xi = 0.2)),
    list(pig1_distrib(),       list(mu = 40, sigma = 0.05)),
    list(vonmises1_distrib(),  list(mu = 0.7, kappa = 3)),
    # the second chart of a family whose first one has an estimate: the
    # inversion is carried across the map rather than derived again
    list(weibull3_distrib(),   list(mean = 300, sigma = 2.5)),
    list(pig2_distrib(),       list(mu = 40, alpha = 3)),
    list(betabinom2_distrib(size = 10), list(alpha = 3, beta = 5)),
    # three moments, and the skew normal reaches its shape from the third
    list(skewnormal2_distrib(), list(mu = 3, sigma = 2, gamma1 = 0.4)),
    list(skewnormal1_distrib(), list(mu = 3, sigma = 2, alpha = 1.5)),
    list(betabinom1_distrib(size = 10), list(mu = 0.35, sigma = 0.25))
  )
  for (cs in cases) {
    d <- cs[[1L]]
    th <- cs[[2L]]
    set.seed(99)
    y <- distrib_rng(d, 200000L, th)
    e <- moment_estimates(d, y)
    expect_false(is.null(e), label = d@distrib_name)
    expect_named(e, d@params)
    rel <- max(abs(unlist(e) - unlist(th)) / pmax(abs(unlist(th)), 1e-3))
    # loose, because it is a starting value and a t's degrees of freedom
    # come from a fourth moment; the point is that it is an ESTIMATE
    expect_lt(rel, 0.1)
  }

  # A family whose moments cannot be inverted in closed form falls back to the
  # interpretation route rather than failing. Five do: the skew t and the
  # elastic net, whose systems are four equations in quantities carrying a
  # distribution function and a Mills ratio, the two generalized gammas, whose
  # moments are ratios of gamma functions in two shapes, and the pseudo-Huber,
  # whose moments have no elementary form.
  for (d in list(skewt_distrib(), pseudohuber_distrib(), gengamma1_distrib(),
                 gengamma2_distrib(), enet_distrib())) {
    expect_null(moment_estimates(d, stats::rnorm(100)), label = d@distrib_name)
    s <- distrib_start(d, stats::rnorm(100), 2L)
    expect_length(s, 2L)
    expect_named(s[[1L]], d@params)
  }

  # A family carrying a fixed constant spells it in its name, which must not
  # reach the lookup key: without dropping it, no beta-binomial of any size
  # would ever match its own entry.
  expect_false(is.null(moment_estimates(betabinom1_distrib(size = 7),
                                        stats::rbinom(200, 7, 0.4))))
})

test_that("a fit recovers its parameters across four decades of scale", {
  # what the moment start is FOR: the old random draw ignored the data, so
  # a response of order 1000 sent the scale to the largest double
  skip_on_cran()
  for (mult in c(1, 10, 100, 1000, 10000)) {
    set.seed(7)
    z <- stats::rnorm(300, 5 * mult, 2 * mult)
    f <- suppressWarnings(fit_distrib(gaussian1_distrib(), z, n_start = 5L))
    expect_true(f@converged)
    expect_equal(coef(f)[["mu"]], 5 * mult, tolerance = 0.05)
    expect_equal(coef(f)[["sigma"]], 2 * mult, tolerance = 0.05)
  }
})
