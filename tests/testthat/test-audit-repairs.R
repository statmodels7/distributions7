# Repairs found by the documentation audit of 2026-10-11. Each block fails on
# 0.76.1 and passes on the repaired code.

test_that("the pairs panels of a fit read the data by coordinate", {
  set.seed(2)
  X <- cbind(rnorm(100), rnorm(100, 100), rnorm(100, 5000))
  fm <- fit_distrib(mvgaussian1_distrib(3), X)
  rec <- new.env()
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  suppressMessages(trace(graphics::points, bquote(assign(
    "xs", c(get0("xs", envir = .(rec)), range(x)), envir = .(rec))),
    print = FALSE))
  on.exit(suppressMessages(untrace(graphics::points)), add = TRUE)
  plot(fm, mv_which = c(2, 3))
  # the one scatter panel has coordinate 2 on its horizontal axis
  expect_true(all(rec$xs > 90 & rec$xs < 110))
})

test_that("a fit over four coordinates names mv_which in its error", {
  f4 <- fit_distrib(mvgaussian1_distrib(4), as.matrix(iris[, 1:4]))
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  expect_error(plot(f4), "mv_which = c(1, 2, 3)", fixed = TRUE)
})

test_that("the elapsed time never reaches the carry of its unit", {
  expect_identical(fit_format_elapsed(0.9996), "1 s")
  expect_identical(fit_format_elapsed(59.96), "1 min 00 s")
  expect_identical(fit_format_elapsed(119.7), "2 min 00 s")
  expect_identical(fit_format_elapsed(125), "2 min 05 s")
  expect_identical(fit_format_elapsed(0.0123), "12.3 ms")
})

test_that("the discrete numerical quantile inverts the exact cdf", {
  PoisD <- S7::new_class("PoisDAudit", parent = discrete_distrib)
  S7::method(distrib_pdf, PoisD) <- function(distrib, y, theta, log = FALSE, ...) {
    stats::dpois(y, theta[[1]], log = log)
  }
  pd <- PoisD(distrib_name = "poisd", dimension = "univariate",
              bounds = c(0, Inf), params = "mu",
              params_interpretation = c(mu = "mean"), n_params = 1L,
              params_bounds = list(mu = c(0, Inf)),
              link_params = list(mu = linkfunctions7::log_link()),
              params_smooth = c(mu = TRUE))
  for (lam in c(3, 20, 200)) {
    k <- 0:stats::qpois(1 - 1e-9, lam)
    p <- stats::ppois(k, lam)
    expect_identical(distrib_quantile(pd, p, list(mu = lam)), as.numeric(k))
    expect_identical(distrib_quantile(pd, p, list(mu = rep(lam, length(p)))),
                     as.numeric(k))
    mid <- (p[-1] + p[-length(p)]) / 2
    expect_identical(distrib_quantile(pd, mid, list(mu = lam)),
                     as.numeric(k[-1]))
  }
})

test_that("cdf derivatives are zero at and below zero on a positive support", {
  q <- c(-1, 0, 0.5)
  cases <- list(list(lognormal1_distrib(), list(mu = 0.2, sigma2 = 0.7)),
                list(weibull1_distrib(), list(mu = 2, sigma = 3)))
  gens <- list(distrib_grad_cdf, distrib_hess_cdf, distrib_deriv3_cdf,
               distrib_deriv4_cdf)
  for (cs in cases) for (g in gens) for (lt in c(TRUE, FALSE)) {
    v <- expect_no_warning(g(cs[[1]], q, cs[[2]], lower.tail = lt, log = FALSE))
    for (nm in names(v)) {
      expect_identical(v[[nm]][1:2], c(0, 0))
      expect_true(is.finite(v[[nm]][3]))
    }
  }
})

test_that("numerical_deriv_y() rejects an order it does not compute", {
  d <- gaussian1_distrib()
  th <- list(mu = 0, sigma = 1)
  expect_error(numerical_deriv_y(d, 0, th, order = 2), "must be 3 or 4")
  expect_error(numerical_deriv_y(d, 0, th, order = 5), "must be 3 or 4")
  expect_length(numerical_deriv_y(d, c(-1, 1), th, order = 3), 2L)
})

test_that("a confidence level outside (0, 1) is rejected", {
  set.seed(1)
  y <- rnorm(50)
  expect_error(fit_distrib(gaussian1_distrib(), y, level = 95), "'level'")
  f <- fit_distrib(gaussian1_distrib(), y)
  expect_error(confint(f, level = 95), "'level'")
  expect_error(confint(f, level = NA), "'level'")
  expect_identical(dim(confint(f, level = 0.9)), c(2L, 2L))
})

test_that("a negative variance gives a NaN standard error, not zero", {
  # an expected information wrong in sign for sigma, on a subclass so that
  # the shipped family is left alone
  Bad <- S7::new_class("Gaussian1BadInfo", parent = Gaussian1Distrib)
  S7::method(distrib_expected_hessian, Bad) <- function(distrib, y, theta,
                                                      scale = c("parameter", "link"),
                                                      approx = c("opg", "bartlett",
                                                                 "integrate", "mc"),
                                                      nsim = 10000,
                                                      ...) {
    s <- theta[[2]]
    n <- length(y)
    list(mu_mu = rep(-1 / s^2, n), sigma_sigma = rep(2 / s^2, n),
         mu_sigma = rep(0, n))
  }
  g <- gaussian1_distrib()
  bad <- Bad(distrib_name = g@distrib_name, dimension = g@dimension,
             bounds = g@bounds, params = g@params,
             params_interpretation = g@params_interpretation,
             n_params = g@n_params, params_bounds = g@params_bounds,
             link_params = g@link_params, params_smooth = g@params_smooth)
  set.seed(3)
  f <- suppressWarnings(fit_distrib(bad, rnorm(40), method = optimizers7::bfgs()))
  expect_true(f@vcov["sigma", "sigma"] < 0)
  expect_true(is.nan(f@se[["sigma"]]))
  expect_true(is.nan(f@se_eta[["sigma"]]))
  expect_true(is.finite(f@se[["mu"]]) && f@se[["mu"]] > 0)
})

test_that("log-cdf orders 3 and 4 build on the exact orders 1 and 2", {
  # references: mpmath at 60 digits, differentiating log F directly
  rel <- function(a, b) max(abs(a - b) / abs(b))
  g2 <- gamma2_distrib()
  q <- c(1, 2, 4)
  th <- list(mu = 2, sigma2 = 1)
  d3 <- distrib_deriv3_cdf(g2, q, th)
  d4 <- distrib_deriv4_cdf(g2, q, th)
  expect_lt(rel(d3$mu_mu_mu, c(-1.0424952775183376556, -0.72864495625879279127,
                               -0.12378024546043424565)), 1e-11)
  expect_lt(rel(d3$mu_sigma2_sigma2, c(-3.7021899106247687858,
                                       -0.58268338836363319962,
                                       0.061944240678006997029)), 1e-11)
  expect_lt(rel(d4$mu_mu_mu_mu, c(0.30758840263818972535,
                                  -0.15781790301604106378,
                                  -0.17125356265521436966)), 1e-10)
  st <- student_t1_distrib()
  q <- c(-2, 0.5, 3)
  th <- list(mu = 0.3, sigma = 1.2, nu = 4)
  d3 <- distrib_deriv3_cdf(st, q, th)
  d4 <- distrib_deriv4_cdf(st, q, th)
  expect_lt(rel(d3$mu_nu_nu, c(0.031177213759187651855, 0.0037169514680037025438,
                               -0.00058265203884100452443)), 1e-11)
  expect_lt(rel(d4$mu_mu_nu_nu, c(-0.0030956660194520713851,
                                  0.0024521443003214568315,
                                  0.00062226734704890439522)), 1e-10)
  expect_lt(rel(d4$sigma_nu_nu_nu, c(0.022337617493581422054,
                                     -0.00043814958493376845066,
                                     0.00037467533776178754188)), 1e-10)
})
