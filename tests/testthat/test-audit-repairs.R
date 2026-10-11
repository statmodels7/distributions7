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

# Repairs found by the audit of the continuous families (tranche T2a). Each
# block fails on 0.76.1 and passes on the repaired code. The references marked
# mpmath are 80-digit values of the derivatives of the log-density.

rel_norm <- function(got, ref) max(abs(got - ref)) / max(abs(ref))

test_that("the Gaussian fifth derivatives are closed forms at a small scale", {
  y <- 0.5e-6
  d1 <- distrib_deriv5(gaussian1_distrib(), y, list(mu = 0, sigma = 1e-6))
  expect_equal(d1$sigma_sigma_sigma_sigma_sigma, -24 / 1e-30 + 360 * y^2 / 1e-42,
               tolerance = 1e-13)
  expect_equal(d1$mu_sigma_sigma_sigma_sigma, 120 * y / 1e-36, tolerance = 1e-13)
  v <- 1e-12
  d2 <- distrib_deriv5(gaussian2_distrib(), y, list(mu = 0, sigma2 = v))
  expect_equal(d2$sigma2_sigma2_sigma2_sigma2_sigma2, -12 / v^5 + 60 * y^2 / v^6,
               tolerance = 1e-13)
  tau <- 1e-3
  d3 <- distrib_deriv5(gaussian3_distrib(), y, list(mu = 0, tau = tau))
  expect_equal(d3$tau_tau_tau_tau_tau, 12 / tau^5, tolerance = 1e-13)
  expect_identical(d3$mu_tau_tau_tau_tau, 0)
})

test_that("response derivatives are finite at the boundary with a unit shape", {
  dy <- function(d, y, th) c(distrib_grad_y(d, y, th), distrib_hess_y(d, y, th),
                             distrib_deriv3_y(d, y, th), distrib_deriv4_y(d, y, th))
  # gamma with shape 1 is the exponential: log f = -y/b - log b
  expect_equal(dy(gamma1_distrib(), 0, list(mu = 2, phi = 1)), c(-0.5, 0, 0, 0))
  expect_equal(dy(gamma2_distrib(), 0, list(mu = 3, sigma2 = 9)), c(-1 / 3, 0, 0, 0))
  expect_equal(dy(weibull1_distrib(), 0, list(mu = 2, sigma = 1)), c(-0.5, 0, 0, 0))
  expect_equal(dy(weibull3_distrib(), 0, list(mean = 2, sigma = 1)), c(-0.5, 0, 0, 0))
  # beta(1, 1) is uniform; beta(1, 3) has log f = log 3 + 2 log(1 - y)
  expect_equal(dy(beta1_distrib(), 0, list(mu = 0.5, phi = 2)), c(0, 0, 0, 0))
  expect_equal(dy(beta1_distrib(), 1, list(mu = 0.5, phi = 2)), c(0, 0, 0, 0))
  expect_equal(dy(beta2_distrib(), 0, list(alpha = 1, beta = 3)), c(-2, -2, -4, -12))
  expect_equal(dy(beta2_distrib(), 1, list(alpha = 3, beta = 1)), c(2, -2, 4, -12))
})

test_that("the gamma1 score in mu is exact next to the mean", {
  y <- 3 + 1e-14
  g <- distrib_gradient(gamma1_distrib(), y, list(mu = 3, phi = 0.5))$mu
  expect_equal(g / ((y - 3) / (0.5 * 9)), 1, tolerance = 1e-14)
})

test_that("gamma2 keeps its digits at a large shape", {
  d <- gamma2_distrib()
  # mpmath, mu = 1000, sigma2 = 0.01 (shape 1e8), y = 998
  th <- list(mu = 1000, sigma2 = 0.01)
  expect_lt(rel_norm(unlist(distrib_gradient(d, 998, th)),
                     c(-200.3995341346138, 19976.70673069018)), 1e-12)
  expect_lt(rel_norm(unlist(distrib_hessian(d, 998, th)),
                     c(-100.40040153413462, -4000341.346154703, 20040.053413461712)), 1e-12)
  expect_lt(rel_norm(unlist(distrib_deriv3(d, 998, th)),
                     c(-0.19999999799999998, 10040.04005341346, -4008010.6826923094,
                       1200602403.846411)), 1e-12)
  # mpmath, mu = 1, sigma2 = 1e-6 (shape 1e6), y = 1.0005
  th <- list(mu = 1, sigma2 = 1e-6)
  expect_lt(rel_norm(unlist(distrib_deriv4(d, 1.0005, th)),
                     c(1999993.9999900002, 2000000000002.0002, -1.9980004998333962e+18,
                       -2.9985004998122453e+21, 9.996251505978658e+20)), 1e-9)
  # expected information at shape 1e8: v^-2 (a psi'(a) - 1) a cancels in the
  # direct form
  e <- distrib_expected_hessian(d, 998, list(mu = 1000, sigma2 = 0.01))
  expect_equal(e$sigma2_sigma2, -5000.0000166666665, tolerance = 1e-13)
  expect_equal(e$mu_sigma2, 0.10000000033333333, tolerance = 1e-13)
})

test_that("order-two components follow hess_names() in every route", {
  fams <- list(
    list(gaussian2_distrib(), list(mu = 1, sigma2 = 2), 0.4),
    list(gaussian3_distrib(), list(mu = 1, tau = 2), 0.4),
    list(gamma1_distrib(), list(mu = 1, phi = 0.5), 0.7),
    list(invgauss1_distrib(), list(mu = 1, phi = 0.5), 0.7),
    list(invgauss2_distrib(), list(mu = 1, lambda = 2), 0.7),
    list(weibull1_distrib(), list(mu = 1, sigma = 2), 0.7))
  for (f in fams) {
    d <- f[[1]]; th <- f[[2]]; y <- f[[3]]
    hn <- hess_names(d@params)
    expect_identical(names(distrib_hessian(d, y, th)), hn)
    expect_identical(names(distrib_expected_hessian(d, y, th)), hn)
    expect_identical(names(distrib_hess_cdf(d, y, th)), hn)
  }
})

test_that("lognormal2 orders three to five keep their digits at a small cv", {
  d <- lognormal2_distrib()
  # mpmath, mean = 100, var = 0.09 (cv 0.3 per cent), y = 100.15
  th <- list(mean = 100, var = 0.09)
  expect_lt(rel_norm(unlist(distrib_deriv3(d, 100.15, th)),
                     c(-0.3329957482851836, 122.90192818839209, 410.5979407686597,
                       -344.4766189731582)), 1e-12)
  expect_lt(rel_norm(unlist(distrib_deriv4(d, 100.15, th)),
                     c(0.001107635232198501, 3.7000027746278064, -2731.1539595567274,
                       -13686.598025620879, 68.49292677533268)), 1e-12)
  expect_lt(rel_norm(unlist(distrib_deriv5(d, 100.15, th)),
                     c(7.271839064738202e-08, -0.012308669704404706, -82.22228388895815,
                       91038.46531855753, 608293.2455831502, 504247.4716542791)), 1e-12)
})

test_that("lognormal skewness and kurtosis keep their digits at a small spread", {
  e <- 1e-12
  expect_equal(skewness(lognormal2_distrib(), list(mean = 1, var = e)),
               (3 + e) * sqrt(e), tolerance = 1e-14)
  expect_equal(kurtosis(lognormal2_distrib(), list(mean = 1, var = e)),
               16 * e + 15 * e^2, tolerance = 1e-14)
  w <- expm1(e)
  expect_equal(kurtosis(lognormal1_distrib(), list(mu = 0, sigma2 = e)),
               w * (16 + w * (15 + w * (6 + w))), tolerance = 1e-14)
})

test_that("weibull3 is finite at a shape where its scale underflows", {
  d <- weibull3_distrib()
  s <- 0.002
  th <- list(mean = 1, sigma = s)
  lb <- -lgamma(1 + 1 / s)
  lz <- log(0.5) - lb
  expect_equal(distrib_pdf(d, 0.5, th, log = TRUE),
               log(s) - lb + (s - 1) * lz - exp(s * lz), tolerance = 1e-13)
  expect_equal(distrib_cdf(d, 0.5, th, lower.tail = FALSE, log.p = TRUE),
               -exp(s * lz), tolerance = 1e-13)
  # the median is near exp(-2600), so the quantile is checked in the upper
  # tail on the log scale, where the values are representable
  lp <- distrib_cdf(d, c(0.5, 2), th, lower.tail = FALSE, log.p = TRUE)
  expect_equal(distrib_quantile(d, lp, th, lower.tail = FALSE, log.p = TRUE),
               c(0.5, 2), tolerance = 1e-12)
  expect_true(all(is.finite(distrib_rng(d, 5, th))))
})

test_that("gengamma_components() rejects an order it does not compute", {
  th <- list(a = 2, d = 3, p = 1.5)
  expect_error(gengamma_components(c(1, 2), th, 5), "must be 1, 2, 3 or 4")
  expect_error(gengamma_components(c(1, 2), th, 0), "must be 1, 2, 3 or 4")
})

test_that("pdf, cdf and quantile reject parameters of a non-conformable length", {
  y <- c(1, 2, 3)
  expect_error(distrib_pdf(gengamma1_distrib(), y, list(a = c(1, 2), d = 3, p = 1.5)),
               "dimension mismatch")
  expect_error(distrib_cdf(gamma1_distrib(), y, list(mu = c(1, 2), phi = 1)),
               "dimension mismatch")
  expect_error(distrib_quantile(gaussian1_distrib(), c(0.1, 0.5),
                                list(mu = c(0, 1, 2), sigma = 1)),
               "'p' must have length 1 or 3")
  expect_length(distrib_pdf(gaussian1_distrib(), 0, list(mu = c(0, 1, 2), sigma = 1)), 3L)
  expect_length(distrib_pdf(gaussian1_distrib(), numeric(0), list(mu = 0, sigma = 1)), 0L)
})
