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

test_that("the compiled cdf derivatives keep their digits below a skewed family's center", {
  # gengamma1 at a = 0.3, d = 8, p = 0.6 has its center at a, far below the
  # mass: F(0.45) = 5.3e-10. Reference: mpmath at 50 digits
  d <- gengamma1_distrib()
  th <- list(a = 0.3, d = 8, p = 0.6)
  F0 <- 5.2743263006301626861e-10
  g <- distrib_grad_cdf(d, 0.45, th, log = FALSE)
  expect_equal(unlist(g), c(a = -1.2821114668315061085e-8,
                            d = -2.101558501494046745e-9,
                            p = 3.0620037342448392449e-8), tolerance = 1e-12)
  gl <- distrib_grad_cdf(d, 0.45, th, log = TRUE)
  expect_equal(gl$a, -1.2821114668315061085e-8 / F0, tolerance = 1e-12)
  h <- distrib_hess_cdf(d, 0.45, th, log = FALSE)
  expect_equal(unlist(h), c(a_a = 3.5192871250284835526e-7,
                            d_d = 8.2692148707183725488e-9,
                            p_p = 1.6561318003955667355e-6,
                            a_d = 4.934034206120274724e-8,
                            a_p = -7.4192296758546975624e-7,
                            d_p = -1.1675649433772458047e-7), tolerance = 1e-12)
  d3 <- distrib_deriv3_cdf(d, 0.45, th, log = FALSE)
  expect_equal(c(d3$a_a_a, d3$d_d_d, d3$p_p_p, d3$a_d_p),
               c(-0.000010767830623535425806, -3.2104070943110229934e-8,
                 0.000082709171285004732319, 2.7273099542232438651e-6),
               tolerance = 1e-11)
})

test_that("the compiled cdf derivatives of a small shape reach the mass below the smallest double", {
  # References: mpmath at 40 digits. At shape 0.005 part of the mass lies
  # below 4.9e-324, where the rule cannot place a node
  cdf_all <- function(d, q, th) c(distrib_grad_cdf(d, q, th, log = FALSE),
                                   distrib_hess_cdf(d, q, th, log = FALSE),
                                   distrib_deriv4_cdf(d, q, th, log = FALSE))
  g <- cdf_all(gamma1_distrib(), 1e-100, list(mu = 2, phi = 200))
  expect_equal(unlist(g[c("mu", "phi", "mu_mu", "mu_phi", "phi_phi",
                          "mu_phi_phi_phi", "phi_phi_phi_phi")]),
               c(mu = -0.00076943902331987, phi = 0.001805726832182122,
                 mu_mu = 0.0003866431092182347, mu_phi = -6.671219638559544e-07,
                 phi_phi = -7.424868150084919e-06,
                 mu_phi_phi_phi = -4.157699974890023e-10,
                 phi_phi_phi_phi = 7.497486244573357e-10), tolerance = 1e-10)
  # a quantile below the cut of the rule
  g <- cdf_all(gamma1_distrib(), 1e-300, list(mu = 2, phi = 200))
  expect_equal(unlist(g[c("mu", "phi", "phi_phi", "phi_phi_phi_phi")]),
               c(mu = -7.6943902331987e-05, phi = 0.0005349124482310583,
                 phi_phi = 3.951442524532925e-06,
                 phi_phi_phi_phi = -1.252944720052137e-10), tolerance = 1e-10)
  g <- cdf_all(chisq_distrib(), 1e-200, list(mu = 0.05))
  expect_equal(c(g$mu, g$mu_mu, g$mu_mu_mu_mu),
               c(-0.0022955343752240616, 0.5287417888585819, 28051.33434587288),
               tolerance = 1e-10)
  # gengamma1, where (y/a)^p underflows below the cut
  th <- list(a = 1, d = 0.02, p = 1.5)
  g <- cdf_all(gengamma1_distrib(), 1e-200, th)
  expect_equal(unlist(g[c("a", "d", "p", "a_a", "a_d", "a_p", "a_d_d_p", "p_p_p_p")]),
               c(a = -2.015159028006222e-06, d = -0.04636343772873302,
                 p = -4.975153368261382e-07, a_a = 2.0554622085663467e-06,
                 a_d = 0.0008265108031743492, a_p = 9.950306736522764e-09,
                 a_d_d_p = 0.0012547487652387612, p_p_p_p = 3.3715653580485747e-06),
               tolerance = 1e-10)
  # the upper side, where the density's cutoff at a was crossed by a wide panel
  g <- cdf_all(gengamma1_distrib(), 1e-5, th)
  expect_equal(distrib_deriv3_cdf(gengamma1_distrib(), 1e-5, th, log = FALSE)$a_p_p,
               -0.00010366663713644697, tolerance = 1e-10)
  expect_equal(unlist(g[c("a", "d", "p", "a_a", "a_d", "a_p", "a_d_d_p", "p_p_p_p")]),
               c(a = -0.01600697662786376, d = -8.917963886473382,
                 p = -0.0039519007377655435, a_a = 0.016327115401143465,
                 a_d = -0.6219895595583572, a_p = 7.903226567976277e-05,
                 a_d_d_p = 0.20176308876144255, p_p_p_p = 0.026773836482916757),
               tolerance = 1e-10)
})
# --- T2b: skewnormal, skewt, laplace, logistic, cauchy, gumbel, gpd, chisq,
# exponential, vonmises, pseudohuber, enet (2026-10-11). References are
# closed forms or mpmath at 40 to 80 digits.

test_that("the generalized Pareto distribution and quantile functions keep both tails", {
  d <- gpd_distrib()
  # below zero the lower tail is 0 whatever the sign of the shape
  expect_equal(distrib_cdf(d, c(-10, -5.1, -1), list(sigma = 1.5, xi = 0.3)), c(0, 0, 0))
  expect_equal(distrib_cdf(d, -10, list(sigma = 1.5, xi = 0.3), lower.tail = FALSE), 1)
  expect_equal(distrib_cdf(d, 1e-10, list(sigma = 1, xi = 0)), pexp(1e-10), tolerance = 1e-14)
  expect_equal(distrib_cdf(d, 1e-16, list(sigma = 1, xi = 0.3)), 1e-16, tolerance = 1e-14)
  expect_equal(distrib_cdf(d, 1000, list(sigma = 1, xi = 0), lower.tail = FALSE, log.p = TRUE),
               -1000)
  expect_equal(distrib_quantile(d, 1e-10, list(sigma = 1, xi = 0.3)),
               expm1(-0.3 * log1p(-1e-10)) / 0.3, tolerance = 1e-14)
  expect_equal(distrib_quantile(d, 1e-20, list(sigma = 1, xi = 0), lower.tail = FALSE),
               -log(1e-20), tolerance = 1e-14)
  expect_equal(distributions7:::gpd_endpoint(c(1, 2, 3), -0.5), c(2, 4, 6))
  expect_named(distrib_hess_cdf(d, 0.3, list(sigma = 1, xi = 0.2)),
               c("sigma_sigma", "xi_xi", "sigma_xi"))
})

test_that("the Laplace distribution and quantile functions keep both tails", {
  d <- laplace_distrib()
  th <- list(mu = 0.4, sigma = 1.5)
  expect_equal(distrib_cdf(d, 0.4 + 75, th, lower.tail = FALSE), 0.5 * exp(-50), tolerance = 1e-13)
  expect_equal(distrib_cdf(d, 0.4 - 1500, th, log.p = TRUE), log(0.5) - 1000, tolerance = 1e-14)
  expect_equal(distrib_quantile(d, c(1e-16, 1e-20), th), 0.4 + 1.5 * log(2 * c(1e-16, 1e-20)),
               tolerance = 1e-14)
  expect_equal(distrib_quantile(d, 1e-20, th, lower.tail = FALSE), 0.4 - 1.5 * log(2e-20),
               tolerance = 1e-14)
  expect_equal(distrib_quantile(d, -50, th, log.p = TRUE), 0.4 + 1.5 * (log(2) - 50),
               tolerance = 1e-14)
  expect_equal(distrib_pdf(d, 0, list(mu = 0, sigma = 1e308), log = TRUE),
               -(log(2) + log(1e308)))
  d2 <- laplace2_distrib()
  expect_equal(distrib_cdf(d2, 25.4, list(mu = 0.4, lambda = 2), lower.tail = FALSE),
               0.5 * exp(-50), tolerance = 1e-13)
  # a tiny scale no longer underflows an intermediate power of it
  expect_equal(distrib_hessian(d, 2e-120, list(mu = 0, sigma = 1e-120))$sigma_sigma, -3e240,
               tolerance = 1e-14)
  expect_equal(distrib_deriv4(d, 2e-65, list(mu = 0, sigma = 1e-65))$sigma_sigma_sigma_sigma,
               -4.2e261, tolerance = 1e-14)
})

test_that("seven families compute their fifth derivatives in closed form", {
  # at a fixed standardized residual, sigma^5 times each component does not
  # depend on sigma; the numerical fallback gave a relative error of 5 at 1e-6
  for (fam in c("laplace", "logistic", "cauchy", "gumbel")) {
    d <- get(paste0(fam, "_distrib"))()
    a <- unlist(distrib_deriv5(d, 0.4 + 3e-6, list(mu = 0.4, sigma = 1e-6))) * 1e-30
    b <- unlist(distrib_deriv5(d, 3.4, list(mu = 0.4, sigma = 1)))
    expect_lt(max(abs(a - b)) / max(abs(b)), 1e-10, label = fam)
  }
  expect_equal(distrib_deriv5(laplace2_distrib(), 1, list(mu = 0.4, lambda = 1e-3))[[6]],
               24 / 1e-15, tolerance = 1e-14)
  expect_equal(distrib_deriv5(exponential_distrib(), 3e-4, list(mu = 1e-4))[[1]],
               (120 * 3e-4 - 24e-4) / 1e-24, tolerance = 1e-13)
  expect_equal(distrib_deriv5(chisq_distrib(), 3, list(mu = 1e-4))[[1]],
               -psigamma(5e-5, 4) / 32, tolerance = 1e-14)
})

test_that("the Cauchy derivatives are finite and exact far in the tail", {
  d <- cauchy_distrib()
  th <- list(mu = 0, sigma = 1)
  expect_equal(distrib_pdf(d, 1e160, th, log = TRUE), -log(pi) - 2 * log(1e160))
  h <- distrib_hessian(d, c(1e80, 1e160), th)
  expect_equal(h$sigma_sigma, c(-1, -1))
  expect_equal(distrib_gradient(d, 1e160, th)$sigma, 1)
  # the fourth derivative of F in mu is f''' at z = 1e5, -24 z (z^2 - 1)/(pi (1 + z^2)^4)
  F4 <- distrib_deriv4_cdf(d, 1e5, th, log = FALSE)
  expect_equal(F4$mu_mu_mu_mu, 24 * 1e5 * (1 - 1e10) / (pi * (1 + 1e10)^4), tolerance = 1e-13)
  S3 <- distrib_deriv3_cdf(d, 1e5, th, lower.tail = FALSE)
  expect_true(all(is.finite(unlist(S3))))
  expect_equal(distrib_grad_cdf(d, 1e8, th, lower.tail = FALSE)$mu, 1e-8, tolerance = 1e-7)
})

test_that("the log-scale cdf derivatives read the family's own tails", {
  th <- list(mu = 0, sigma = 1)
  d <- logistic_distrib()
  # d log S/dmu = F for the logistic
  expect_equal(distrib_grad_cdf(d, 40, th, lower.tail = FALSE)$mu, plogis(40))
  expect_equal(distrib_grad_cdf(d, 30, th, lower.tail = FALSE)$mu, plogis(30), tolerance = 1e-14)
  expect_true(all(is.finite(unlist(distrib_deriv4_cdf(d, 300, th, lower.tail = FALSE)))))
  g <- gumbel_distrib()
  # d log F/dmu = -exp(-y) for the Gumbel
  expect_equal(distrib_grad_cdf(g, -10, th)$mu, -exp(10), tolerance = 1e-14)
  expect_equal(distrib_grad_cdf(g, 50, th, lower.tail = FALSE)$mu, 1, tolerance = 1e-14)
  expect_equal(distrib_cdf(g, 800, th, lower.tail = FALSE, log.p = TRUE), -800)
})

test_that("order-two names follow hess_names() in five more places", {
  expect_named(distrib_hessian(vonmises1_distrib(), 0.3, list(mu = 0, kappa = 2)),
               c("mu_mu", "kappa_kappa", "mu_kappa"))
  expect_named(distrib_expected_hessian(vonmises1_distrib(), 0.3, list(mu = 0, kappa = 2)),
               c("mu_mu", "kappa_kappa", "mu_kappa"))
  expect_named(distrib_hess_cdf(skewnormal1_distrib(), 0.3, list(mu = 0, sigma = 1, alpha = 2)),
               hess_names(c("mu", "sigma", "alpha")))
  expect_named(distrib_hess_cdf(skewnormal2_distrib(), 0.3, list(mu = 0, sigma = 1, gamma1 = 0.3)),
               hess_names(c("mu", "sigma", "gamma1")))
})

test_that("small defects: logistic tails, chi-squared at y = 0, overflowing log-densities", {
  d <- logistic_distrib()
  expect_equal(distrib_hess_y(d, 40, list(mu = 0, sigma = 1))[[1]],
               -2 * exp(-40) / (1 + exp(-40))^2, tolerance = 1e-14)
  expect_equal(distrib_hessian(d, 40, list(mu = 0, sigma = 1))$mu_mu,
               -2 * exp(-40) / (1 + exp(-40))^2, tolerance = 1e-14)
  expect_equal(distrib_grad_y(chisq_distrib(), 0, list(mu = 2))[[1]], -0.5)
  expect_equal(distrib_hess_y(chisq_distrib(), 0, list(mu = 2))[[1]], 0)
  expect_true(is.finite(distrib_pdf(pseudohuber_distrib(), 1e160,
                                    list(mu = 0.4, sigma = 1.2, nu = 2), log = TRUE)))
  expect_true(is.finite(distrib_pdf(pseudohuber2_distrib(), 1e300,
                                    list(mu = 0, sigma = 1, nu = 2), log = TRUE)))
  st <- skewt_distrib()
  th <- list(mu = 0, sigma = 1, alpha = 3, nu = 6)
  expect_equal(distrib_hess_y(st, 1e90, th)[[1]], distrib_hessian(st, 1e90, th)$mu_mu)
})

test_that("the skew normal and von Mises draws cover their whole parameter domain", {
  set.seed(1)
  x <- distrib_rng(skewnormal1_distrib(), 5000, list(mu = 0, sigma = 1, alpha = 1e200))
  expect_true(all(x >= 0))
  set.seed(1)
  expect_false(anyNA(distrib_rng(vonmises1_distrib(), 50, list(mu = 0.5, kappa = 1e-8))))
  set.seed(1)
  expect_false(anyNA(distrib_rng(vonmises1_distrib(), 3, list(mu = 0.5, kappa = c(1e-9, 1, 1e-310)))))
  set.seed(1)
  expect_false(anyNA(distrib_rng(vonmises2_distrib(), 5, list(mu = 0.5, rho = 1e-9))))
})

test_that("the elastic net's derivatives keep their digits as alpha approaches one", {
  d <- enet_distrib()
  # mpmath, 60 digits
  expect_equal(distrib_gradient(d, 0.4, list(mu = 0, lambda = 1, alpha = 1 - 1e-6))$alpha,
               -0.3199980000149999, tolerance = 1e-12)
  expect_equal(distrib_hessian(d, 3, list(mu = 0, lambda = 3, alpha = 0.999))$alpha_alpha,
               -0.22214866375342673, tolerance = 1e-12)
  expect_equal(variance(d, list(mu = 0, lambda = 2, alpha = 0.9999)), 0.49997501124519006,
               tolerance = 1e-13)
  th <- list(mu = 0, lambda = 2, alpha = 0.5)
  expect_equal(distrib_cdf(d, 12, th, lower.tail = FALSE), distrib_cdf(d, -12, th), tolerance = 1e-14)
  expect_equal(distrib_quantile(d, 1e-30, th, lower.tail = FALSE), -distrib_quantile(d, 1e-30, th),
               tolerance = 1e-14)
  expect_true(is.finite(distrib_cdf(d, -40, th, log.p = TRUE)))
})

test_that("the skew normal's kernels, distribution function and centered information", {
  d <- skewnormal1_distrib()
  # mpmath, 80 digits: d^4/dalpha^4 log f at y = 1, alpha = -38, and at y = -30, alpha = 5
  expect_equal(distrib_deriv4(d, 1, list(mu = 0, sigma = 1, alpha = -38))$alpha_alpha_alpha_alpha,
               2.8381325760967472e-06, tolerance = 1e-12)
  expect_equal(distrib_deriv4(d, -30, list(mu = 0, sigma = 1, alpha = 5))$alpha_alpha_alpha_alpha,
               0.009591473298469228, tolerance = 1e-12)
  th <- list(mu = 0, sigma = 1, alpha = 3)
  # mpmath, 40 digits
  expect_equal(distrib_cdf(d, c(-2, -3), th), c(5.0891259751793012e-12, 3.2642743096552352e-23),
               tolerance = 1e-12)
  expect_equal(distrib_cdf(d, 8, th, lower.tail = FALSE), 1.2441921148543568e-15, tolerance = 1e-12)
  expect_equal(distrib_cdf(d, 10, th, lower.tail = FALSE, log.p = TRUE), -52.538137969952525,
               tolerance = 1e-13)
  expect_equal(distrib_cdf(d, c(-0.1, 1), list(mu = 0, sigma = 1, alpha = 1e200)),
               c(0, 2 * pnorm(1) - 1), tolerance = 1e-14)
  # near the skewness ceiling the centered information is negative definite
  e <- distrib_expected_hessian(skewnormal2_distrib(), 0, list(mu = 0, sigma = 1, gamma1 = 0.99527))
  expect_lt(e$gamma1_gamma1, 0)
})

test_that("tails below double precision resolution: compiled cdf, quantiles, von Mises", {
  expect_true(is.finite(distrib_cdf(pseudohuber_distrib(), -2000, list(mu = 0.4, sigma = 1.2, nu = 2),
                                    log.p = TRUE)))
  expect_true(is.finite(distrib_cdf(skewt_distrib(), -1e50,
                                    list(mu = 0.3, sigma = 1.2, alpha = 2, nu = 6), log.p = TRUE)))
  d <- skewt_distrib()
  th <- list(mu = 0.3, sigma = 1.2, alpha = 2, nu = 6)
  q <- distrib_quantile(d, 1e-100, th)
  expect_equal(distrib_cdf(d, q, th, log.p = TRUE), log(1e-100), tolerance = 1e-12)
  expect_equal(distrib_quantile(d, distrib_cdf(d, c(-2, 0.5, 3), th), th), c(-2, 0.5, 3),
               tolerance = 1e-10)
  ph <- pseudohuber_distrib()
  thp <- list(mu = 0.4, sigma = 1.2, nu = 2)
  expect_equal(distrib_quantile(ph, 1e-30, thp, lower.tail = FALSE),
               0.8 - distrib_quantile(ph, 1e-30, thp), tolerance = 1e-12)
  vm <- vonmises1_distrib()
  # mpmath, 40 digits
  expect_equal(distrib_cdf(vm, -1, list(mu = 0, kappa = 100)), 5.1016782326113297e-22,
               tolerance = 1e-11)
  expect_equal(distrib_cdf(vm, -1.5, list(mu = 0, kappa = 100), log.p = TRUE),
               -96.147161269784047, tolerance = 1e-12)
})

test_that("the pseudo-Huber derivatives in nu at both ends of its domain", {
  # mpmath, 60 to 80 digits
  expect_equal(distrib_hessian(pseudohuber_distrib(), 3, list(mu = 0, sigma = 1, nu = 1e-14))$nu_nu,
               -24999999999934.105, tolerance = 1e-12)
  th <- list(mu = 2, sigma = 2, nu = 1e8)
  expect_equal(distrib_hessian(pseudohuber2_distrib(), 3, th)$nu_nu, 1.4644199203098685e-21,
               tolerance = 1e-10)
  expect_equal(distrib_deriv3(pseudohuber2_distrib(), 3, th)$nu_nu_nu, -3.660837885159833e-29,
               tolerance = 1e-10)
})

test_that("the skew t's derivatives in nu at large nu", {
  d <- skewt_distrib()
  # mpmath, 80 digits
  h <- distrib_hessian(d, 200, list(mu = 0, sigma = 1, alpha = -5, nu = 1000))
  expect_equal(h$nu_nu, 0.0004990377703385572, tolerance = 1e-8)
  expect_equal(h$alpha_nu, 0.19193841843404105, tolerance = 1e-8)
})

test_that("the skew t's components free of nu at large nu, the second table", {
  # mpmath, 60 digits; the reduced table lost 7e-5 in sigma^4 at nu = 1e6
  d4 <- distrib_deriv4(skewt_distrib(), 1, list(mu = 0, sigma = 1, alpha = 3, nu = 1e6))
  expect_equal(d4$sigma_sigma_sigma_sigma, -53.00466070727036, tolerance = 1e-13)
  expect_equal(d4$mu_mu_sigma_sigma, -7.504473645143255, tolerance = 1e-13)
})

test_that("the centered skew normal's kernels far below the step", {
  # mpmath, 50 digits, at y = -9, gamma1 = 0.8 (4e-7 normwise at order 5 before)
  d <- skewnormal2_distrib()
  th <- list(mu = 0, sigma = 1, gamma1 = 0.8)
  d5 <- distrib_deriv5(d, -9, th)
  expect_equal(d5$gamma1_gamma1_gamma1_gamma1_gamma1, -95247826.17814274, tolerance = 1e-12)
  expect_equal(d5$mu_mu_mu_mu_mu, -0.0007906296480186424, tolerance = 1e-9)
  expect_equal(distrib_deriv3(d, -9, th)$gamma1_gamma1_gamma1, -181653.25348333703,
               tolerance = 1e-12)
})

test_that("the pseudo-Huber nu-only derivatives in the middle of nu's range", {
  # mpmath, 60 digits, at nu = 200 (1e-10 at order five by the Bessel ratios)
  o <- distributions7:::ph_nu_derivs_cpp(200, 5L)
  expect_equal(o$B[5, 1], -1.2840692384019815e-10, tolerance = 1e-13)
  expect_equal(o$logR[5, 1], 2.9029818456005618e-11, tolerance = 1e-12)
})

test_that("a probability outside the unit interval gives NaN in every T2b quantile", {
  fams <- list(
    list(laplace_distrib(), list(mu = 0, sigma = 1)),
    list(laplace2_distrib(), list(mu = 0, lambda = 1)),
    list(gpd_distrib(), list(sigma = 1, xi = 0.3)),
    list(gpd_distrib(), list(sigma = 1, xi = 0)),
    list(enet_distrib(), list(mu = 0, lambda = 2, alpha = 0.5)),
    list(pseudohuber_distrib(), list(mu = 0, sigma = 1, nu = 2)),
    list(pseudohuber2_distrib(), list(mu = 0, sigma = 1, nu = 2)))
  for (f in fams) {
    d <- f[[1]]
    th <- f[[2]]
    q <- suppressWarnings(c(
      distrib_quantile(d, c(-0.1, 1.2), th),
      distrib_quantile(d, c(-0.1, 1.2), th, lower.tail = FALSE),
      distrib_quantile(d, 0.1, th, log.p = TRUE),
      distrib_quantile(d, 0.1, th, log.p = TRUE, lower.tail = FALSE)))
    expect_true(all(is.nan(q)), label = d@distrib_name)
  }
})

test_that("the von Mises tails at a large concentration and across the antimode", {
  # mpmath, 40 digits (vm_tail_ref.py); +Inf and -5329.67 before
  d <- vonmises1_distrib()
  expect_equal(distrib_cdf(d, -0.1, list(mu = 0, kappa = 1e6), log.p = TRUE),
               -5001.357263403596, tolerance = 1e-12)
  # the upper tail [-2.5, pi) holds both ends of the antimode's arc
  expect_equal(distrib_cdf(d, -2.5, list(mu = -3, kappa = 1e6), lower.tail = FALSE,
                           log.p = TRUE),
               -10013.371998567147, tolerance = 1e-12)
})

test_that("the von Mises tail between 1e-10 and 1e-6 comes from the quadrature", {
  # mpmath, 40 digits (vm_band_ref.py); the series left up to 3e-5 here
  d <- vonmises1_distrib()
  expect_equal(distrib_cdf(d, 2.895, list(mu = 2.9, kappa = 1e6), log.p = TRUE),
               -15.064968139398426, tolerance = 1e-10)
})

test_that("the logistic expected third and fourth derivatives are closed form", {
  # mpmath, 40 digits (logis_exp_ref.py), at sigma = 1 and scaled by sigma^-k
  d <- logistic_distrib()
  s <- c(0.7, 1.4)
  th <- list(mu = c(0.5, -2), sigma = s)
  e3 <- distrib_deriv3(d, c(0, 0), th, expected = TRUE)
  e4 <- distrib_deriv4(d, c(0, 0), th, expected = TRUE)
  expect_identical(e3$mu_mu_mu, c(0, 0))
  expect_equal(e3$sigma_sigma_sigma, 5.9348022005446793094 / s^3, tolerance = 1e-14)
  expect_equal(e4$mu_mu_sigma_sigma, -1.1140087910869031418 / s^4, tolerance = 1e-14)
  expect_equal(e4$sigma_sigma_sigma_sigma, -26.888440015180436047 / s^4,
               tolerance = 1e-14)
  # and the numerical route agrees
  num <- distributions7:::expected_derivative(d, 0, list(mu = 0.5, sigma = 0.7),
                                              order = 4L, approx = "integrate")
  expect_equal(unlist(lapply(e4, `[`, 1)), unlist(num)[names(e4)], tolerance = 1e-8)
})

# relative error against a reference, whatever the reference's magnitude:
# expect_equal() compares absolutely where the expected value is below the
# tolerance, which passes any value of order 1e-28 against any other
expect_rel <- function(object, expected, tolerance) {
  expect_lt(max(abs(object - expected) / abs(expected)), tolerance)
}

test_that("the binomial recycles a vector size in every derivative", {
  d <- binomial_distrib(size = c(5, 10, 20))
  th <- list(mu = 0.3)
  expect_rel(distrib_gradient(d, 2, th)$mu, (2 - c(5, 10, 20) * 0.3) / 0.21,
               tolerance = 1e-14)
  expect_rel(distrib_expected_hessian(d, 2, th)$mu_mu, -c(5, 10, 20) / 0.21,
               tolerance = 1e-14)
  expect_length(distrib_deriv4(d, 2, th)$mu_mu_mu_mu, 3L)
  # a size of length 2 against three responses was read out of bounds
  expect_error(distrib_gradient(binomial_distrib(size = c(5, 10)), 1:3, th),
               "'size' has length 2")
})

test_that("the binomial size is validated", {
  for (s in list(2.5, 0, -3, c(5, NA), Inf)) {
    expect_error(binomial_distrib(size = s), "positive integers")
  }
  expect_no_error(binomial_distrib(size = c(3, 7)))
})

test_that("the count families' derivatives at a tiny mean and a zero count", {
  expect_equal(distrib_hessian(geometric_distrib(), 0, list(mu = 1e-170))$mu_mu, 1)
  expect_equal(distrib_deriv3(geometric_distrib(), 0, list(mu = 1e-160))[[1]], -2)
  expect_equal(distrib_hessian(bernoulli_distrib(), 0, list(mu = 1e-170))$mu_mu, -1)
  expect_equal(distrib_deriv3(bernoulli_distrib(), 0, list(mu = 1e-110))[[1]], -2)
  expect_equal(distrib_hessian(binomial_distrib(size = 10), 0,
                               list(mu = 1e-170))$mu_mu, -10)
  expect_identical(distrib_hessian(poisson_distrib(), 0, list(mu = 1e-170))$mu_mu, 0)
})

test_that("the geometric derivatives at a count near a large mean", {
  m <- 1e14
  q <- m * (1 + m)
  d <- geometric_distrib()
  expect_rel(distrib_hessian(d, m, list(mu = m))$mu_mu, -1 / q, tolerance = 1e-14)
  expect_rel(distrib_deriv3(d, m, list(mu = m))[[1]], 2 * (1 + 2 * m) / q^2,
               tolerance = 1e-14)
  expect_rel(distrib_deriv3(d, 0, list(mu = m), expected = TRUE)[[1]],
               2 * (1 + 2 * m) / q^2, tolerance = 1e-14)
})

test_that("a missing count has a missing mass in every discrete family", {
  fams <- list(list(betabinom1_distrib(size = 10), list(mu = 0.3, sigma = 0.2)),
               list(betabinom2_distrib(size = 10), list(alpha = 2, beta = 3)),
               list(pig1_distrib(), list(mu = 3, sigma = 0.8)),
               list(pig2_distrib(), list(mu = 3, alpha = 1)))
  for (f in fams) {
    p <- distrib_pdf(f[[1]], c(NA, 1), f[[2]])
    expect_true(is.na(p[1]), label = f[[1]]@distrib_name)
    expect_true(p[2] > 0)
  }
})

test_that("the betabinom1 quantile gives NaN outside the unit interval", {
  d <- betabinom1_distrib(size = 10)
  th <- list(mu = 0.3, sigma = 0.2)
  expect_identical(distrib_quantile(d, c(-0.5, 0, 1, 1.5), th), c(NaN, 0, 10, NaN))
  expect_identical(distrib_quantile(d, c(-1, 2), th, lower.tail = FALSE), c(NaN, NaN))
  expect_identical(distrib_quantile(d, 0.5, th, log.p = TRUE), NaN)
})

test_that("the upper tails of the PIG and the beta-binomials are summed", {
  # mpmath, 40 digits, the sum of the mass over 201:1199
  d <- pig1_distrib()
  th <- list(mu = 3, sigma = 0.8)
  expect_rel(distrib_cdf(d, 200, th, lower.tail = FALSE), 1.7461679759672038e-19,
               tolerance = 1e-12)
  expect_rel(distrib_cdf(d, 200, th, lower.tail = FALSE, log.p = TRUE),
               log(1.7461679759672038e-19), tolerance = 1e-14)
  q <- distrib_quantile(d, 1e-17, th, lower.tail = FALSE)
  expect_true(is.finite(q))
  expect_lte(distrib_cdf(d, q, th, lower.tail = FALSE), 1e-17)
  expect_gt(distrib_cdf(d, q - 1, th, lower.tail = FALSE), 1e-17)
  d <- pig2_distrib()
  th <- list(mu = 3, alpha = 10)
  expect_rel(distrib_cdf(d, 60, th, lower.tail = FALSE),
               sum(distrib_pdf(d, 61:5000, th)), tolerance = 1e-12)
  for (d in list(betabinom1_distrib(size = 200), betabinom2_distrib(size = 200))) {
    th <- if (d@params[1] == "mu") list(mu = 0.01, sigma = 1e-4) else list(alpha = 100, beta = 9900)
    up <- sum(distrib_pdf(d, 151:200, th))
    expect_rel(distrib_cdf(d, 150, th, lower.tail = FALSE), up, tolerance = 1e-12)
    q <- distrib_quantile(d, 1e-20, th, lower.tail = FALSE)
    expect_lt(q, 150)
    expect_lte(distrib_cdf(d, q, th, lower.tail = FALSE), 1e-20)
  }
})

test_that("the NB1 mass and derivatives at the Poisson boundary", {
  # mpmath, 80 digits (vd_T3_ref.py)
  d <- negbin1_distrib()
  expect_rel(distrib_pdf(d, 1, list(mu = 500, theta = 1e-17), log = TRUE),
               -493.78539190157780577, tolerance = 1e-14)
  expect_rel(distrib_pdf(d, 0, list(mu = 500, theta = 1e-15), log = TRUE),
               -499.99999999999975, tolerance = 1e-14)
  th <- list(mu = 4, theta = 1e-8)
  expect_rel(distrib_hessian(d, 3, th)$theta_theta, 0.020833336145833253411,
               tolerance = 1e-12)
  expect_rel(distrib_hessian(d, 0, th)$theta_theta, -2.6666666066666676267,
               tolerance = 1e-13)
  expect_rel(distrib_gradient(d, 3, th)$theta, -0.2499999997916666526,
               tolerance = 1e-14)
  expect_rel(distrib_hessian(d, 3, th)$mu_theta, 0.31249999489583339779,
               tolerance = 1e-14)
})

test_that("the NB1 distribution function, quantile and draws at a tiny dispersion", {
  d <- negbin1_distrib()
  th <- list(mu = 3, theta = 1e-16)
  expect_rel(distrib_cdf(d, 3, th), 0.64723188878223125873, tolerance = 1e-13)
  expect_identical(distrib_quantile(d, 0.5, th), 3)
  set.seed(1)
  expect_rel(mean(distrib_rng(d, 1e4, list(mu = 3, theta = 1e-17))), 3, tolerance = 0.05)
})

test_that("the NB2 derivatives in the mean at a small size", {
  # mpmath, 80 digits
  d <- negbin2_distrib()
  th <- list(mu = 5, theta = 1e-8)
  expect_rel(distrib_hessian(d, 3, th)$mu_mu, -8.000000015999999904e-11, tolerance = 1e-14)
  expect_rel(distrib_deriv3(d, 3, th)$mu_mu_mu, 1.27999999808e-10, tolerance = 1e-14)
  expect_rel(distrib_deriv4(d, 3, th)$mu_mu_mu_mu, -1.3439999961600000077e-10,
               tolerance = 1e-14)
})

test_that("the beta-binomial derivatives at a large concentration", {
  # mpmath, 80 digits (vd_T3_ref.py, vd_T3_bbref.py); the betabinom2 value is
  # the closed form in mpmath's polygammas at 60 digits
  d <- betabinom1_distrib(size = 10)
  th <- list(mu = 0.3, sigma = 1e-4)
  expect_rel(distrib_hessian(d, 3, th)$sigma_sigma, 43.649080886004494728, tolerance = 1e-13)
  expect_rel(distrib_deriv3(d, 3, th)$sigma_sigma_sigma, -809.65350400224666619,
               tolerance = 1e-13)
  expect_rel(distrib_deriv4(d, 3, th)$sigma_sigma_sigma_sigma, 22471.126081180217613,
               tolerance = 1e-13)
  expect_rel(distrib_gradient(d, 3, list(mu = 0.4, sigma = 1e-15))$sigma,
               -2.4999999999999990278, tolerance = 1e-14)
  expect_rel(distrib_expected_hessian(d, 3, list(mu = 0.3, sigma = 1e-8))$sigma_sigma,
               -44.999990657144356683, tolerance = 1e-12)
  d100 <- betabinom1_distrib(size = 100)
  expect_rel(distrib_deriv4(d100, 37, list(mu = 0.4, sigma = 1e-3))$sigma_sigma_sigma_sigma,
               -83637651.88739765558899421, tolerance = 1e-12)
  d2 <- betabinom2_distrib(size = 10)
  expect_rel(distrib_deriv4(d2, 3, list(alpha = 4e11, beta = 6e11))$alpha_alpha_alpha_alpha,
               -6.4312499999404875e-46, tolerance = 1e-13)
})

test_that("the betabinom1 expected-information derivatives at a large concentration", {
  # mpmath, 60 digits (vd_T3_bbexp.py)
  d <- betabinom1_distrib(size = 10)
  th <- list(mu = 0.3, sigma = 1e-6)
  e1 <- distrib_dexpected_hessian(d, 1, th)
  e2 <- distrib_d2expected_hessian(d, 1, th)
  expect_rel(e1$sigma_sigma_sigma, 934.25572412801040184, tolerance = 1e-10)
  expect_rel(e1$mu_mu_mu, 90.702131529692123424, tolerance = 1e-10)
  expect_rel(e2$sigma_sigma_sigma_sigma, -29989.499093604765334, tolerance = 1e-9)
})

test_that("the betabinom1 expected higher derivatives are exact sums", {
  d <- betabinom1_distrib(size = 10)
  th <- list(mu = 0.3, sigma = 0.2)
  w <- distrib_pdf(d, 0:10, th)
  o3 <- distrib_deriv3(d, 0:10, th)
  e3 <- distrib_deriv3(d, 0, th, expected = TRUE, approx = "mc", nsim = 10)
  for (k in names(o3)) expect_rel(e3[[k]], sum(w * o3[[k]]), tolerance = 1e-13)
})

test_that("the beta-binomial moments at a small dispersion and a large size", {
  # mpmath, 120 digits (vd_T3_ref.py)
  d <- betabinom1_distrib(size = 1e6)
  th <- list(mu = 0.3, sigma = 1e-6)
  expect_rel(kurtosis(d, th), -2.0476192380843810055e-6, tolerance = 1e-13)
  expect_rel(skewness(betabinom1_distrib(size = 1e8), list(mu = 0.3, sigma = 1e-8)),
               0.00018516401717705004623, tolerance = 1e-13)
  d2 <- betabinom2_distrib(size = 1e6)
  expect_rel(kurtosis(d2, list(alpha = 3e5, beta = 7e5)), -2.0476192380843810055e-6,
               tolerance = 1e-12)
})

test_that("the skew normal and skew t moments at a huge alpha and a huge nu", {
  sn <- skewnormal1_distrib()
  expect_rel(mean(sn, list(mu = 0, sigma = 1, alpha = c(1e160, -1e200))),
               c(1, -1) * sqrt(2 / pi), tolerance = 1e-15)
  # mpmath, 80 digits
  expect_rel(variance(skewt_distrib(), list(mu = 0, sigma = 1, alpha = 3, nu = 1e14)),
               0.42704220486918819687, tolerance = 1e-14)
  expect_equal(mean(skewt_distrib(), list(mu = 0, sigma = 1, alpha = 1e160, nu = 5)),
               mean(skewt_distrib(), list(mu = 0, sigma = 1, alpha = 1e100, nu = 5)))
})

test_that("the skew t's divergent moments are infinite, as the Student t's", {
  d <- skewt_distrib()
  expect_identical(variance(d, list(mu = 0, sigma = 1, alpha = 2, nu = c(0.5, 1.5, 2))),
                   c(NaN, Inf, Inf))
  expect_identical(kurtosis(d, list(mu = 0, sigma = 1, alpha = 2, nu = c(1.5, 3, 4))),
                   c(NaN, Inf, Inf))
})

test_that("the Weibull moments at a large and a small shape", {
  # mpmath, 150 digits
  d <- weibull1_distrib()
  expect_rel(kurtosis(d, list(mu = 1, sigma = 1e5)), 2.3997108383775288848, tolerance = 1e-13)
  expect_rel(skewness(d, list(mu = 1, sigma = 1e5)), -1.1394874345084644293, tolerance = 1e-13)
  expect_rel(kurtosis(d, list(mu = 1, sigma = 0.02)), 9.0548514656103281165e+58,
               tolerance = 1e-12)
})

test_that("the beta moments at extreme precisions", {
  expect_rel(skewness(beta1_distrib(), list(mu = 0.3, phi = 1e-170)),
               0.87287156094396953, tolerance = 1e-15)
  expect_rel(kurtosis(beta1_distrib(), list(mu = 0.3, phi = 1e-170)),
               -26 / 21, tolerance = 1e-15)
  expect_equal(variance(beta2_distrib(), list(alpha = 1e-300, beta = 1e-300)), 0.25)
  expect_rel(variance(beta2_distrib(), list(alpha = 1e200, beta = 1e200)), 1.25e-201,
               tolerance = 1e-14)
  expect_rel(skewness(beta2_distrib(), list(alpha = 1e-200, beta = 2e-200)),
               1 / sqrt(2), tolerance = 1e-15)
})

test_that("the NB1 expected information and its derivatives at a small dispersion", {
  # mpmath, 60 digits: G = E[l_theta,theta] = -E[D] - mu q'(theta) and its
  # derivatives by mp.diff; E[l_mu,theta] = -theta G/mu
  d <- negbin1_distrib()
  ref <- list(c(4, 1e-12, -0.499999999998875, 1.1249999999961042, -3.8958333333149113),
              c(4, 1e-8, -0.4999999887500002, 1.1249999610416677, -3.895833149114589),
              c(50, 1e-8, -0.49999998990000016, 1.009999969329334, -3.06706654239443),
              c(4, 0.6, -0.18060194723592765, 0.23761269355675782, -0.47449809149075506))
  for (r in ref) {
    th <- list(mu = r[1], theta = r[2])
    e <- distrib_expected_hessian(d, 0, th)
    expect_rel(e$theta_theta, r[3], tolerance = 1e-12)
    expect_rel(e$mu_theta, -r[2] * r[3] / r[1], tolerance = 1e-12)
    expect_rel(distrib_dexpected_hessian(d, 0, th)$theta_theta_theta, r[4], tolerance = 1e-12)
    expect_rel(distrib_d2expected_hessian(d, 0, th)$theta_theta_theta_theta, r[5],
               tolerance = 1e-12)
  }
})

test_that("the NB1 observed curvature in theta below theta = 1e-9", {
  # mpmath, 40 digits: (y - mu)/(1 + theta)^2 - D_y - mu q'(theta)
  d <- negbin1_distrib()
  expect_rel(distrib_hessian(d, 3, list(mu = 4, theta = 1e-10))$theta_theta,
             0.020833333361458333323, tolerance = 1e-13)
})

test_that("the cdf derivatives read a vector size per observation", {
  d3 <- binomial_distrib(size = c(5, 10, 20))
  th <- list(mu = 0.3)
  # the closed form of the gradient of F, -n dbinom(k, n - 1, p), per size
  expect_rel(distrib_grad_cdf(d3, c(1, 4, 9), th, log = FALSE)$mu,
             -c(5, 10, 20) * dbinom(c(1, 4, 9), c(4, 9, 19), 0.3), tolerance = 1e-14)
  for (f in list(distrib_hess_cdf, distrib_deriv3_cdf, distrib_deriv4_cdf)) {
    v <- unlist(f(d3, 2, th))
    s <- vapply(c(5, 10, 20), function(n) unlist(f(binomial_distrib(size = n), 2, th)), 0)
    expect_rel(v, s, tolerance = 1e-12)
  }
  z3 <- zero_inflated(binomial_distrib(size = c(5, 10, 20)))
  thz <- list(mu = 0.3, zi = 0.2)
  v <- distrib_hess_cdf(z3, 2, thz)
  s <- lapply(c(5, 10, 20), function(n) distrib_hess_cdf(zero_inflated(binomial_distrib(size = n)), 2, thz))
  for (k in names(v)) expect_rel(v[[k]], vapply(s, `[[`, 0, k), tolerance = 1e-12)
  # a size whose length is neither 1 nor that of the response is rejected
  expect_error(distrib_grad_cdf(binomial_distrib(size = c(5, 10)), 1:3, th), "'size' has length 2")
  # a missing quantile gives a missing derivative
  expect_true(all(is.na(unlist(lapply(distrib_hess_cdf(negbin2_distrib(), NA, list(mu = 3, theta = 2)),
                                       identity)))))
})

# T4 (2026-10-12): the multivariate families. References are mpmath at 60 to
# 400 digits (vd_T4_testref.py, vd_T4_4.py, vd_T4_dirref2.py).

test_that("the higher orders of a t on the precision read the precision", {
  # references: mpmath at 60 digits on the log-density (vd_T4_testref.py)
  d <- mvstudent_t2_distrib(2)
  th <- list(mu1 = 0.5, mu2 = -0.3, omega_log_L1 = 0.1, omega_log_L2 = -0.2,
             omega_L2.1 = 0.4, nu = 6)
  y <- rbind(c(0.7, -0.1))
  d3 <- distrib_deriv3(d, y, th)
  expect_rel(c(d3$omega_log_L1_omega_log_L1_omega_log_L1, d3$mu1_mu1_omega_log_L1),
             c(-0.26892228239176166, -2.9977214752055362), tolerance = 1e-12)
  d4 <- distrib_deriv4(d, y, th)
  expect_rel(c(d4$omega_log_L1_omega_log_L1_omega_log_L1_omega_log_L1,
               d4$mu1_omega_log_L1_omega_log_L1_nu),
             c(-0.49377410560316784, -0.043808656253851526), tolerance = 1e-12)
})

test_that("the derivatives at a nearly singular covariance read the chart's inverse", {
  y <- rbind(c(0.7, -0.1))
  g <- mvgaussian1_distrib(2)
  th <- list(mu1 = 0.5, mu2 = -0.3, sigma_log_L1 = 0.1, sigma_log_L2 = -11.5,
             sigma_L2.1 = 0.4)
  expect_rel(c(distrib_hessian(g, y, th)$sigma_L2.1_sigma_L2.1,
               distrib_deriv3(g, y, th)$sigma_log_L1_sigma_L2.1_sigma_L2.1,
               distrib_deriv4(g, y, th)$sigma_log_L1_sigma_log_L1_sigma_L2.1_sigma_L2.1),
             c(-319134810.56577108, 638269621.13154215, -1276539242.2630843),
             tolerance = 1e-12)
  th9 <- list(mu1 = 0, mu2 = 0, sigma_log_L1 = 0, sigma_log_L2 = -9, sigma_L2.1 = 1)
  expect_rel(distrib_hessian(g, rbind(c(0.3, -0.8)), th9)$sigma_log_L1_sigma_log_L1,
             -27577187.217678815, tolerance = 1e-12)
  t1 <- mvstudent_t1_distrib(2)
  tt <- c(th, list(nu = 6))
  expect_rel(c(distrib_hessian(t1, y, tt)$sigma_L2.1_sigma_L2.1,
               distrib_deriv3(t1, y, tt)$sigma_L2.1_sigma_L2.1_sigma_L2.1),
             c(16.087966413079724, 45.628553775988601), tolerance = 1e-10)
})

test_that("a nearly singular covariance is sampled through the chart's factor", {
  th <- list(mu1 = 0, mu2 = 0, sigma_log_L1 = 0, sigma_log_L2 = -20, sigma_L2.1 = 1)
  set.seed(1)
  y <- distrib_rng(mvgaussian1_distrib(2), 50, th)
  # L = [[1, 0], [1, e^-20]], so the two coordinates differ by e^-20 z
  expect_lt(max(abs(y[, 2] - y[, 1])), 1e-7)
  set.seed(1)
  yt <- distrib_rng(mvstudent_t1_distrib(2), 50, c(th, list(nu = 5)))
  expect_true(all(is.finite(yt)))
  # the precision side: the sample covariance approaches the inverse
  om <- list(mu1 = 0, mu2 = 0, omega_log_L1 = 0.3, omega_log_L2 = -0.4, omega_L2.1 = 0.5)
  d2 <- mvgaussian2_distrib(2, parameters7::log_cholesky(2))
  set.seed(2)
  expect_lt(max(abs(var(distrib_rng(d2, 2e5, om)) - mv_sigma(d2, om))), 0.02)
})

test_that("a response with no rows gives empty results without a warning", {
  th <- list(mu1 = 0.5, mu2 = -0.3, sigma_log_L1 = 0.1, sigma_log_L2 = -0.2,
             sigma_L2.1 = 0.4)
  y0 <- matrix(numeric(0), 0, 2)
  expect_silent(distrib_cross_y(mvgaussian1_distrib(2), y0, th))
  expect_silent(distrib_grad_y_hess(mvgaussian1_distrib(2), y0, th))
  t1 <- mvstudent_t1_distrib(2)
  tt <- c(th, list(nu = 6))
  expect_silent(distrib_cross_y(t1, y0, tt))
  expect_silent(distrib_hess_y_hess(t1, y0, tt))
})

test_that("the arguments of the multivariate families are validated", {
  d <- mvgaussian1_distrib(3)
  th <- list(mu1 = 1, mu2 = 0, mu3 = 0, sigma_log_L1 = 0, sigma_log_L2 = 0,
             sigma_log_L3 = 0, sigma_L2.1 = 0, sigma_L3.1 = 0, sigma_L3.2 = 0)
  expect_error(as_mv_matrix(mvgaussian1_distrib(2), c("a", "b")), "must be numeric")
  expect_error(mv_marginal(d, th, 1.7), "distinct coordinates")
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  expect_error(plot(d, th, which = integer(0)), "distinct coordinates")
  expect_error(multinomial_distrib(3.5, size = 4), "'n_dim' must be")
  expect_error(multinomial_distrib(3.5, size = 4, probs = parameters7::simplex(3)),
               "'n_dim' must be")
  expect_error(dirichlet_distrib(3.7), "'n_dim' must be")
  expect_error(dirichlet_distrib(3.7, mean = parameters7::simplex(3)), "'n_dim' must be")
})

test_that("the Dirichlet reads a missing row as missing and draws no empty row", {
  d <- dirichlet_distrib(3)
  th <- list(mean_alr1 = 0.3, mean_alr2 = -0.2, phi = 12)
  v <- distrib_pdf(d, rbind(c(NaN, NaN, NaN), c(NA, 0.5, 0.5), c(0.2, 0.3, 0.5)), th)
  expect_true(all(is.na(v[1:2])) && is.finite(v[3]))
  set.seed(1)
  r <- distrib_rng(d, 2000, list(mean_alr1 = 0.3, mean_alr2 = -0.2, phi = 0.001))
  expect_false(anyNA(r))
  expect_lt(max(abs(rowSums(r) - 1)), 1e-14)
})

test_that("the Dirichlet's expected higher orders are exact", {
  # references: E log y fixed at psi(alpha) - psi(phi), mpmath diff in theta
  d <- dirichlet_distrib(3)
  y <- rbind(c(0.2, 0.3, 0.5))
  e <- distrib_deriv3(d, y, list(mean_alr1 = 0.3, mean_alr2 = -0.2, phi = 0.02),
                      expected = TRUE)
  expect_rel(c(e$mean_alr1_mean_alr1_mean_alr1, e$mean_alr1_mean_alr2_phi, e$phi_phi_phi),
             c(-0.23846356819017141, 1.1992994370226738, 499998.01781172543),
             tolerance = 1e-10)
  e8 <- distrib_deriv3(d, y, list(mean_alr1 = 0.3, mean_alr2 = -0.2, phi = 8),
                       expected = TRUE)
  expect_rel(c(e8$mean_alr1_mean_alr1_mean_alr1, e8$mean_alr1_mean_alr2_phi, e8$phi_phi_phi),
             c(-0.72168791215910121, 0.10917417299768075, 0.0048747124378896325),
             tolerance = 1e-10)
})

test_that("the Dirichlet's derivatives at a tiny shape", {
  d <- dirichlet_distrib(3)
  th <- list(mean_alr1 = -300, mean_alr2 = -0.2, phi = 10)
  y <- rbind(c(0.2, 0.3, 0.5))
  d3 <- distrib_deriv3(d, y, th)
  d4 <- distrib_deriv4(d, y, th)
  expect_true(all(is.finite(unlist(d3))) && all(is.finite(unlist(d4))))
  # references at 400 digits: the components of the rare coordinate are of
  # the order of its shape, which the chain rule on log Gamma(alpha) lost
  expect_rel(c(d3$mean_alr1_mean_alr1_mean_alr1, d4$mean_alr1_mean_alr1_mean_alr2_phi),
             c(3.6853475222666919e-130, -2.9797559445734233e-131), tolerance = 1e-13)
  h <- distrib_hessian(d, y, list(mean_alr1 = -30, mean_alr2 = -0.2, phi = 10))
  expect_rel(h$mean_alr1_mean_alr1, 6.6986696578612316e-13, tolerance = 1e-13)
})

test_that("the multinomial at a probability that underflows", {
  d <- multinomial_distrib(3, size = 5)
  th <- list(probs_alr1 = -800, probs_alr2 = 0)
  y <- rbind(c(0, 2, 3))
  expect_rel(c(distrib_pdf(d, y, th, log = TRUE), distrib_gradient(d, y, th)$probs_alr2,
               distrib_hessian(d, y, th)$probs_alr2_probs_alr2),
             c(-1.1631508098056809, -0.5, -1.25), tolerance = 1e-14)
  d10 <- multinomial_distrib(3, size = 10)
  v <- unlist(c(distrib_deriv3(d10, rbind(c(0, 4, 6)), list(probs_alr1 = -400, probs_alr2 = 0)),
                distrib_deriv4(d10, rbind(c(0, 4, 6)), list(probs_alr1 = -200, probs_alr2 = 0))))
  expect_true(all(is.finite(v)))
  # the expected third order in closed form, -size times the cumulant
  e <- distrib_deriv3(d10, rbind(c(0, 4, 6)), list(probs_alr1 = 0.3, probs_alr2 = -0.2),
                      expected = TRUE)
  p <- exp(c(0.3, -0.2, 0)) / sum(exp(c(0.3, -0.2, 0)))
  expect_rel(e$probs_alr1_probs_alr1_probs_alr1, -10 * p[1] * (1 - p[1]) * (1 - 2 * p[1]),
             tolerance = 1e-13)
})

test_that("check_distrib() on a t without a covariance has no failed row", {
  d <- mvstudent_t1_distrib(2)
  th <- list(mu1 = 0.5, mu2 = -0.3, sigma_log_L1 = 0.1, sigma_log_L2 = -0.2,
             sigma_L2.1 = 0.4, nu = 1.5)
  set.seed(1)
  r <- check_distrib(d, theta = th, verbose = FALSE)
  expect_true(all(r$status == "OK"))
  expect_true("rng matches the first two moments" %in% attr(r, "skipped")$check)
})

test_that("mv_summary() reads the simplex of a Dirichlet and a multinomial", {
  set.seed(1)
  fd <- fit_distrib(dirichlet_distrib(3),
                    distrib_rng(dirichlet_distrib(3), 200,
                                list(mean_alr1 = 0.3, mean_alr2 = -0.2, phi = 8)))
  expect_true(all(c("p1", "p2", "p3") %in% rownames(mv_summary(fd))))
  dm <- multinomial_distrib(3, size = 10)
  fm <- fit_distrib(dm, distrib_rng(dm, 200, list(probs_alr1 = 0.3, probs_alr2 = -0.2)))
  expect_true(all(c("p1", "p2", "p3") %in% rownames(mv_summary(fm))))
})

test_that("the t's quantities of nu are exact at many degrees of freedom", {
  d3 <- mvstudent_t1_distrib(3)
  z3 <- list(mu1 = 0, mu2 = 0, mu3 = 0, sigma_log_L1 = 0, sigma_log_L2 = 0,
             sigma_log_L3 = 0, sigma_L2.1 = 0, sigma_L3.1 = 0, sigma_L3.2 = 0)
  expect_rel(distrib_pdf(d3, c(0, 0, 0), c(z3, list(nu = 1e12)), log = TRUE),
             -2.7568155996132682, tolerance = 1e-15)
  expect_rel(-distrib_expected_hessian(d3, matrix(0, 1, 3), c(z3, list(nu = 1e8)))$nu_nu,
             1.349999889000007e-31, tolerance = 1e-13)
  d2 <- mvstudent_t1_distrib(2)
  z2 <- list(mu1 = 0, mu2 = 0, sigma_log_L1 = 0, sigma_log_L2 = 0, sigma_L2.1 = 0)
  y <- rbind(c(1, -0.5))
  # d^k/dnu^k of the log-density at q = 1.25, p = 2 (vd_T4_4.py)
  expect_rel(distrib_hessian(d2, y, c(z2, list(nu = 1e100)))$nu_nu, -1.71875e-300,
             tolerance = 1e-13)
  expect_rel(distrib_deriv3(d2, y, c(z2, list(nu = 1e9)))$nu_nu_nu,
             5.156249989062500020751953e-36, tolerance = 1e-13)
  expect_rel(distrib_deriv4(d2, y, c(z2, list(nu = 1e9)))$nu_nu_nu_nu,
             -2.062499994531250012451172e-44, tolerance = 1e-13)
  # d/dnu of the response gradient, -y (q - p)/(nu + q)^2
  expect_rel(as.numeric(distrib_cross_y(d2, y, c(z2, list(nu = 1e12)))$nu),
             -c(1, -0.5) * (1.25 - 2) / (1e12 + 1.25)^2, tolerance = 1e-14)
  # the odd-dimension pair at the crossovers the old series used
  expect_rel(mvt_A(150, 3), -6.6078492311602689e-5, tolerance = 1e-14)
  expect_rel(2 * mvt_cp_d(100, 3, 2L), 2.9408900886328065e-6, tolerance = 1e-14)
})

test_that("the t's response derivatives at a distant observation", {
  d <- mvstudent_t1_distrib(1)
  th <- list(mu1 = 0, sigma_log_L1 = 0, nu = 5)
  # l_yy = 6 (y^2 - 5)/(5 + y^2)^2, l_yyy = 12 y (15 - y^2)/(5 + y^2)^3
  expect_rel(as.numeric(distrib_hess_y(d, matrix(1e100), th)), 6e-200, tolerance = 1e-13)
  expect_rel(as.numeric(distrib_deriv3_y(d, matrix(1e60), th)), -1.2e-179, tolerance = 1e-13)
  expect_rel(as.numeric(distrib_deriv4_y(d, matrix(1e60), th)), 3.6e-239, tolerance = 1e-13)
  expect_true(is.finite(as.numeric(distrib_hess_y(d, matrix(1e300), th))))
})

test_that("the response Hessian of an elliptical family is exactly symmetric", {
  # the inverse of the matrix is symmetrized once, in mv_matrix_pieces()
  for (d in list(mvgaussian1_distrib(3), mvstudent_t1_distrib(3))) {
    set.seed(4)
    th <- generate_random_theta(d)
    y <- distrib_rng(d, 6, th)
    h <- distrib_hess_y(d, y, th)
    expect_identical(h, if (length(dim(h)) == 2L) t(h) else aperm(h, c(2L, 1L, 3L)))
  }
})
