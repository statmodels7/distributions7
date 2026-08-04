# The two extreme-value families. Both are checked against routes the
# implementation does not take: the density against base R or against the
# formula written out, the derivatives against finite differences, the expected
# information against a large sample, and the moments against the generator.
# The two are also checked against each other, since exp(-Gumbel) is Weibull.

fd_grad <- function(f, x, h = 1e-6) {
  vapply(seq_along(x), function(k) {
    up <- dn <- x
    up[k] <- x[k] + h
    dn[k] <- x[k] - h
    (f(up) - f(dn)) / (2 * h)
  }, numeric(1))
}

# One stencil, never two nested first differences.
fd_hess <- function(f, x, k, l, h = 1e-4) {
  hk <- h * max(1, abs(x[k]))
  hl <- h * max(1, abs(x[l]))
  if (k == l) {
    up <- dn <- x
    up[k] <- x[k] + hk
    dn[k] <- x[k] - hk
    return((f(up) - 2 * f(x) + f(dn)) / hk^2)
  }
  pp <- pm <- mp <- mm <- x
  pp[k] <- pp[k] + hk; pp[l] <- pp[l] + hl
  pm[k] <- pm[k] + hk; pm[l] <- pm[l] - hl
  mp[k] <- mp[k] - hk; mp[l] <- mp[l] + hl
  mm[k] <- mm[k] - hk; mm[l] <- mm[l] - hl
  (f(pp) - f(pm) - f(mp) + f(mm)) / (4 * hk * hl)
}

expect_derivatives_ok <- function(d, th, y, tol_g = 1e-8, tol_h = 1e-5) {
  v0 <- unlist(th)
  ll <- function(v) {
    sum(distrib_pdf(d, y, as.list(stats::setNames(v, d@params)), log = TRUE))
  }
  # Richardson extrapolation, not the plain central difference of fd_grad():
  # a three-point stencil is good to about the two-thirds power of machine
  # precision, which is 1e-11 relative at best and worse than that once the
  # log-likelihood is summed over a sample. Asking a plain difference to agree
  # to 1e-8 is asking it for digits it does not have, and the same comparison
  # passed on Windows and failed on macOS.
  g <- vapply(distrib_gradient(d, y, th), sum, numeric(1))
  expect_equal(unname(g), numDeriv::grad(ll, v0), tolerance = tol_g)

  h <- distrib_hessian(d, y, th)
  expect_named(h, hess_names(d@params))
  pr <- hess_pairs(d@params)
  for (nm in hess_names(d@params)) {
    k <- match(pr[[nm]], d@params)
    expect_equal(sum(h[[nm]]), fd_hess(ll, v0, k[1], k[2]),
      tolerance = tol_h, label = nm
    )
  }
}


# --- Weibull ---------------------------------------------------------------

test_that("the Weibull density is base R's, on the same parametrisation", {
  d <- weibull_distrib()
  th <- list(mu = 2, sigma = 1.5)
  y <- c(0.4, 1, 2, 4)

  # mu is the SCALE, sigma the SHAPE; getting this backwards is the one
  # mistake this family invites.
  expect_equal(distrib_pdf(d, y, th), stats::dweibull(y, shape = 1.5, scale = 2))
  expect_equal(distrib_cdf(d, y, th), stats::pweibull(y, shape = 1.5, scale = 2))

  # and the log-density is the formula written out
  z <- y / 2
  want <- log(1.5) - log(2) + 0.5 * log(z) - z^1.5
  expect_equal(distrib_pdf(d, y, th, log = TRUE), want)

  # shape 1 is the exponential distribution with mean mu
  expect_equal(
    distrib_pdf(d, y, list(mu = 2, sigma = 1)),
    stats::dexp(y, rate = 0.5)
  )
})

test_that("the Weibull scale is not its mean", {
  d <- weibull_distrib()
  th <- list(mu = 2, sigma = 1.5)

  expect_equal(mean(d, th), 2 * gamma(1 + 1 / 1.5))
  expect_false(isTRUE(all.equal(mean(d, th), 2)))
  expect_equal(variance(d, th), 4 * (gamma(1 + 2 / 1.5) - gamma(1 + 1 / 1.5)^2))

  # shape 1: the exponential, whose mean and sd are both the scale
  expect_equal(mean(d, list(mu = 3, sigma = 1)), 3)
  expect_equal(variance(d, list(mu = 3, sigma = 1)), 9)
  expect_equal(skewness(d, list(mu = 3, sigma = 1)), 2)
  expect_equal(kurtosis(d, list(mu = 3, sigma = 1)), 6)

  # the shape alone fixes skewness and kurtosis
  expect_equal(
    skewness(d, list(mu = 1, sigma = 2)), skewness(d, list(mu = 100, sigma = 2))
  )
})

test_that("the Weibull derivatives agree with finite differences", {
  d <- weibull_distrib()
  set.seed(41)
  for (th in list(list(mu = 2, sigma = 1.5), list(mu = 0.3, sigma = 4),
                  list(mu = 10, sigma = 0.7))) {
    y <- distrib_rng(d, 25, th)
    expect_derivatives_ok(d, th, y)
  }
})

test_that("the Weibull expected information is the closed form it claims", {
  d <- weibull_distrib()
  th <- list(mu = 2, sigma = 1.5)
  eg <- -digamma(1)

  eh <- distrib_expected_hessian(d, 1, th)
  expect_equal(eh$mu_mu[1], -1.5^2 / 4)
  expect_equal(eh$sigma_sigma[1], -((1 - eg)^2 + pi^2 / 6) / 1.5^2)
  expect_equal(eh$mu_sigma[1], (1 - eg) / 2)

  # against the sample average of the observed Hessian, which is a different
  # quantity that happens to coincide because the family is regular
  set.seed(42)
  big <- distrib_rng(d, 3e5, th)
  hb <- distrib_hessian(d, big, th)
  for (nm in names(eh)) {
    se <- stats::sd(hb[[nm]]) / sqrt(length(hb[[nm]]))
    expect_lt(abs(eh[[nm]][1] - mean(hb[[nm]])), 4 * se)
  }

  # the location and the scale are NOT orthogonal here: the mixed block is a
  # fixed non-zero number, unlike in a symmetric location-scale family
  expect_gt(abs(eh$mu_sigma[1]), 0.1)
})

test_that("the Weibull response derivatives are closed form", {
  d <- weibull_distrib()
  th <- list(mu = 2, sigma = 1.5)
  y <- c(0.4, 1, 2, 4)
  u <- (y / 2)^1.5

  expect_equal(distrib_grad_y(d, y, th), (1.5 - 1 - 1.5 * u) / y)
  expect_equal(distrib_hess_y(d, y, th), -(1.5 - 1) * (1 + 1.5 * u) / y^2)

  # shape 1: the exponential has a linear log-density, so the second
  # derivative in y vanishes identically
  expect_equal(distrib_hess_y(d, y, list(mu = 2, sigma = 1)), rep(0, length(y)))
})


# --- Gumbel ----------------------------------------------------------------

test_that("the Gumbel density and distribution function are the formulas", {
  d <- gumbel_distrib()
  th <- list(mu = 1, sigma = 2)
  y <- c(-2, 0, 1, 4)
  z <- (y - 1) / 2

  expect_equal(distrib_pdf(d, y, th, log = TRUE), -log(2) - z - exp(-z))
  expect_equal(distrib_cdf(d, y, th), exp(-exp(-z)))
  expect_equal(distrib_quantile(d, c(0.1, 0.5, 0.9), th),
               1 - 2 * log(-log(c(0.1, 0.5, 0.9))))

  # the lower tail is exact on the log scale, where exp(-exp(-z)) would
  # underflow to zero
  expect_equal(distrib_cdf(d, -400, th, log.p = TRUE), -exp(200.5))
  expect_equal(distrib_cdf(d, -400, th), 0)
})

test_that("the Gumbel has a fixed shape", {
  d <- gumbel_distrib()
  eg <- -digamma(1)

  expect_equal(mean(d, list(mu = 1, sigma = 2)), 1 + 2 * eg)
  expect_equal(variance(d, list(mu = 1, sigma = 2)), pi^2 * 4 / 6)

  # the third and fourth standardised moments do not depend on either
  # parameter, which is the substantive statement about this family
  s1 <- skewness(d, list(mu = 0, sigma = 1))
  s2 <- skewness(d, list(mu = -50, sigma = 7))
  expect_equal(s1, s2)
  expect_equal(s1, 1.1395, tolerance = 1e-4)
  expect_equal(kurtosis(d, list(mu = 0, sigma = 1)), 12 / 5)
  expect_equal(kurtosis(d, list(mu = -50, sigma = 7)), 12 / 5)
})

test_that("the Gumbel derivatives agree with finite differences", {
  d <- gumbel_distrib()
  set.seed(43)
  for (th in list(list(mu = 1, sigma = 2), list(mu = -3, sigma = 0.5),
                  list(mu = 0, sigma = 5))) {
    y <- distrib_rng(d, 25, th)
    expect_derivatives_ok(d, th, y)
  }
})

test_that("the Gumbel expected information is the closed form it claims", {
  d <- gumbel_distrib()
  th <- list(mu = 1, sigma = 2)
  eg <- -digamma(1)

  eh <- distrib_expected_hessian(d, 0, th)
  expect_equal(eh$mu_mu[1], -1 / 4)
  expect_equal(eh$sigma_sigma[1], -((1 - eg)^2 + pi^2 / 6) / 4)
  expect_equal(eh$mu_sigma[1], (1 - eg) / 4)

  set.seed(44)
  big <- distrib_rng(d, 3e5, th)
  hb <- distrib_hessian(d, big, th)
  for (nm in names(eh)) {
    se <- stats::sd(hb[[nm]]) / sqrt(length(hb[[nm]]))
    expect_lt(abs(eh[[nm]][1] - mean(hb[[nm]])), 4 * se)
  }

  # the mixed block does not vanish: the density is skewed, so the location and
  # the scale are not orthogonal
  expect_gt(abs(eh$mu_sigma[1]), 0.05)
})

test_that("the Gumbel response derivatives are minus the location ones", {
  d <- gumbel_distrib()
  th <- list(mu = 1, sigma = 2)
  y <- c(-2, 0, 1, 4)

  # a location family: d/dy = -d/dmu, exactly
  expect_equal(distrib_grad_y(d, y, th), -distrib_gradient(d, y, th)$mu)
  expect_equal(distrib_hess_y(d, y, th), distrib_hessian(d, y, th)$mu_mu)
})


# --- the two together ------------------------------------------------------

test_that("exp(-Gumbel) is Weibull", {
  # The families are one another on the log scale, so the transformation
  # wrapper must reproduce the Weibull exactly rather than approximately.
  dg <- gumbel_distrib()
  dw <- weibull_distrib()
  thg <- list(mu = 1, sigma = 2)
  thw <- list(mu = exp(-1), sigma = 1 / 2)

  x <- c(0.05, 0.2, 0.6, 1.5)
  # P(exp(-Y) <= x) = P(Y >= -log x) = 1 - F_G(-log x)
  expect_equal(
    distrib_cdf(dw, x, thw),
    distrib_cdf(dg, -log(x), thg, lower.tail = FALSE)
  )

  set.seed(45)
  draws <- exp(-distrib_rng(dg, 20000, thg))
  expect_gt(
    suppressWarnings(stats::ks.test(
      draws, function(q) distrib_cdf(dw, q, thw)
    )$p.value),
    0.01
  )
})

test_that("check_distrib passes on both, and catches a wrong gradient", {
  for (cfg in list(list(d = weibull_distrib(), th = list(mu = 2, sigma = 1.5)),
                   list(d = gumbel_distrib(), th = list(mu = 1, sigma = 2)))) {
    set.seed(46)
    res <- check_distrib(cfg$d, theta = cfg$th, nsim = 3e4, verbose = FALSE)
    expect_true(all(res$status == "OK"), label = cfg$d@distrib_name)
  }

  # a gradient 5% wrong must still be caught, or the agreements above prove
  # nothing
  Wrong <- S7::new_class("WrongWeibull", parent = WeibullDistrib, package = NULL)
  gen <- distrib_gradient
  S7::method(gen, Wrong) <- function(distrib, y, theta,
                                     scale = c("parameter", "link"), ...) {
    g <- S7::method(distrib_gradient, WeibullDistrib)(distrib, y, theta)
    g[["mu"]] <- 1.05 * g[["mu"]]
    g
  }
  good <- weibull_distrib()
  bad <- Wrong(
    distrib_name = "wrong", dimension = "univariate", bounds = good@bounds,
    params = good@params, params_interpretation = good@params_interpretation,
    n_params = good@n_params, params_bounds = good@params_bounds,
    link_params = good@link_params
  )
  set.seed(47)
  res <- check_distrib(bad, theta = list(mu = 2, sigma = 1.5), nsim = 1e4,
                       verbose = FALSE)
  expect_identical(
    res$status[res$check == "gradient vs finite differences"], "FAIL"
  )
})

test_that("both families are fitted by maximum likelihood", {
  for (cfg in list(list(d = weibull_distrib(), th = list(mu = 2, sigma = 1.5)),
                   list(d = gumbel_distrib(), th = list(mu = 1, sigma = 2)))) {
    set.seed(48)
    y <- distrib_rng(cfg$d, 3000, cfg$th)
    fit <- fit_distrib(cfg$d, y)

    expect_true(fit@converged)
    est <- coef(fit)
    expect_equal(unname(est), unname(unlist(cfg$th)), tolerance = 0.1)

    # the optimum is a stationary point on the link scale
    sc <- vapply(
      distrib_gradient(cfg$d, y, as.list(est), scale = "link"), sum, numeric(1)
    )
    expect_lt(max(abs(sc)) / length(y), 1e-6)

    # and the three methods agree on the maximised log-likelihood
    lls <- vapply(c("fisher", "newton", "bfgs"), function(m) {
      as.numeric(logLik(fit_distrib(cfg$d, y, method = m)))
    }, numeric(1))
    expect_equal(diff(range(lls)), 0, tolerance = 1e-6)
  }
})
