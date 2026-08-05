# gaussian1/2/3 and gamma1/2 are the same two laws in different coordinates.
# That is what makes them testable without choosing a tolerance: at
# corresponding parameters they must agree on the density, on the maximized
# log-likelihood and on AIC/BIC exactly, and their standard errors must
# correspond through the delta method.
#
# The comparison is only evidence because the implementations share nothing:
# each family has its own class, its own kernels and its own moments.

test_that("the three gaussians are one law", {
  set.seed(7)
  mu <- 1.3
  sg <- 2.1
  y <- rnorm(400, mu, sg)

  d1 <- gaussian1_distrib(); t1 <- list(mu = mu, sigma = sg)
  d2 <- gaussian2_distrib(); t2 <- list(mu = mu, sigma2 = sg^2)
  d3 <- gaussian3_distrib(); t3 <- list(mu = mu, tau = 1 / sg^2)

  expect_identical(distrib_pdf(d2, y, t2), distrib_pdf(d1, y, t1))
  expect_identical(distrib_pdf(d3, y, t3), distrib_pdf(d1, y, t1))
  expect_identical(distrib_cdf(d2, y, t2), distrib_cdf(d1, y, t1))
  expect_identical(distrib_cdf(d3, y, t3), distrib_cdf(d1, y, t1))

  p <- c(0.01, 0.25, 0.5, 0.75, 0.99)
  expect_identical(distrib_quantile(d2, p, t2), distrib_quantile(d1, p, t1))
  expect_identical(distrib_quantile(d3, p, t3), distrib_quantile(d1, p, t1))

  # the moments are the same numbers in each parametrization's coordinates
  expect_equal(variance(d2, t2), sg^2)
  expect_equal(variance(d3, t3), sg^2)
  expect_identical(skewness(d2, t2), 0)
  expect_identical(kurtosis(d3, t3), 0)
})


test_that("the three gaussians reach the same maximum", {
  set.seed(7)
  y <- rnorm(400, 1.3, 2.1)
  f1 <- fit_distrib(gaussian1_distrib(), y)
  f2 <- fit_distrib(gaussian2_distrib(), y)
  f3 <- fit_distrib(gaussian3_distrib(), y)

  for (f in list(f1, f2, f3)) expect_true(f@converged)

  expect_equal(as.numeric(logLik(f2)), as.numeric(logLik(f1)), tolerance = 1e-10)
  expect_equal(as.numeric(logLik(f3)), as.numeric(logLik(f1)), tolerance = 1e-10)
  expect_equal(AIC(f2), AIC(f1), tolerance = 1e-10)
  expect_equal(BIC(f3), BIC(f1), tolerance = 1e-10)

  # The estimates agree only as closely as the rule that stopped each run
  # promises, so the tolerance follows crit_grad()'s default rather than being
  # chosen: three optimizations of three different objectives land at three
  # nearby points, not at one.
  tol <- sqrt(eval(formals(optimizers7::crit_grad)$tol))
  s <- unname(coef(f1)[2])
  expect_equal(unname(coef(f2)[2]), s^2, tolerance = tol)
  expect_equal(unname(coef(f3)[2]), 1 / s^2, tolerance = tol)
  expect_equal(unname(coef(f2)[1]), unname(coef(f1)[1]), tolerance = tol)

  # and the standard errors correspond by the delta method, |dg/dsigma| se
  se <- unname(sqrt(diag(vcov(f1)))[2])
  expect_equal(unname(sqrt(diag(vcov(f2)))[2]), 2 * s * se, tolerance = tol)
  expect_equal(unname(sqrt(diag(vcov(f3)))[2]), 2 / s^3 * se, tolerance = tol)
})


test_that("gamma1 and gamma2 are one law", {
  set.seed(9)
  m <- 3
  ph <- 0.4
  y <- rgamma(500, shape = 1 / ph, rate = 1 / (ph * m))

  g1 <- gamma1_distrib(); t1 <- list(mu = m, phi = ph)
  g2 <- gamma2_distrib(); t2 <- list(mu = m, sigma2 = ph * m^2)

  expect_identical(distrib_pdf(g1, y, t1), distrib_pdf(g2, y, t2))
  expect_identical(distrib_cdf(g1, y, t1), distrib_cdf(g2, y, t2))

  expect_equal(variance(g1, t1), ph * m^2)
  expect_equal(skewness(g1, t1), skewness(g2, t2))
  expect_equal(kurtosis(g1, t1), kurtosis(g2, t2))

  f1 <- fit_distrib(g1, y)
  f2 <- fit_distrib(g2, y)
  expect_true(f1@converged)
  expect_equal(as.numeric(logLik(f1)), as.numeric(logLik(f2)), tolerance = 1e-8)
  expect_equal(unname(coef(f2)[2]),
               unname(coef(f1)[2] * coef(f1)[1]^2),
               tolerance = sqrt(eval(formals(optimizers7::crit_grad)$tol)))
})


test_that("each new family passes the whole validator", {
  cases <- list(
    list(gaussian2_distrib(), list(mu = 1, sigma2 = 4)),
    list(gaussian3_distrib(), list(mu = 1, tau = 0.25)),
    list(gamma1_distrib(), list(mu = 3, phi = 0.5))
  )
  for (cs in cases) {
    set.seed(1)
    res <- check_distrib(cs[[1]], cs[[2]], verbose = FALSE)
    expect_identical(nrow(res), 13L, label = cs[[1]]@distrib_name)
    expect_true(all(res$status == "OK"),
                info = paste(cs[[1]]@distrib_name,
                             paste(res$check[res$status != "OK"], collapse = ", ")))
  }
})


test_that("a broken gradient is still caught", {
  # The paired injection, so that passing the validator means something.
  good <- S7::method(distrib_gradient, Gaussian2Distrib)
  on.exit(S7::method(distrib_gradient, Gaussian2Distrib) <- good, add = TRUE)
  suppressMessages(
    S7::method(distrib_gradient, Gaussian2Distrib) <- function(distrib, y, theta,
                                                               scale = c("parameter", "link"), ...) {
      g <- good(distrib, y, theta, scale = scale, ...)
      g[[2]] <- g[[2]] * 1.05
      g
    }
  )
  set.seed(1)
  res <- check_distrib(gaussian2_distrib(), list(mu = 1, sigma2 = 4),
                       verbose = FALSE)
  expect_equal(res$status[res$check == "gradient vs finite differences"], "FAIL")
})


test_that("the third and fourth derivatives are the analytical ones", {
  # Against ONE Richardson differentiation of the analytical order below, which
  # is the rule the package follows everywhere. Note the component order: the
  # Hessian list is lexicographic while hess_names() puts the diagonal first,
  # so the components are looked up BY NAME and never by position.
  skip_if_not_installed("numDeriv")
  cases <- list(
    list(gaussian2_distrib(), list(mu = 1, sigma2 = 4)),
    list(gaussian3_distrib(), list(mu = 1, tau = 0.25)),
    list(gamma1_distrib(), list(mu = 3, phi = 0.5))
  )
  for (cs in cases) {
    d <- cs[[1]]
    th <- cs[[2]]
    set.seed(3)
    y <- distrib_rng(d, 30, th)
    p <- d@params
    d3 <- distrib_deriv3(d, y, th)
    d4 <- distrib_deriv4(d, y, th)

    for (hn in hess_names(p)) {
      for (k in seq_along(p)) {
        key <- paste(sort(c(hess_pairs(p)[[hn]], p[k])), collapse = "_")
        if (is.null(d3[[key]])) next
        f <- function(v) {
          q <- unlist(th)
          q[k] <- v
          sum(distrib_hessian(d, y, as.list(stats::setNames(q, p)))[[hn]])
        }
        expect_equal(sum(d3[[key]]), numDeriv::grad(f, th[[k]]),
                     tolerance = 1e-5,
                     label = paste(d@distrib_name, key))
      }
    }

    # order 4 against one differentiation of the analytical order 3
    for (n3 in names(d3)) {
      for (k in seq_along(p)) {
        key <- paste(sort(c(strsplit(n3, "(?<=\\D)_(?=\\D)", perl = TRUE)[[1]], p[k])),
                     collapse = "_")
        if (is.null(d4[[key]])) next
        f <- function(v) {
          q <- unlist(th)
          q[k] <- v
          sum(distrib_deriv3(d, y, as.list(stats::setNames(q, p)))[[n3]])
        }
        expect_equal(sum(d4[[key]]), numDeriv::grad(f, th[[k]]),
                     tolerance = 1e-4,
                     label = paste(d@distrib_name, key))
      }
    }
  }
})


test_that("the expected derivatives match a Monte Carlo of the observed ones", {
  cases <- list(
    list(gaussian2_distrib(), list(mu = 1, sigma2 = 4)),
    list(gaussian3_distrib(), list(mu = 1, tau = 0.25)),
    list(gamma1_distrib(), list(mu = 3, phi = 0.5))
  )
  for (cs in cases) {
    d <- cs[[1]]
    th <- cs[[2]]
    set.seed(11)
    y <- distrib_rng(d, 2e5, th)
    for (ord in 3:4) {
      a <- if (ord == 3) distrib_deriv3(d, y[1], th, expected = TRUE) else
        distrib_deriv4(d, y[1], th, expected = TRUE)
      o <- if (ord == 3) distrib_deriv3(d, y, th) else distrib_deriv4(d, y, th)
      scl <- max(abs(vapply(a, function(z) z[1], numeric(1))))
      for (nm in names(a)) {
        expect_lt(abs(a[[nm]][1] - mean(o[[nm]])), 0.02 * max(scl, 1e-8),
                  label = paste(d@distrib_name, ord, nm))
      }
    }
  }
})
