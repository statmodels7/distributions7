# The two skew families. The skew normal is closed form throughout; the skew t
# is closed form except where the density's tilting factor is a Student t
# distribution function whose derivative in its degrees of freedom does not
# exist in elementary terms. Both are checked against routes the implementation
# does not take, and the special cases they contain are checked against the
# families already in the package.

# numDeriv's Richardson extrapolation is used as the reference for the skew t,
# rather than a plain difference: the components of its score that are
# themselves finite differences must be compared against something better than
# the same arithmetic.

test_that("the skew normal contains the gaussian", {
  d <- skewnormal1_distrib()
  y <- c(-2, -0.5, 0, 1, 3)

  expect_equal(
    distrib_pdf(d, y, list(mu = 0, sigma = 1, alpha = 0)), stats::dnorm(y)
  )
  expect_equal(
    distrib_cdf(d, y, list(mu = 2, sigma = 3, alpha = 0)),
    stats::pnorm(y, mean = 2, sd = 3)
  )
  # and the derivatives reduce too
  g <- distrib_gradient(d, y, list(mu = 0, sigma = 1, alpha = 0))
  gg <- distrib_gradient(gaussian1_distrib(), y, list(mu = 0, sigma = 1))
  expect_equal(g$mu, gg$mu)
  expect_equal(g$sigma, gg$sigma)
})

test_that("the skew normal density is the formula and integrates to one", {
  d <- skewnormal1_distrib()
  th <- list(mu = 0.5, sigma = 2, alpha = -2.5)
  y <- c(-2, -0.5, 0, 1, 3)
  z <- (y - 0.5) / 2

  expect_equal(
    distrib_pdf(d, y, th, log = TRUE),
    log(2) - log(2) + stats::dnorm(z, log = TRUE) + stats::pnorm(-2.5 * z, log.p = TRUE)
  )
  expect_equal(
    stats::integrate(function(t) distrib_pdf(d, t, th), -Inf, Inf)$value, 1,
    tolerance = 1e-8
  )
})

test_that("Owen's T gives the same distribution function as quadrature", {
  d <- skewnormal1_distrib()
  th <- list(mu = 0.5, sigma = 2, alpha = -2.5)
  y <- c(-2, -0.5, 0, 1, 3)

  # the implementation uses Owen's T; the reference integrates the density,
  # which is a different computation
  num <- vapply(y, function(q) {
    stats::integrate(function(t) distrib_pdf(d, t, th), -Inf, q)$value
  }, numeric(1))
  expect_equal(distrib_cdf(d, y, th), num, tolerance = 1e-8)

  # the two identities the implementation relies on
  expect_equal(owen_t(c(0.5, -1.2), 0), c(0, 0))
  expect_equal(owen_t(1.3, Inf), 0.5 * stats::pnorm(-1.3))
  expect_equal(owen_t(1.3, -2), -owen_t(1.3, 2))
})

test_that("the inverse Mills ratio survives the tail both densities underflow in", {
  # phi(t) and Phi(t) both underflow below about t = -38, while their ratio is
  # finite there and close to -t. Formed directly it would be 0/0.
  m <- mills_ratio(c(-60, -40, -10, 0, 10))
  expect_true(all(is.finite(m$r)))
  expect_true(all(is.finite(m$dr)))
  expect_equal(m$r[1], 60, tolerance = 1e-3)
  expect_equal(m$r[4], sqrt(2 / pi))

  # R' = -R(t + R), checked against a difference of R itself
  h <- 1e-5
  num <- (mills_ratio(-2 + h)$r - mills_ratio(-2 - h)$r) / (2 * h)
  expect_equal(mills_ratio(-2)$dr, num, tolerance = 1e-7)
})

test_that("the skew normal derivatives agree with an independent reference", {
  d <- skewnormal1_distrib()
  set.seed(51)
  for (th in list(list(mu = 0, sigma = 1, alpha = 3),
                  list(mu = -1, sigma = 2.5, alpha = -4),
                  list(mu = 2, sigma = 0.4, alpha = 0.001))) {
    y <- distrib_rng(d, 25, th)
    v0 <- unlist(th)
    ll <- function(v) {
      sum(distrib_pdf(d, y, as.list(stats::setNames(v, d@params)), log = TRUE))
    }
    g <- vapply(distrib_gradient(d, y, th), sum, numeric(1))
    expect_equal(unname(g), numDeriv::grad(ll, v0), tolerance = 1e-7)

    # the Hessian against ONE differentiation of the analytic gradient
    gfun <- function(v) {
      vapply(distrib_gradient(d, y, as.list(stats::setNames(v, d@params))),
             sum, numeric(1))
    }
    j <- numDeriv::jacobian(gfun, v0)
    h <- distrib_hessian(d, y, th)
    pr <- hess_pairs(d@params)
    for (nm in hess_names(d@params)) {
      k <- match(pr[[nm]], d@params)
      expect_equal(sum(h[[nm]]), (j[k[1], k[2]] + j[k[2], k[1]]) / 2,
        tolerance = 1e-6, label = paste(nm, th$alpha)
      )
    }
  }
})

test_that("the skew normal's moments are bounded, as the family requires", {
  d <- skewnormal1_distrib()
  th <- list(mu = 0, sigma = 1, alpha = 3)

  set.seed(52)
  big <- distrib_rng(d, 3e5, th)
  expect_equal(mean(d, th), mean(big), tolerance = 0.01)
  expect_equal(variance(d, th), stats::var(big), tolerance = 0.01)
  expect_equal(skewness(d, th),
    mean((big - mean(big))^3) / stats::sd(big)^3, tolerance = 0.05)

  # the bounds are the reason the skew t exists
  expect_lt(skewness(d, list(mu = 0, sigma = 1, alpha = 1e6)), 0.9953)
  expect_lt(kurtosis(d, list(mu = 0, sigma = 1, alpha = 1e6)), 0.8692)
  expect_equal(skewness(d, list(mu = 0, sigma = 1, alpha = 0)), 0)
})

test_that("the skew normal response derivatives are the location ones", {
  d <- skewnormal1_distrib()
  th <- list(mu = 0, sigma = 1, alpha = 3)
  y <- c(-2, 0, 1)
  expect_equal(distrib_grad_y(d, y, th), -distrib_gradient(d, y, th)$mu)
  expect_equal(distrib_hess_y(d, y, th), distrib_hessian(d, y, th)$mu_mu)
})


# --- the skew t ------------------------------------------------------------

test_that("the skew t contains the Student t and the skew normal", {
  d <- skewt_distrib()
  y <- c(-3, -0.5, 0, 1.2, 4)

  expect_equal(
    distrib_pdf(d, y, list(mu = 0, sigma = 1, alpha = 0, nu = 5)),
    stats::dt(y, df = 5)
  )
  expect_equal(
    distrib_pdf(d, y, list(mu = 0, sigma = 1, alpha = 3, nu = 1e7)),
    distrib_pdf(skewnormal1_distrib(), y, list(mu = 0, sigma = 1, alpha = 3)),
    tolerance = 1e-6
  )
  expect_equal(
    stats::integrate(function(t) {
      distrib_pdf(d, t, list(mu = 0, sigma = 1, alpha = 3, nu = 5))
    }, -Inf, Inf)$value,
    1, tolerance = 1e-9
  )
})

test_that("the skew t score agrees with an independent reference in every component", {
  d <- skewt_distrib()
  th <- list(mu = 0.5, sigma = 2, alpha = -3, nu = 6)
  y <- c(-3, -0.5, 0, 1.2, 4)
  v0 <- unlist(th)
  ll <- function(v) {
    sum(distrib_pdf(d, y, as.list(stats::setNames(v, d@params)), log = TRUE))
  }
  g <- vapply(distrib_gradient(d, y, th), sum, numeric(1))

  # Richardson extrapolation, not a plain difference: the nu component of the
  # score IS a difference, so a plain one would be the same arithmetic twice.
  expect_equal(unname(g), numDeriv::grad(ll, v0), tolerance = 1e-7)

  # the three closed-form components are exact, and are worth separating from
  # the fourth
  expect_equal(g[["mu"]], numDeriv::grad(ll, v0)[1], tolerance = 1e-9)
  expect_equal(g[["alpha"]], numDeriv::grad(ll, v0)[3], tolerance = 1e-9)
})

test_that("the skew t Hessian agrees with one differentiation of its score", {
  d <- skewt_distrib()
  th <- list(mu = 0.5, sigma = 2, alpha = -3, nu = 6)
  y <- c(-3, -0.5, 0, 1.2, 4)
  v0 <- unlist(th)
  gfun <- function(v) {
    vapply(distrib_gradient(d, y, as.list(stats::setNames(v, d@params))),
           sum, numeric(1))
  }
  j <- numDeriv::jacobian(gfun, v0)
  h <- distrib_hessian(d, y, th)
  pr <- hess_pairs(d@params)

  expect_named(h, hess_names(d@params))
  for (nm in hess_names(d@params)) {
    k <- match(pr[[nm]], d@params)
    # the components involving nu carry the stencil's error, which is what the
    # constructor's table states; the rest are exact
    tol <- if (any(pr[[nm]] == "nu")) 1e-5 else 1e-8
    expect_equal(sum(h[[nm]]), (j[k[1], k[2]] + j[k[2], k[1]]) / 2,
      tolerance = tol, label = nm
    )
  }
})

test_that("the skew t moments exist only up to nu", {
  d <- skewt_distrib()
  th <- list(mu = 0, sigma = 1, alpha = 3, nu = 8)

  set.seed(53)
  big <- distrib_rng(d, 4e5, th)
  expect_equal(mean(d, th), mean(big), tolerance = 0.02)
  expect_equal(variance(d, th), stats::var(big), tolerance = 0.02)
  expect_equal(skewness(d, th),
    mean((big - mean(big))^3) / stats::sd(big)^3, tolerance = 0.1)

  # each moment is NaN below its own threshold, while the density is finite
  expect_true(is.finite(distrib_pdf(d, 0, list(mu = 0, sigma = 1, alpha = 3, nu = 0.5))))
  expect_true(is.nan(mean(d, list(mu = 0, sigma = 1, alpha = 3, nu = 0.5))))
  expect_true(is.finite(mean(d, list(mu = 0, sigma = 1, alpha = 3, nu = 1.5))))
  expect_true(is.nan(variance(d, list(mu = 0, sigma = 1, alpha = 3, nu = 1.5))))
  expect_true(is.nan(skewness(d, list(mu = 0, sigma = 1, alpha = 3, nu = 2.5))))
  expect_true(is.nan(kurtosis(d, list(mu = 0, sigma = 1, alpha = 3, nu = 3.5))))

  # and it reaches skewness the skew normal cannot
  expect_gt(skewness(d, list(mu = 0, sigma = 1, alpha = 8, nu = 5)), 0.9953)
})

test_that("the skew t generator matches its own density", {
  d <- skewt_distrib()
  th <- list(mu = 0, sigma = 1, alpha = 3, nu = 5)
  set.seed(54)
  r <- distrib_rng(d, 5000, th)
  expect_gt(
    suppressWarnings(stats::ks.test(r, function(q) distrib_cdf(d, q, th))$p.value),
    0.01
  )
})

test_that("check_distrib passes on both skew families", {
  set.seed(55)
  res <- check_distrib(skewnormal1_distrib(),
    theta = list(mu = 0, sigma = 1, alpha = 3), nsim = 2e4, verbose = FALSE)
  expect_true(all(res$status == "OK"))

  set.seed(56)
  res <- check_distrib(skewt_distrib(),
    theta = list(mu = 0, sigma = 1, alpha = 3, nu = 6), nsim = 1e4,
    verbose = FALSE)
  expect_true(all(res$status == "OK"))
})

test_that("both skew families are fitted by maximum likelihood", {
  d <- skewnormal1_distrib()
  set.seed(57)
  y <- distrib_rng(d, 3000, list(mu = 1, sigma = 2, alpha = 4))
  fit <- fit_distrib(d, y)
  expect_true(fit@converged)
  expect_equal(unname(coef(fit)), c(1, 2, 4), tolerance = 0.3)
  sc <- vapply(
    distrib_gradient(d, y, as.list(coef(fit)), scale = "link"), sum, numeric(1)
  )
  # The bound follows the stopping rule's own default, which is a tolerance on
  # the score per observation; naming the constant here instead would go stale
  # the moment the default moved.
  expect_lt(max(abs(sc)) / length(y), eval(formals(optimizers7::crit_grad)$tol))

  # The skew t's score in nu carries a finite difference, whose error on the
  # score per observation is around 1e-11, so the default tolerance of 1e-6 is
  # comfortably above it and the run stops where it should.
  dt4 <- skewt_distrib()
  set.seed(58)
  y4 <- distrib_rng(dt4, 2000, list(mu = 1, sigma = 2, alpha = 4, nu = 6))
  fit4 <- fit_distrib(dt4, y4)
  expect_true(fit4@converged)
  expect_equal(unname(coef(fit4))[1:3], c(1, 2, 4), tolerance = 0.4)

  # and at the default tolerance the fit is still RETURNED rather than raising:
  # a run that reached a point is not a run that failed
  fit_default <- fit_distrib(dt4, y4, method = fisher_scoring(maxit = 40))
  expect_s3_class(coef(fit_default), NA)
  expect_equal(as.numeric(logLik(fit_default)), as.numeric(logLik(fit4)),
    tolerance = 1e-4
  )
})
