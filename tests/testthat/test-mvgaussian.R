# The multivariate gaussian, and with it the contract multivariate
# distributions follow. Every closed form is checked against a route the
# implementation does not take: the density against the formula written out by
# hand, the derivatives against finite differences, the expected information
# against Monte Carlo, and the fit against the maximum likelihood estimator in
# closed form.

fd_grad <- function(f, x, h = 1e-5) {
  vapply(seq_along(x), function(k) {
    up <- dn <- x
    up[k] <- x[k] + h
    dn[k] <- x[k] - h
    (f(up) - f(dn)) / (2 * h)
  }, numeric(1))
}

# A second derivative from ONE stencil. Composing two first differences is the
# obvious route and it is the one this toolkit forbids everywhere else: the
# reference becomes the error of an error, and whether it lands inside a
# tolerance is then a fact about the platform's arithmetic. Written out here
# because a test must not reach into the package for the thing it is checking.
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

test_that("the constructor validates its arguments", {
  expect_error(mvgaussian1_distrib(0), "positive integer")
  expect_error(mvgaussian1_distrib(2.5), "positive integer")
  # The side is the family's, so neither constructor takes the other's
  # argument: what used to be an over-determined call is now unwritable.
  expect_error(
    mvgaussian1_distrib(2, omega = parameters7::log_cholesky(2)),
    "unused argument"
  )
  expect_error(
    mvgaussian2_distrib(2, sigma = parameters7::log_cholesky(2)),
    "unused argument"
  )
  expect_error(mvgaussian1_distrib(2, diag(2)), "parameter")
  expect_error(
    mvgaussian1_distrib(3, parameters7::log_cholesky(2)),
    "dimension 2 but the distribution has dimension 3"
  )

  # A rank-deficient structure is a penalty, not a density, and the refusal is
  # the same whichever side it parametrizes.
  pen <- parameters7::scaled_matrix(crossprod(diff(diag(5), differences = 2)))
  expect_error(mvgaussian1_distrib(5, pen), "rank deficient")
  expect_error(mvgaussian2_distrib(5, pen), "rank deficient")
})

test_that("the parameters are the mean components and the structure's free values", {
  d <- mvgaussian1_distrib(2)

  expect_identical(d@params, c("mu1", "mu2", "sigma_log_L1", "sigma_log_L2", "sigma_L2.1"))
  expect_identical(d@n_params, 5L)
  expect_identical(d@n_dim, 2L)
  expect_identical(d@dimension, "multivariate")
  expect_true(S7::S7_inherits(d, multivariate_distrib))
  # Not a continuous_distrib: the one-dimensional defaults registered there --
  # a cdf by quadrature, a quantile by root finding -- have no counterpart.
  expect_false(S7::S7_inherits(d, continuous_distrib))

  # Every parameter is unconstrained, so every link is the identity: the
  # constraint lives inside the structure.
  expect_true(all(vapply(
    d@link_params, function(l) l@link_name == "identity", logical(1)
  )))
  expect_true(all(vapply(d@params_bounds, function(b) all(!is.finite(b)), logical(1))))

  # a diagonal structure contributes fewer of them
  dd <- mvgaussian1_distrib(3, parameters7::diagonal_matrix(3))
  expect_identical(dd@params,
    c("mu1", "mu2", "mu3", "sigma_log_d1", "sigma_log_d2", "sigma_log_d3"))
})

test_that("the density is the formula, written out by hand", {
  d <- mvgaussian1_distrib(2)
  th <- list(mu1 = 0.5, mu2 = -0.3, sigma_log_L1 = 0.1, sigma_log_L2 = -0.2, sigma_L2.1 = 0.4)
  y <- rbind(c(0, 0), c(1, -1), c(-0.5, 0.8))

  mu <- mv_location(d, th)
  s <- mv_sigma(d, th)
  want <- apply(y, 1, function(row) {
    r <- row - mu
    -0.5 * (2 * log(2 * pi) + log(det(s)) + drop(t(r) %*% solve(s) %*% r))
  })

  expect_equal(distrib_pdf(d, y, th, log = TRUE), unname(want))
  expect_equal(distrib_pdf(d, y, th), unname(exp(want)))

  # a vector of length p is one observation
  expect_equal(distrib_pdf(d, c(0, 0), th, log = TRUE), unname(want[1]))
  # and one of the wrong length is refused rather than recycled
  expect_error(distrib_pdf(d, c(0, 0, 0), th), "read as one observation")
  expect_error(distrib_pdf(d, cbind(1, 2, 3), th), "3 column")

  expect_equal(distrib_pdf(d, matrix(numeric(0), 0, 2), th), numeric(0))
})

test_that("the mean and the covariance come back in the shapes they belong to", {
  d <- mvgaussian1_distrib(2)
  th <- list(mu1 = 1, mu2 = -1, sigma_log_L1 = 0, sigma_log_L2 = 0, sigma_L2.1 = 0.5)

  expect_equal(unname(mv_location(d, th)), c(1, -1))
  s <- mv_sigma(d, th)
  expect_identical(dim(s), c(2L, 2L))
  # Sigma = L L' with the factor the structure assembles: L is LOWER
  # triangular, so the free value sigma_L2.1 sits below the diagonal and the variance
  # it inflates is the second one.
  l <- matrix(c(1, 0, 0.5, 1), 2, 2, byrow = TRUE)
  expect_equal(unname(s), tcrossprod(l))
  expect_equal(unname(s), matrix(c(1, 0.5, 0.5, 1.25), 2, 2))

  expect_equal(mean(d, th), mv_location(d, th))
  expect_equal(variance(d, th), s)
})

test_that("a parameter that varies by observation is refused", {
  d <- mvgaussian1_distrib(2)
  th <- list(mu1 = c(0, 1), mu2 = 0, sigma_log_L1 = 0, sigma_log_L2 = 0, sigma_L2.1 = 0)
  expect_error(distrib_pdf(d, rbind(c(0, 0), c(1, 1)), th), "belong to a model")
})

test_that("the score and the Hessian agree with finite differences", {
  d <- mvgaussian1_distrib(2)
  th <- list(mu1 = 0.5, mu2 = -0.3, sigma_log_L1 = 0.1, sigma_log_L2 = -0.2, sigma_L2.1 = 0.4)
  set.seed(21)
  y <- distrib_rng(d, 40, th)
  v0 <- unlist(th)
  ll <- function(v) {
    sum(distrib_pdf(d, y, as.list(stats::setNames(v, d@params)), log = TRUE))
  }

  g <- vapply(distrib_gradient(d, y, th), sum, numeric(1))
  expect_equal(unname(g), fd_grad(ll, v0), tolerance = 1e-6)
  expect_named(distrib_gradient(d, y, th), d@params)

  h <- distrib_hessian(d, y, th)
  expect_named(h, hess_names(d@params))
  pr <- hess_pairs(d@params)
  for (nm in hess_names(d@params)) {
    k <- match(pr[[nm]], d@params)
    fd <- fd_hess(ll, v0, k[1], k[2])
    expect_equal(sum(h[[nm]]), fd, tolerance = 1e-4, label = nm)
  }
})

test_that("the expected information is exact and block diagonal", {
  d <- mvgaussian1_distrib(2)
  th <- list(mu1 = 0.5, mu2 = -0.3, sigma_log_L1 = 0.1, sigma_log_L2 = -0.2, sigma_L2.1 = 0.4)
  set.seed(22)
  big <- distrib_rng(d, 2e5, th)

  eh <- distrib_expected_hessian(d, big[1, , drop = FALSE], th)
  hb <- distrib_hessian(d, big, th)
  for (nm in hess_names(d@params)) {
    expect_equal(eh[[nm]][1], mean(hb[[nm]]),
      tolerance = 0.05, label = nm
    )
  }

  # The mean and the matrix parameters are orthogonal: E[w] = 0 kills the mixed
  # block exactly, not approximately.
  expect_equal(eh[["mu1_sigma_log_L1"]][1], 0)
  expect_equal(eh[["mu2_sigma_L2.1"]][1], 0)

  # and the mean block is -Sigma^{-1}
  si <- solve(mv_sigma(d, th))
  expect_equal(eh[["mu1_mu1"]][1], -si[1, 1])
  expect_equal(eh[["mu1_mu2"]][1], -si[1, 2])
})

test_that("the generator matches the first two moments", {
  d <- mvgaussian1_distrib(2)
  th <- list(mu1 = 1.5, mu2 = -0.5, sigma_log_L1 = 0.2, sigma_log_L2 = -0.1, sigma_L2.1 = 0.6)
  set.seed(23)
  r <- distrib_rng(d, 2e5, th)

  expect_identical(dim(r), c(200000L, 2L))
  expect_equal(unname(colMeans(r)), unname(mv_location(d, th)), tolerance = 0.02)
  expect_equal(unname(stats::cov(r)), unname(mv_sigma(d, th)), tolerance = 0.02)
})

test_that("the response derivatives are closed form", {
  d <- mvgaussian1_distrib(2)
  th <- list(mu1 = 0.5, mu2 = -0.3, sigma_log_L1 = 0.1, sigma_log_L2 = -0.2, sigma_L2.1 = 0.4)
  y <- rbind(c(0, 0), c(1, -1))
  si <- solve(mv_sigma(d, th))

  gy <- distrib_grad_y(d, y, th)
  want <- -sweep(y, 2L, mv_location(d, th)) %*% si
  expect_equal(unname(gy), unname(want))
  expect_equal(unname(distrib_hess_y(d, y, th)), unname(-si))
})

test_that("the link scale is the parameter scale", {
  d <- mvgaussian1_distrib(2)
  th <- list(mu1 = 0.5, mu2 = -0.3, sigma_log_L1 = 0.1, sigma_log_L2 = -0.2, sigma_L2.1 = 0.4)
  set.seed(24)
  y <- distrib_rng(d, 20, th)

  expect_equal(
    distrib_gradient(d, y, th, scale = "link"),
    distrib_gradient(d, y, th)
  )
  expect_equal(
    distrib_hessian(d, y, th, scale = "link"),
    distrib_hessian(d, y, th)
  )
})

test_that("the one-dimensional quantities are refused rather than approximated", {
  d <- mvgaussian1_distrib(2)
  th <- list(mu1 = 0, mu2 = 0, sigma_log_L1 = 0, sigma_log_L2 = 0, sigma_L2.1 = 0)

  expect_error(distrib_cdf(d, rbind(c(0, 0)), th), "orthant")
  expect_error(distrib_quantile(d, 0.5, th), "ordering")
})

test_that("the precision form describes the same law as the covariance form", {
  # Omega = Sigma^{-1} through a structure on the other side: the density, the
  # score and the Hessian must agree once the two are matched.
  ds <- mvgaussian1_distrib(2)
  do <- mvgaussian2_distrib(2, parameters7::log_cholesky(2))

  th_s <- list(mu1 = 0.3, mu2 = -0.4, sigma_log_L1 = 0.1, sigma_log_L2 = -0.2, sigma_L2.1 = 0.5)
  sigma <- mv_sigma(ds, th_s)
  # the free values of the precision that give this covariance
  eta_o <- parameters7::param_free(do@param, solve(sigma))
  th_o <- as.list(stats::setNames(c(0.3, -0.4, unname(eta_o)), do@params))

  expect_equal(mv_sigma(do, th_o), sigma, tolerance = 1e-10)

  set.seed(25)
  y <- distrib_rng(ds, 30, th_s)
  expect_equal(
    distrib_pdf(do, y, th_o, log = TRUE),
    distrib_pdf(ds, y, th_s, log = TRUE)
  )
  # the score in the mean is the same quantity in both parametrizations
  expect_equal(
    distrib_gradient(do, y, th_o)[["mu1"]],
    distrib_gradient(ds, y, th_s)[["mu1"]]
  )
})

test_that("the precision form has correct derivatives of its own", {
  do <- mvgaussian2_distrib(2, parameters7::log_cholesky(2))
  th <- list(mu1 = 0.2, mu2 = -0.1, omega_log_L1 = 0.15, omega_log_L2 = -0.05,
             omega_L2.1 = 0.3)
  set.seed(26)
  y <- distrib_rng(do, 40, th)
  v0 <- unlist(th)
  ll <- function(v) {
    sum(distrib_pdf(do, y, as.list(stats::setNames(v, do@params)), log = TRUE))
  }

  g <- vapply(distrib_gradient(do, y, th), sum, numeric(1))
  expect_equal(unname(g), fd_grad(ll, v0), tolerance = 1e-6)

  h <- distrib_hessian(do, y, th)
  pr <- hess_pairs(do@params)
  for (nm in hess_names(do@params)) {
    k <- match(pr[[nm]], do@params)
    fd <- fd_hess(ll, v0, k[1], k[2])
    expect_equal(sum(h[[nm]]), fd, tolerance = 1e-4, label = nm)
  }
})

test_that("fit_distrib recovers the closed-form maximum likelihood estimate", {
  d <- mvgaussian1_distrib(2)
  true <- list(mu1 = 1.5, mu2 = -0.5, sigma_log_L1 = log(1.2), sigma_log_L2 = log(0.8),
               sigma_L2.1 = 0.6)
  set.seed(27)
  y <- distrib_rng(d, 2000, true)

  fit <- fit_distrib(d, y)
  expect_true(fit@converged)
  expect_identical(fit@n, 2000L)

  # The maximum likelihood estimator of a gaussian is the sample mean and the
  # sample second moment about it, both in closed form, so the maximum of the
  # log-likelihood is known without optimizing anything.
  mu_hat <- colMeans(y)
  s_hat <- crossprod(sweep(y, 2L, mu_hat)) / nrow(y)
  th_hat <- as.list(stats::setNames(
    c(mu_hat, parameters7::param_free(d@param, s_hat)), d@params
  ))
  ll_hat <- sum(distrib_pdf(d, y, th_hat, log = TRUE))

  # What a fit promises is the value of the objective, and that is what is
  # asked of it here. The point it stops AT is a different question: the
  # objective is flat near the maximum, so the last steps move the point much
  # more than they move the value, and how far they get is a fact about the
  # platform's arithmetic. The comparison on the value is the stable one.
  est <- coef(fit)
  expect_equal(as.numeric(logLik(fit)), ll_hat, tolerance = 1e-8)
  expect_equal(unname(mv_location(d, est)), unname(mu_hat), tolerance = 1e-4)
  expect_equal(unname(mv_sigma(d, est)), unname(s_hat), tolerance = 1e-4)

  sc <- vapply(distrib_gradient(d, y, as.list(est)), sum, numeric(1))
  expect_lt(max(abs(sc)) / nrow(y), 1e-4)

  # and the reported log-likelihood is the one the density gives there
  expect_equal(
    as.numeric(logLik(fit)),
    sum(distrib_pdf(d, y, as.list(est), log = TRUE))
  )
})

test_that("the three fitting methods reach the same optimum", {
  d <- mvgaussian1_distrib(2)
  set.seed(28)
  y <- distrib_rng(d, 800, list(mu1 = 0, mu2 = 1, sigma_log_L1 = 0, sigma_log_L2 = 0,
                                sigma_L2.1 = 0.4))
  lls <- vapply(c("fisher", "newton", "bfgs"), function(m) {
    as.numeric(logLik(fit_distrib(d, y, method = m)))
  }, numeric(1))
  expect_equal(diff(range(lls)), 0, tolerance = 1e-6)
})

test_that("a diagonal covariance is fitted with fewer parameters", {
  dd <- mvgaussian1_distrib(3, parameters7::diagonal_matrix(3))
  set.seed(29)
  y <- distrib_rng(dd, 1500, list(mu1 = 0, mu2 = 1, mu3 = -1,
                                  sigma_log_d1 = log(2),
                                  sigma_log_d2 = log(0.5),
                                  sigma_log_d3 = 0))
  fit <- fit_distrib(dd, y)

  expect_true(fit@converged)
  expect_identical(length(coef(fit)), 6L)
  s <- mv_sigma(dd, coef(fit))
  # the fitted covariance is diagonal by construction, not by luck
  expect_equal(max(abs(s[upper.tri(s)])), 0)
  expect_equal(unname(diag(s)), unname(apply(y, 2L, function(z) mean((z - mean(z))^2))),
    tolerance = 1e-4
  )
})

test_that("the fit prints the quantities a reader reads, and each once", {
  d <- mvgaussian1_distrib(2)
  set.seed(30)
  y <- distrib_rng(d, 400, list(mu1 = 0, mu2 = 0, sigma_log_L1 = 0, sigma_log_L2 = 0,
                                sigma_L2.1 = 0.3))
  fit <- fit_distrib(d, y)

  out <- utils::capture.output(print(fit))
  expect_true(any(grepl("^Location:", out)))
  expect_true(any(grepl("^Standard deviations:", out)))
  expect_true(any(grepl("^Correlations:", out)))

  # The mean and the covariance used to be printed a second time as bare
  # matrices, next to the same numbers with their standard errors. One table
  # per quantity.
  expect_false(any(grepl("^Mean:", out)))
  expect_false(any(grepl("^Covariance:", out)))

  # and the coordinates of the structure are not a parameter table anybody
  # needs to see
  expect_false(any(grepl("^Parameter scale:", out)))
  expect_false(any(grepl("^sigma_log_L1", out)))

  # Every link is the identity, so the link-scale table would repeat the
  # numbers above line for line.
  expect_true(any(grepl("link scale is the parameter scale", out)))
  expect_false(any(grepl("^Link scale", out)))

  # what the optimizer did is information, and is kept
  expect_true(any(grepl("^Method: .*iterations:", out)))
  expect_true(any(grepl("^Converged: yes \\(", out)))
})

test_that("check_distrib runs the multivariate battery", {
  set.seed(31)
  d <- mvgaussian1_distrib(2)
  th <- list(mu1 = 0.4, mu2 = -0.2, sigma_log_L1 = 0.1, sigma_log_L2 = -0.1, sigma_L2.1 = 0.35)
  res <- check_distrib(d, theta = th, nsim = 5e4, verbose = FALSE)

  expect_true(all(res$status == "OK"))
  # the battery is the multivariate one, not the univariate one run anyway
  expect_true("density integrates to 1" %in% res$check)
  expect_true("information matches the score variance" %in% res$check)
  expect_false("quantile/cdf round-trip" %in% res$check)
})

test_that("check_distrib catches a deliberately wrong component", {
  # The battery must be able to see a 5% error, or its agreements prove
  # nothing. Injected into the score, which several of the checks read.
  Wrong <- S7::new_class("WrongMv", parent = MvGaussianDistrib, package = NULL)
  gen <- distrib_gradient
  S7::method(gen, Wrong) <- function(distrib, y, theta,
                                     scale = c("parameter", "link"), ...) {
    g <- S7::method(distrib_gradient, MvGaussianDistrib)(distrib, y, theta)
    g[["mu1"]] <- 1.05 * g[["mu1"]]
    g
  }
  good <- mvgaussian1_distrib(2)
  bad <- Wrong(
    distrib_name = "wrong", dimension = "multivariate", n_dim = good@n_dim,
    bounds = good@bounds, params = good@params,
    params_interpretation = good@params_interpretation,
    n_params = good@n_params, params_bounds = good@params_bounds,
    link_params = good@link_params, param = good@param,
    inverted = good@inverted
  )

  th <- list(mu1 = 0.4, mu2 = -0.2, sigma_log_L1 = 0.1, sigma_log_L2 = -0.1, sigma_L2.1 = 0.35)
  set.seed(32)
  res_bad <- check_distrib(bad, theta = th, nsim = 5e4, verbose = FALSE)
  expect_identical(
    res_bad$status[res_bad$check == "gradient vs finite differences"], "FAIL"
  )
  set.seed(32)
  res_good <- check_distrib(good, theta = th, nsim = 5e4, verbose = FALSE)
  expect_identical(
    res_good$status[res_good$check == "gradient vs finite differences"], "OK"
  )
})

test_that("n_obs counts observations rather than entries", {
  expect_identical(n_obs(gaussian1_distrib(), c(1, 2, 3)), 3L)
  expect_identical(n_obs(mvgaussian1_distrib(2), matrix(0, 5, 2)), 5L)
  # a bare vector is one observation
  expect_identical(n_obs(mvgaussian1_distrib(3), c(1, 2, 3)), 1L)
  expect_identical(n_obs(mvgaussian1_distrib(2), matrix(numeric(0), 0, 2)), 0L)
})


test_that("closed-form third and fourth derivatives match one stencil", {
  for (inv in c(FALSE, TRUE)) {
    d <- if (inv) mvgaussian2_distrib(2, parameters7::log_cholesky(2))
         else mvgaussian1_distrib(2)
    set.seed(2)
    th <- generate_random_theta(d)
    y <- distrib_rng(d, 6, th)
    a3 <- distrib_deriv3(d, y, th)
    n3 <- numerical_deriv3(d, y, th)
    expect_setequal(names(a3), deriv_names(d@params, 3))
    for (k in names(a3)) {
      expect_equal(a3[[k]], n3[[k]], tolerance = 1e-5,
        label = paste(if (inv) "omega" else "sigma", "d3", k))
    }
    a4 <- distrib_deriv4(d, y, th)
    n4 <- numerical_deriv4(d, y, th)
    for (k in names(a4)) {
      expect_equal(a4[[k]], n4[[k]], tolerance = 1e-3,
        label = paste(if (inv) "omega" else "sigma", "d4", k))
    }
  }
})


test_that("the mixed response-parameter block agrees with numDeriv", {
  skip_if_not_installed("numDeriv")
  # the reference differentiates the ANALYTIC response gradient in theta, so
  # the two routes share no arithmetic; p starts at 2 and runs to 4 because a
  # single matrix coordinate says nothing about the cross terms
  set.seed(3)
  for (p in 2:4) {
    for (inv in c(FALSE, TRUE)) {
      d <- if (inv) mvgaussian2_distrib(p, parameters7::log_cholesky(p))
           else mvgaussian1_distrib(p)
      nm <- d@params
      v <- c(stats::rnorm(p, 0, 0.5), stats::rnorm(length(nm) - p, 0, 0.3))
      th <- stats::setNames(as.list(v), nm)
      y <- matrix(stats::rnorm(3 * p), ncol = p)
      got <- distrib_cross_y(d, y, th)
      expect_named(got, nm)
      for (k in seq_along(nm)) {
        ref <- matrix(numDeriv::jacobian(function(z) {
          vv <- v; vv[k] <- z
          as.vector(distrib_grad_y(d, y, stats::setNames(as.list(vv), nm)))
        }, v[k]), nrow = nrow(y), ncol = p)
        expect_equal(got[[k]], ref, tolerance = 1e-6,
                     info = paste("p =", p, "inverted =", inv, nm[k]))
      }
    }
  }
})


test_that("the mean block of the mixed derivative is Sigma^-1 and constant", {
  d <- mvgaussian1_distrib(3)
  th <- stats::setNames(as.list(c(0.5, -1, 2, 0.1, -0.2, 0.3, 0.4, 0.1, -0.1)),
                        d@params)
  y <- matrix(stats::rnorm(12), ncol = 3)
  got <- distrib_cross_y(d, y, th)
  si <- solve(mv_sigma(d, th))
  for (j in 1:3) {
    expect_equal(unname(got[[j]]),
                 matrix(unname(si)[j, ], nrow = 4, ncol = 3, byrow = TRUE))
  }
})


test_that("the families that do not implement it still refuse", {
  dd <- dirichlet_distrib(3)
  th <- stats::setNames(as.list(c(0, 0, 5)), dd@params)
  expect_error(distrib_cross_y(dd, matrix(c(0.2, 0.3, 0.5), ncol = 3), th),
               "closed form")
})


test_that("the higher mixed response derivatives agree with one difference", {
  skip_if_not_installed("numDeriv")
  # ONE difference on an ANALYTIC quantity, never two in a row: the second
  # theta derivative of the response Hessian is checked against a difference
  # of distrib_cross2_y, which is closed form, and the second theta
  # derivative of the response gradient against one of distrib_cross_y. A
  # nested reference reports gaps of 0.3 on correct code.
  set.seed(11)
  for (p in 2:3) {
    d <- mvgaussian1_distrib(p)
    nm <- d@params
    v <- c(stats::rnorm(p, 0, 0.4), stats::rnorm(length(nm) - p, 0, 0.25))
    th <- stats::setNames(as.list(v), nm)
    y <- matrix(stats::rnorm(3 * p), ncol = p)
    set_th <- function(z) stats::setNames(as.list(z), nm)

    c2 <- distrib_cross2_y(d, y, th)
    expect_named(c2, nm)
    for (k in seq_along(nm)) {
      ref <- matrix(numDeriv::jacobian(function(z) {
        vv <- v; vv[k] <- z
        as.vector(distrib_hess_y(d, y, set_th(vv)))
      }, v[k]), p, p)
      expect_equal(c2[[k]], ref, tolerance = 1e-6,
                   info = paste("cross2_y p =", p, nm[k]))
    }
    # the mean does not enter the response Hessian at all
    for (j in seq_len(p)) expect_equal(max(abs(c2[[j]])), 0)

    hh <- distrib_hess_y_hess(d, y, th)
    gh <- distrib_grad_y_hess(d, y, th)
    expect_named(hh, hess_names(nm))
    expect_named(gh, hess_names(nm))
    for (k in hess_names(nm)) {
      ij <- hess_pairs(nm)[[k]]
      a <- match(ij[1L], nm)
      b <- match(ij[2L], nm)
      ref2 <- matrix(numDeriv::jacobian(function(z) {
        vv <- v; vv[b] <- z
        as.vector(distrib_cross2_y(d, y, set_th(vv))[[a]])
      }, v[b]), p, p)
      expect_equal(hh[[k]], ref2, tolerance = 1e-5,
                   info = paste("hess_y_hess p =", p, k))
      ref3 <- matrix(numDeriv::jacobian(function(z) {
        vv <- v; vv[b] <- z
        as.vector(distrib_cross_y(d, y, set_th(vv))[[a]])
      }, v[b]), nrow(y), p)
      expect_equal(gh[[k]], ref3, tolerance = 1e-5,
                   info = paste("grad_y_hess p =", p, k))
    }
  }
})


test_that("the higher mixed derivatives do not depend on the order", {
  # a component collects the same terms whichever index is differentiated
  # first, which a tolerance cannot see from one direction alone
  skip_if_not_installed("numDeriv")
  set.seed(12)
  p <- 3
  d <- mvgaussian1_distrib(p)
  nm <- d@params
  v <- c(stats::rnorm(p, 0, 0.4), stats::rnorm(length(nm) - p, 0, 0.25))
  th <- stats::setNames(as.list(v), nm)
  y <- matrix(stats::rnorm(3 * p), ncol = p)
  hh <- distrib_hess_y_hess(d, y, th)
  for (k in hess_names(nm)) {
    ij <- hess_pairs(nm)[[k]]
    a <- match(ij[1L], nm)
    b <- match(ij[2L], nm)
    if (a == b) next
    alt <- matrix(numDeriv::jacobian(function(z) {
      vv <- v; vv[a] <- z
      as.vector(distrib_cross2_y(d, y, stats::setNames(as.list(vv), nm))[[b]])
    }, v[a]), p, p)
    expect_equal(hh[[k]], alt, tolerance = 1e-5, info = k)
  }
})


test_that("the two families are two classes under one parent", {
  d1 <- mvgaussian1_distrib(2)
  d2 <- mvgaussian2_distrib(2)
  expect_true(S7::S7_inherits(d1, MvGaussian1Distrib))
  expect_true(S7::S7_inherits(d2, MvGaussian2Distrib))
  expect_true(S7::S7_inherits(d1, MvGaussianDistrib))
  expect_true(S7::S7_inherits(d2, MvGaussianDistrib))
  expect_false(S7::S7_inherits(d1, MvGaussian2Distrib))
  expect_false(S7::S7_inherits(d2, MvGaussian1Distrib))
})

test_that("each family says which matrix its free values describe", {
  d1 <- mvgaussian1_distrib(2)
  d2 <- mvgaussian2_distrib(2)
  expect_false(d1@inverted)
  expect_true(d2@inverted)
  expect_identical(d1@params[3:5],
                   c("sigma_log_L1", "sigma_log_L2", "sigma_L2.1"))
  expect_identical(d2@params[3:5],
                   c("omega_log_L1", "omega_log_L2", "omega_L2.1"))
  expect_identical(unname(d1@params_interpretation[3]), "covariance")
  expect_identical(unname(d2@params_interpretation[3]), "precision")
})

test_that("a chart closed under inversion gives the two families one law", {
  # log_cholesky is closed, so the same data reach the same maximum from
  # either side and the split is a change of coordinates there
  set.seed(11)
  p <- 3
  n <- 300
  s <- unclass(parameters7::param_value(parameters7::ar1(p), c(0, atanh(0.5))))
  y <- matrix(stats::rnorm(n * p), n, p) %*% chol(s)
  y <- sweep(y, 2, c(1, -1, 0.5), "+")

  f1 <- fit_distrib(mvgaussian1_distrib(p), y)
  f2 <- fit_distrib(mvgaussian2_distrib(p), y)
  expect_equal(as.numeric(logLik(f1)), as.numeric(logLik(f2)),
               tolerance = 1e-6)
  # and the fitted matrices agree, one being the other's inverse
  expect_equal(mv_sigma(f1@distrib, as.list(coef(f1))), mv_sigma(f2@distrib, as.list(coef(f2))),
               tolerance = 1e-5)
})

test_that("a chart NOT closed under inversion gives two different models", {
  # THE NEGATIVE CONTROL: the inverse of an AR(1) covariance is tridiagonal
  # and is not an AR(1) at any parameters, so imposing the pattern on the two
  # sides is two models and the maximized likelihoods differ. Without this the
  # split would be cosmetic.
  set.seed(11)
  p <- 4
  n <- 400
  s <- unclass(parameters7::param_value(parameters7::ar1(p), c(0, atanh(0.7))))
  y <- matrix(stats::rnorm(n * p), n, p) %*% chol(s)
  y <- sweep(y, 2, c(1, -1, 0.5, 0), "+")

  a <- fit_distrib(mvgaussian1_distrib(p, parameters7::ar1(p)), y)
  b <- fit_distrib(mvgaussian2_distrib(p, parameters7::ar1(p)), y)
  expect_gt(abs(as.numeric(logLik(a)) - as.numeric(logLik(b))), 10)
  # the covariance side is the one the data were drawn from, so it wins
  expect_gt(as.numeric(logLik(a)), as.numeric(logLik(b)))
})

test_that("ar1_inv on the precision side is the AR(1) process again", {
  # the same law as an AR(1) covariance, written on the other side
  set.seed(3)
  p <- 4
  n <- 300
  s <- unclass(parameters7::param_value(parameters7::ar1(p), c(0.2, atanh(0.6))))
  y <- matrix(stats::rnorm(n * p), n, p) %*% chol(s)

  a <- fit_distrib(mvgaussian1_distrib(p, parameters7::ar1(p)), y)
  b <- fit_distrib(mvgaussian2_distrib(p, parameters7::ar1_inv(p)), y)
  expect_equal(as.numeric(logLik(a)), as.numeric(logLik(b)), tolerance = 1e-5)
})

test_that("neither family comments on the chart it is given", {
  # a pattern on the side it is not named for is a legitimate law, so it is
  # built in silence: no error, no warning, no message
  expect_silent(mvgaussian2_distrib(4, parameters7::ar1(4)))
  expect_silent(mvgaussian1_distrib(4, parameters7::ar1_inv(4)))
})

test_that("check_distrib passes on both families", {
  for (d in list(mvgaussian1_distrib(2), mvgaussian2_distrib(2),
                 mvgaussian1_distrib(3, parameters7::ar1(3)),
                 mvgaussian2_distrib(3, parameters7::ar1_inv(3)))) {
    r <- check_distrib(d, verbose = FALSE)
    expect_identical(sum(r$status == "FAIL"), 0L, info = d@distrib_name)
  }
})

test_that("the marginal keeps the parametrization it came from", {
  d1 <- mvgaussian1_distrib(3)
  d2 <- mvgaussian2_distrib(3)
  th1 <- as.list(stats::setNames(c(0, 1, -1, 0, 0, 0, 0.4, -0.2, 0.1),
                                 d1@params))
  th2 <- as.list(stats::setNames(unlist(th1), d2@params))

  m1 <- mv_marginal(d1, th1, c(1, 3))
  m2 <- mv_marginal(d2, th2, c(1, 3))
  expect_true(S7::S7_inherits(m1$distrib, MvGaussian1Distrib))
  expect_true(S7::S7_inherits(m2$distrib, MvGaussian2Distrib))

  # and both report the right submatrix of the covariance
  expect_equal(mv_sigma(m1$distrib, m1$theta),
               mv_sigma(d1, th1)[c(1, 3), c(1, 3)], ignore_attr = TRUE)
  expect_equal(mv_sigma(m2$distrib, m2$theta),
               mv_sigma(d2, th2)[c(1, 3), c(1, 3)], ignore_attr = TRUE)
})
