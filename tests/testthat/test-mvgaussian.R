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
  expect_error(mvgaussian_distrib(0), "positive integer")
  expect_error(mvgaussian_distrib(2.5), "positive integer")
  expect_error(
    mvgaussian_distrib(2,
      sigma = parameters7::log_cholesky(2),
      omega = parameters7::log_cholesky(2)
    ),
    "at most one"
  )
  expect_error(mvgaussian_distrib(2, sigma = diag(2)), "parameter")
  expect_error(
    mvgaussian_distrib(3, sigma = parameters7::log_cholesky(2)),
    "dimension 2 but the distribution has dimension 3"
  )

  # A rank-deficient structure is a penalty, not a density, and the refusal is
  # the same whichever side it parametrises.
  pen <- parameters7::scaled_matrix(crossprod(diff(diag(5), differences = 2)))
  expect_error(mvgaussian_distrib(5, sigma = pen), "rank deficient")
  expect_error(mvgaussian_distrib(5, omega = pen), "rank deficient")
})

test_that("the parameters are the mean components and the structure's free values", {
  d <- mvgaussian_distrib(2)

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
  dd <- mvgaussian_distrib(3, sigma = parameters7::diagonal_matrix(3))
  expect_identical(dd@params, c("mu1", "mu2", "mu3", "sigma_d1", "sigma_d2", "sigma_d3"))
})

test_that("the density is the formula, written out by hand", {
  d <- mvgaussian_distrib(2)
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
  d <- mvgaussian_distrib(2)
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
  d <- mvgaussian_distrib(2)
  th <- list(mu1 = c(0, 1), mu2 = 0, sigma_log_L1 = 0, sigma_log_L2 = 0, sigma_L2.1 = 0)
  expect_error(distrib_pdf(d, rbind(c(0, 0), c(1, 1)), th), "belong to a model")
})

test_that("the score and the Hessian agree with finite differences", {
  d <- mvgaussian_distrib(2)
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
  d <- mvgaussian_distrib(2)
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
  d <- mvgaussian_distrib(2)
  th <- list(mu1 = 1.5, mu2 = -0.5, sigma_log_L1 = 0.2, sigma_log_L2 = -0.1, sigma_L2.1 = 0.6)
  set.seed(23)
  r <- distrib_rng(d, 2e5, th)

  expect_identical(dim(r), c(200000L, 2L))
  expect_equal(unname(colMeans(r)), unname(mv_location(d, th)), tolerance = 0.02)
  expect_equal(unname(stats::cov(r)), unname(mv_sigma(d, th)), tolerance = 0.02)
})

test_that("the response derivatives are closed form", {
  d <- mvgaussian_distrib(2)
  th <- list(mu1 = 0.5, mu2 = -0.3, sigma_log_L1 = 0.1, sigma_log_L2 = -0.2, sigma_L2.1 = 0.4)
  y <- rbind(c(0, 0), c(1, -1))
  si <- solve(mv_sigma(d, th))

  gy <- distrib_grad_y(d, y, th)
  want <- -sweep(y, 2L, mv_location(d, th)) %*% si
  expect_equal(unname(gy), unname(want))
  expect_equal(unname(distrib_hess_y(d, y, th)), unname(-si))
})

test_that("the link scale is the parameter scale", {
  d <- mvgaussian_distrib(2)
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
  d <- mvgaussian_distrib(2)
  th <- list(mu1 = 0, mu2 = 0, sigma_log_L1 = 0, sigma_log_L2 = 0, sigma_L2.1 = 0)

  expect_error(distrib_cdf(d, rbind(c(0, 0)), th), "orthant")
  expect_error(distrib_quantile(d, 0.5, th), "ordering")
})

test_that("the precision form describes the same law as the covariance form", {
  # Omega = Sigma^{-1} through a structure on the other side: the density, the
  # score and the Hessian must agree once the two are matched.
  ds <- mvgaussian_distrib(2)
  do <- mvgaussian_distrib(2, omega = parameters7::log_cholesky(2))

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
  # the score in the mean is the same quantity in both parametrisations
  expect_equal(
    distrib_gradient(do, y, th_o)[["mu1"]],
    distrib_gradient(ds, y, th_s)[["mu1"]]
  )
})

test_that("the precision form has correct derivatives of its own", {
  do <- mvgaussian_distrib(2, omega = parameters7::log_cholesky(2))
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
  d <- mvgaussian_distrib(2)
  true <- list(mu1 = 1.5, mu2 = -0.5, sigma_log_L1 = log(1.2), sigma_log_L2 = log(0.8),
               sigma_L2.1 = 0.6)
  set.seed(27)
  y <- distrib_rng(d, 2000, true)

  fit <- fit_distrib(d, y)
  expect_true(fit@converged)
  expect_identical(fit@n, 2000L)

  # The maximum likelihood estimator of a gaussian is the sample mean and the
  # sample second moment about it, both in closed form, so the maximum of the
  # log-likelihood is known without optimising anything.
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
  d <- mvgaussian_distrib(2)
  set.seed(28)
  y <- distrib_rng(d, 800, list(mu1 = 0, mu2 = 1, sigma_log_L1 = 0, sigma_log_L2 = 0,
                                sigma_L2.1 = 0.4))
  lls <- vapply(c("fisher", "newton", "bfgs"), function(m) {
    as.numeric(logLik(fit_distrib(d, y, method = m)))
  }, numeric(1))
  expect_equal(diff(range(lls)), 0, tolerance = 1e-6)
})

test_that("a diagonal covariance is fitted with fewer parameters", {
  dd <- mvgaussian_distrib(3, sigma = parameters7::diagonal_matrix(3))
  set.seed(29)
  y <- distrib_rng(dd, 1500, list(mu1 = 0, mu2 = 1, mu3 = -1,
                                  sigma_d1 = log(2), sigma_d2 = log(0.5),
                                  sigma_d3 = 0))
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
  d <- mvgaussian_distrib(2)
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

  # what the optimiser did is information, and is kept
  expect_true(any(grepl("^Method: .*iterations:", out)))
  expect_true(any(grepl("^Converged: yes \\(", out)))
})

test_that("check_distrib runs the multivariate battery", {
  set.seed(31)
  d <- mvgaussian_distrib(2)
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
  good <- mvgaussian_distrib(2)
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
  expect_identical(n_obs(gaussian_distrib(), c(1, 2, 3)), 3L)
  expect_identical(n_obs(mvgaussian_distrib(2), matrix(0, 5, 2)), 5L)
  # a bare vector is one observation
  expect_identical(n_obs(mvgaussian_distrib(3), c(1, 2, 3)), 1L)
  expect_identical(n_obs(mvgaussian_distrib(2), matrix(numeric(0), 0, 2)), 0L)
})


test_that("closed-form third and fourth derivatives match one stencil", {
  for (inv in c(FALSE, TRUE)) {
    d <- if (inv) mvgaussian_distrib(2, omega = parameters7::log_cholesky(2))
         else mvgaussian_distrib(2)
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
