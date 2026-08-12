# The multivariate Student t. Every closed form is checked against a route the
# implementation does not take: the density against the formula written out by
# hand, the derivatives against finite differences, the moments against the
# generator, and the whole family against the gaussian it tends to.

fd_grad <- function(f, x, h = 1e-5) {
  vapply(seq_along(x), function(k) {
    up <- dn <- x
    up[k] <- x[k] + h
    dn[k] <- x[k] - h
    (f(up) - f(dn)) / (2 * h)
  }, numeric(1))
}

# A second derivative from ONE stencil, for the reason given in
# test-mvgaussian.R: composing two first differences makes the reference the
# error of an error.
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

th2 <- function() {
  list(mu1 = 0.4, mu2 = -0.3, sigma_log_L1 = 0.1, sigma_log_L2 = -0.2, sigma_L2.1 = 0.5, nu = 6)
}

test_that("the constructor validates its arguments", {
  expect_error(mvstudent_t_distrib(0), "positive integer")
  expect_error(mvstudent_t_distrib(2, sigma = diag(2)), "parameter")
  expect_error(
    mvstudent_t_distrib(3, sigma = parameters7::log_cholesky(2)),
    "dimension 2 but the distribution has dimension 3"
  )
  expect_error(mvstudent_t_distrib(2, link_nu = "log"), "link object")

  pen <- parameters7::scaled_matrix(crossprod(diff(diag(5), differences = 2)))
  expect_error(mvstudent_t_distrib(5, sigma = pen), "rank deficient")
})

test_that("the degrees of freedom are a parameter of their own", {
  d <- mvstudent_t_distrib(2)

  expect_identical(d@params, c("mu1", "mu2", "sigma_log_L1", "sigma_log_L2", "sigma_L2.1", "nu"))
  expect_identical(d@n_params, 6L)
  expect_true(S7::S7_inherits(d, multivariate_distrib))

  # nu is positive and carries a log link, so this family's link scale is NOT
  # its parameter scale -- unlike the multivariate gaussian, whose every link
  # is the identity.
  expect_identical(d@link_params[["nu"]]@link_name, "log")
  expect_identical(d@params_bounds[["nu"]], c(0, Inf))
  expect_true(all(vapply(
    d@link_params[d@params != "nu"],
    function(l) l@link_name == "identity", logical(1)
  )))
})

test_that("the density is the formula, written out by hand", {
  d <- mvstudent_t_distrib(2)
  th <- th2()
  y <- rbind(c(0, 0), c(1, -1), c(-0.5, 0.8), c(3, 2))

  mu <- mv_location(d, th)
  s <- mv_sigma(d, th)
  want <- apply(y, 1, function(row) {
    r <- row - mu
    q <- drop(t(r) %*% solve(s) %*% r)
    lgamma((6 + 2) / 2) - lgamma(6 / 2) - log(6 * pi) - 0.5 * log(det(s)) -
      ((6 + 2) / 2) * log1p(q / 6)
  })

  expect_equal(distrib_pdf(d, y, th, log = TRUE), unname(want))
  expect_equal(distrib_pdf(d, y, th), unname(exp(want)))
  expect_equal(distrib_pdf(d, c(0, 0), th, log = TRUE), unname(want[1]))
  expect_error(distrib_pdf(d, c(0, 0, 0), th), "read as one observation")
})

test_that("the scale matrix and the covariance are different objects", {
  d <- mvstudent_t_distrib(2)
  th <- th2()
  s <- mv_sigma(d, th)

  # what the parametrization carries
  l <- matrix(c(exp(0.1), 0, 0.5, exp(-0.2)), 2, 2, byrow = TRUE)
  expect_equal(unname(s), tcrossprod(l))
  # and the moment, inflated by nu/(nu-2)
  expect_equal(variance(d, th), s * 6 / 4)
  expect_equal(mean(d, th), mv_location(d, th))

  # below two degrees of freedom the covariance does not exist while the
  # density does, and below one the mean goes with it
  th$nu <- 1.5
  expect_true(is.finite(distrib_pdf(d, c(0, 0), th)))
  expect_true(all(is.infinite(variance(d, th))))
  expect_false(anyNA(mean(d, th)))

  th$nu <- 0.5
  expect_true(all(is.nan(mean(d, th))))
  expect_true(all(is.infinite(variance(d, th))))
  # and the scale matrix is unaffected: it is a parameter, not a moment
  expect_equal(unname(mv_sigma(d, th)), tcrossprod(l))
})

test_that("the score and the Hessian agree with finite differences", {
  d <- mvstudent_t_distrib(2)
  th <- th2()
  set.seed(41)
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
    expect_equal(sum(h[[nm]]), fd_hess(ll, v0, k[1], k[2]),
      tolerance = 1e-4, label = nm
    )
  }
})

test_that("the link scale is the chain rule through the log link on nu", {
  d <- mvstudent_t_distrib(2)
  th <- th2()
  set.seed(42)
  y <- distrib_rng(d, 30, th)

  gp <- distrib_gradient(d, y, th)
  gl <- distrib_gradient(d, y, th, scale = "link")
  # only nu is transformed; dnu/deta = nu for a log link
  expect_equal(gl[["mu1"]], gp[["mu1"]])
  expect_equal(gl[["nu"]], gp[["nu"]] * 6)

  # and the same score, obtained by differentiating the log-likelihood written
  # as a function of the unconstrained parameters
  eta0 <- c(unlist(th)[1:5], log(6))
  lle <- function(v) {
    tv <- as.list(stats::setNames(c(v[1:5], exp(v[6])), d@params))
    sum(distrib_pdf(d, y, tv, log = TRUE))
  }
  expect_equal(
    unname(vapply(gl, sum, numeric(1))), fd_grad(lle, eta0),
    tolerance = 1e-5
  )
})

test_that("the generator matches the moments and the Mahalanobis law", {
  d <- mvstudent_t_distrib(2)
  th <- th2()
  set.seed(43)
  r <- distrib_rng(d, 2e5, th)

  expect_identical(dim(r), c(200000L, 2L))
  expect_equal(unname(colMeans(r)), unname(mv_location(d, th)), tolerance = 0.05)
  expect_equal(unname(stats::cov(r)), unname(variance(d, th)), tolerance = 0.05)

  # The moments alone would not separate this from a gaussian with the same
  # covariance. The distribution of the Mahalanobis distance does: q/p is
  # F(p, nu) for an elliptical t, and F(2, 6) is not what a gaussian gives.
  q <- stats::mahalanobis(r, mv_location(d, th), mv_sigma(d, th))
  expect_gt(suppressWarnings(stats::ks.test(q / 2, "pf", 2, 6)$p.value), 0.01)
  expect_lt(suppressWarnings(stats::ks.test(q, "pchisq", 2)$p.value), 1e-10)
})

test_that("the response derivatives are closed form", {
  d <- mvstudent_t_distrib(2)
  th <- th2()
  y <- rbind(c(0, 0), c(1, -1), c(2.5, -2))

  # -(nu + p)/(nu + q) * Sigma^{-1} (y - mu): the gaussian score scaled by the
  # weight that makes the family resistant.
  si <- solve(mv_sigma(d, th))
  r <- sweep(y, 2L, mv_location(d, th))
  q <- rowSums((r %*% si) * r)
  want <- -((6 + 2) / (6 + q)) * (r %*% si)
  expect_equal(unname(distrib_grad_y(d, y, th)), unname(want))

  gy_num <- t(apply(y, 1, function(row) {
    fd_grad(function(z) distrib_pdf(d, matrix(z, 1L), th, log = TRUE), row)
  }))
  expect_equal(unname(distrib_grad_y(d, y, th)), gy_num, tolerance = 1e-6)
})

test_that("the gaussian is the limit of large degrees of freedom", {
  d <- mvstudent_t_distrib(2)
  g <- mvgaussian_distrib(2)
  th <- th2()
  th$nu <- 1e8
  y <- rbind(c(0, 0), c(1, -1), c(-0.5, 0.8), c(3, 2))

  expect_equal(
    distrib_pdf(d, y, th, log = TRUE),
    distrib_pdf(g, y, th[1:5], log = TRUE),
    tolerance = 1e-5
  )
  expect_equal(
    distrib_grad_y(d, y, th), distrib_grad_y(g, y, th[1:5]),
    tolerance = 1e-5
  )
})

test_that("a marginal is a t with the same degrees of freedom", {
  d <- mvstudent_t_distrib(3)
  th <- as.list(stats::setNames(
    c(0, 1, -1, 0.1, -0.1, 0.2, 0.5, -0.2, 0.3, 7), d@params
  ))

  m <- mv_marginal(d, th, c(1L, 3L))
  expect_true(S7::S7_inherits(m$distrib, MvStudentTDistrib))
  expect_identical(m$distrib@n_dim, 2L)
  # The scale matrix is the corresponding block and nu does not change: a
  # marginal of a t is not a t with fewer degrees of freedom.
  expect_equal(
    unname(mv_sigma(m$distrib, m$theta)),
    unname(mv_sigma(d, th)[c(1, 3), c(1, 3)])
  )
  expect_equal(m$theta$nu, 7)
  expect_equal(unname(unlist(m$theta[c("mu1", "mu2")])), c(0, -1))

  expect_error(mv_marginal(d, th, 4L), "distinct coordinates")
})

test_that("the expected information is exact and needs no strategy", {
  d <- mvstudent_t_distrib(2)
  th <- th2()
  y1 <- matrix(c(0, 0), 1L, 2L)

  # It was sampled until the scale mixture closed it. The check that it is
  # the information and not something else: the average of the OBSERVED
  # Hessian over a large sample is the same matrix, the first Bartlett
  # identity holding for this family away from nothing.
  eh <- distrib_expected_hessian(d, y1, th)
  set.seed(45)
  big <- distrib_rng(d, 1e5, th)
  hb <- distrib_hessian(d, big, th)
  for (nm in hess_names(d@params)) {
    expect_equal(eh[[nm]][1], mean(hb[[nm]]), tolerance = 0.08, label = nm)
  }

  # the answer does not move between calls, which a sampled one does
  expect_identical(distrib_expected_hessian(d, y1, th),
                   distrib_expected_hessian(d, y1, th))
  # and `approx` is now read by nothing here, so naming one is harmless
  # rather than an error: what refuses it is fisher_scoring(), where the
  # argument would silently do nothing
  expect_no_error(distrib_expected_hessian(d, y1, th, approx = "integrate"))
})

test_that("check_distrib runs the multivariate battery on the t", {
  set.seed(46)
  d <- mvstudent_t_distrib(2)
  res <- check_distrib(d, theta = th2(), nsim = 5e4, verbose = FALSE)
  expect_true(all(res$status == "OK"))
})

test_that("fit_distrib recovers the degrees of freedom", {
  d <- mvstudent_t_distrib(2)
  true <- list(mu1 = 1, mu2 = -0.5, sigma_log_L1 = log(1.2), sigma_log_L2 = log(0.8),
               sigma_L2.1 = 0.6, nu = 6)
  set.seed(47)
  y <- distrib_rng(d, 3000, true)

  fit <- fit_distrib(d, y, method = "newton")
  expect_true(fit@converged)
  est <- coef(fit)

  # nu is the parameter the family exists for, and it is the hardest to
  # estimate: the likelihood is flat in it, so the tolerance is wide.
  expect_equal(unname(est[["nu"]]), 6, tolerance = 0.3)
  expect_equal(unname(mv_location(d, est)), c(1, -0.5), tolerance = 0.1)

  # the optimum is a stationary point on the link scale, which is where the
  # optimization happens
  sc <- vapply(
    distrib_gradient(d, y, as.list(est), scale = "link"), sum, numeric(1)
  )
  expect_lt(max(abs(sc)) / nrow(y), 1e-4)
  expect_equal(
    as.numeric(logLik(fit)),
    sum(distrib_pdf(d, y, as.list(est), log = TRUE))
  )

  # a gaussian fitted to the same data is beaten by the t, which is the whole
  # point of the family
  g <- mvgaussian_distrib(2)
  expect_gt(as.numeric(logLik(fit)), as.numeric(logLik(fit_distrib(g, y))))
})

test_that("Fisher scoring inverts the exact information here", {
  d <- mvstudent_t_distrib(2)
  set.seed(48)
  y <- distrib_rng(d, 800, list(mu1 = 0, mu2 = 0, sigma_log_L1 = 0, sigma_log_L2 = 0,
                                sigma_L2.1 = 0.3, nu = 5))
  ff <- fit_distrib(d, y, method = fisher_scoring())
  fn <- fit_distrib(d, y, method = "newton")

  expect_true(ff@converged)
  # the two now agree far more closely than they could when one of them
  # inverted a sampled matrix
  expect_equal(as.numeric(logLik(ff)), as.numeric(logLik(fn)),
               tolerance = 1e-8)

  # naming an approximation is refused, since the family has a closed form
  # and the argument would be ignored
  expect_error(fit_distrib(d, y, method = fisher_scoring(approx = "mc")),
               "would be ignored")
})


test_that("the response Hessian is the reweighted gaussian expression", {
  d <- mvstudent_t_distrib(2)
  th <- list(mu1 = 0.4, mu2 = -0.3, sigma_log_L1 = 0.2, sigma_log_L2 = -0.1,
             sigma_L2.1 = 0.5, nu = 5)
  set.seed(3)
  y <- distrib_rng(d, 6, th)

  # one central stencil on the closed-form response gradient, per coordinate
  H <- distrib_hess_y(d, y, th)
  expect_equal(dim(H), c(2L, 2L, 6L))
  h <- 1e-6
  for (j in 1:2) {
    yp <- y; ym <- y
    yp[, j] <- y[, j] + h
    ym[, j] <- y[, j] - h
    fd <- (distrib_grad_y(d, yp, th) - distrib_grad_y(d, ym, th)) / (2 * h)
    for (i in seq_len(nrow(y))) {
      expect_equal(H[, j, i], fd[i, ], tolerance = 1e-6)
    }
  }

  # and in the nu -> infinity limit it is the gaussian's -Sigma^{-1}
  th_big <- th; th_big$nu <- 1e8
  Hb <- distrib_hess_y(d, y, th_big)
  Hg <- distrib_hess_y(mvgaussian_distrib(2), y, th[names(th) != "nu"])
  expect_equal(Hb[, , 1], Hg, tolerance = 1e-6, ignore_attr = TRUE)
})


test_that("response derivatives without a closed form are refused, not guessed", {
  # The base class refuses rather than serving the univariate fallback, whose
  # stencils difference along a line and would return numbers of the wrong
  # shape for a matrix response.
  d <- mvgaussian_distrib(2)
  th <- list(mu1 = 0, mu2 = 0, sigma_log_L1 = 0, sigma_log_L2 = 0,
             sigma_L2.1 = 0.3)
  y <- distrib_rng(d, 3, th)
  expect_error(distrib_cross_y(d, y, th), "not defined")
  expect_error(
    distrib_cross_y(mvstudent_t_distrib(2), cbind(y, 0)[, 1:2], c(th, nu = 5)),
    "not defined"
  )
})


test_that("third and fourth parameter derivatives work on a matrix response", {
  # The fallbacks stencil along the parameters, so a matrix response passes
  # through untouched; checked against a Richardson jacobian of the summed
  # analytic Hessian, which shares nothing with the single stencil.
  d <- mvgaussian_distrib(2)
  th <- list(mu1 = 0.3, mu2 = -0.2, sigma_log_L1 = 0.1, sigma_log_L2 = -0.1,
             sigma_L2.1 = 0.4)
  set.seed(4)
  y <- distrib_rng(d, 4, th)

  n3 <- distrib_deriv3(d, y, th)
  expect_setequal(names(n3), deriv_names(d@params, 3))
  expect_true(all(lengths(n3) == nrow(y)))

  hess_sum <- function(v) {
    th2 <- as.list(stats::setNames(v, d@params))
    vapply(distrib_hessian(d, y, th2), sum, numeric(1))
  }
  J <- numDeriv::jacobian(hess_sum, unlist(th))
  hn <- hess_names(d@params)
  idx <- deriv_indices(d@params, 3)
  for (t in seq_along(n3)) {
    hcomp <- paste(d@params[idx[[t]][1:2]], collapse = "_")
    expect_equal(sum(n3[[t]]), J[match(hcomp, hn), idx[[t]][3]],
      tolerance = 1e-6, label = names(n3)[t]
    )
  }

  n4 <- distrib_deriv4(d, y, th)
  expect_setequal(names(n4), deriv_names(d@params, 4))
  expect_true(all(lengths(n4) == nrow(y)))
  expect_true(all(is.finite(unlist(n4))))
})

test_that("the expected information is closed and matches the sampling route", {
  # The family's own scale mixture closes it: q/(q+nu) is exactly
  # Beta(p/2, nu/2) and is independent of the direction, so every
  # expectation is a Beta moment or a polygamma. The reference is the
  # definition -- minus the outer product of the score -- estimated from
  # draws, which is the route the closed form replaces and shares no
  # arithmetic with it.
  skip_on_cran()
  set.seed(7)
  for (p in c(2L, 3L)) {
    d <- mvstudent_t_distrib(p)
    y0 <- matrix(stats::rnorm(20 * p), 20, p)
    th <- distrib_start(d, y0, 1L)[[1L]]
    for (nu in c(4, 8, 25)) {
      th$nu <- nu
      cl <- distrib_expected_hessian(d, y0[1, , drop = FALSE], th)
      ys <- distrib_rng(d, 4e5, th)
      s <- distrib_gradient(d, ys, th)
      idx <- mv_hess_indices(d)
      nm <- hess_names(d@params)
      got <- vapply(nm, function(k) cl[[k]][1L], numeric(1))
      mc <- vapply(idx, function(ab) -mean(s[[ab[1L]]] * s[[ab[2L]]]),
                   numeric(1))
      # ON THE SCALE OF THE WHOLE MATRIX: several blocks are exactly zero
      # by symmetry, and a relative measure against a Monte Carlo estimate
      # of zero reports one whatever the closed form says
      expect_lt(max(abs(got - mc)) / max(abs(mc)), 0.02,
                label = sprintf("p = %d, nu = %g", p, nu))
    }
  }
})

test_that("the expected information tends to the gaussian's", {
  p <- 3L
  d <- mvstudent_t_distrib(p)
  g <- mvgaussian_distrib(p)
  set.seed(3)
  y0 <- matrix(stats::rnorm(10 * p), 10, p)
  tht <- distrib_start(d, y0, 1L)[[1L]]
  tht$nu <- 1e7
  thg <- tht[setdiff(names(tht), "nu")]
  et <- distrib_expected_hessian(d, y0[1, , drop = FALSE], tht)
  eg <- distrib_expected_hessian(g, y0[1, , drop = FALSE], thg)
  common <- intersect(names(et), names(eg))
  expect_gt(length(common), 40L)
  expect_equal(vapply(common, function(k) et[[k]][1L], numeric(1)),
               vapply(common, function(k) eg[[k]][1L], numeric(1)),
               tolerance = 1e-5)
  # and the degrees of freedom cease to be identified in the limit
  expect_lt(abs(et$nu_nu[1L]), 1e-12)
})

test_that("the location stays orthogonal to the scale and the shape", {
  # odd in z against even in z: every cross-expectation with the location
  # vanishes, which is why those blocks are exactly zero rather than small
  d <- mvstudent_t_distrib(2)
  set.seed(5)
  y0 <- matrix(stats::rnorm(20), 10, 2)
  th <- distrib_start(d, y0, 1L)[[1L]]
  e <- distrib_expected_hessian(d, y0, th)
  cross <- grep("^mu[0-9]+_(sigma|nu)", names(e), value = TRUE)
  expect_gt(length(cross), 3L)
  for (k in cross) expect_identical(unique(e[[k]]), 0)
})
