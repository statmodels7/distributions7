# Regression tests for the expected information and its derivatives of orders
# 0 to 2 at the points where the routes before distributions7 0.88.0 failed:
# negbin2 and negbin1 at large theta, negbin1 at large mu/theta, lognormal2 at
# a small var/mean^2, and the skew t at a small nu. Each check compares the
# Richardson-extrapolated central difference of order k with the analytic order
# k + 1, relative to the largest component of the same pair; before the
# repairs these read up to 1e7 (negbin2), 1e2 (negbin1), 1e-3 (skew t).

d34_chain_error <- function(d, th, k) {
  P <- d@params
  pn <- length(P)
  gen <- list(distrib_expected_hessian, distrib_dexpected_hessian,
              distrib_d2expected_hessian, distrib_d3expected_hessian)
  y <- distrib_rng(d, 1, th)
  up <- gen[[k + 2]](d, y, th)
  pairs <- which(upper.tri(diag(pn), diag = TRUE), arr.ind = TRUE)
  tuples <- if (k == 0) list(integer(0)) else if (k == 1) as.list(seq_len(pn)) else
    lapply(seq_len(nrow(pairs)), function(r) pairs[r, ])
  hpn <- getFromNamespace("hess_pair_name", "distributions7")
  err <- numeric(0); scl <- numeric(0); grp <- character(0)
  for (cc in seq_len(pn)) {
    h <- 1e-3 * abs(th[[cc]])
    ev <- function(s) { t2 <- th; t2[[cc]] <- th[[cc]] + s; gen[[k + 1]](d, y, t2) }
    p1 <- ev(h); m1 <- ev(-h); p2 <- ev(h / 2); m2 <- ev(-h / 2)
    for (r in seq_len(nrow(pairs))) {
      a <- pairs[r, 1]; b <- pairs[r, 2]
      for (tu in tuples) {
        kk <- switch(k + 1, hpn(P, a, b), dexpected_key(P, a, b, tu[1]),
                     d2expected_key(P, a, b, tu[1], tu[2]))
        ku <- switch(k + 1, dexpected_key(P, a, b, cc),
                     d2expected_key(P, a, b, min(tu[1], cc), max(tu[1], cc)),
                     d3expected_key(P, a, b, tu[1], tu[2], cc))
        D1 <- (p1[[kk]][1] - m1[[kk]][1]) / (2 * h)
        D2 <- (p2[[kk]][1] - m2[[kk]][1]) / h
        err <- c(err, abs((4 * D2 - D1) / 3 - up[[ku]][1]))
        scl <- c(scl, abs(up[[ku]][1]))
        grp <- c(grp, hpn(P, a, b))
      }
    }
  }
  ps <- tapply(scl, grp, max)[grp]
  max(err / pmax(ps, 1e-6 * max(scl)))
}

test_that("negbin2's expected information is consistent across orders at large theta", {
  d <- negbin2_distrib()
  for (th in list(list(mu = 1e6, theta = 1e6), list(mu = 5e6, theta = 1e6),
                  list(mu = 3e4, theta = 1e5))) {
    for (k in 0:2) expect_lt(d34_chain_error(d, th, k), 1e-8,
                             label = paste("negbin2", th$mu, th$theta, k))
  }
})

test_that("negbin1's expected information is consistent across orders at large theta and mu/theta", {
  d <- negbin1_distrib()
  for (th in list(list(mu = 3e5, theta = 1e6), list(mu = 1e7, theta = 1e4),
                  list(mu = 2e3, theta = 1))) {
    for (k in 0:2) expect_lt(d34_chain_error(d, th, k), 1e-8,
                             label = paste("negbin1", th$mu, th$theta, k))
  }
})

test_that("lognormal2's derivatives of the expected information hold at a small var/mean^2", {
  # mpmath at 50 digits: E = J' diag(-1/S, -1/(2 S^2)) J, S = log1p(v/m^2)
  d <- lognormal2_distrib()
  th <- list(mean = 40, var = 0.001)
  d1 <- distrib_dexpected_hessian(d, 1, th)
  d2 <- distrib_d2expected_hessian(d, 1, th)
  expect_equal(d1[["mean_mean_mean"]], 1.40624918619848e-4, tolerance = 1e-8)
  expect_equal(d2[["mean_mean_mean_mean"]], -1.05468648274838e-5, tolerance = 1e-8)
  for (k in 0:2) expect_lt(d34_chain_error(d, th, k), 1e-8)
})

test_that("the skew t's expected information holds at a small nu", {
  d <- skewt_distrib()
  for (nu in c(0.1, 0.15)) {
    th <- list(mu = 0.3, sigma = 1.7, alpha = -0.9, nu = nu)
    for (k in 0:2) expect_lt(d34_chain_error(d, th, k), 1e-8,
                             label = paste("skewt nu", nu, "order", k))
  }
})
