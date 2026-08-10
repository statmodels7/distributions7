# The two shapes the census left that need no new algebra: a family written
# out against a tabulated map, and a family with no location.

t2_ref <- function(d, y, th, order, pr) {
  inner <- if (order == 1L) distrib_cross_y else distrib_cross2_y
  vapply(seq_along(y), function(i) {
    numDeriv::grad(function(v) {
      t2 <- th
      t2[[pr[1L]]] <- v
      inner(d, y[i], t2)[[pr[2L]]]
    }, th[[pr[1L]]])
  }, numeric(1))
}

t2_close <- function(got, ref, tol = 1e-6) {
  expect_lt(max(abs(got - ref)) / max(1e-3, max(abs(ref))), tol)
}

test_that("a written-out reparametrization carries the grid through its map", {
  skip_if_not_installed("numDeriv")
  cases <- list(
    list(d = gaussian2_distrib(), th = list(mu = 0.4, sigma2 = 1.7),
         y = c(-1.6, 0.3, 1.9)),
    list(d = gaussian3_distrib(), th = list(mu = 0.4, tau = 0.6),
         y = c(-1.6, 0.3, 1.9)),
    list(d = laplace2_distrib(), th = list(mu = 0.3, lambda = 1.4),
         y = c(-1.6, 0.35, 1.9))
  )
  for (cs in cases) {
    d <- cs$d
    th <- cs$th
    y <- cs$y
    # one theta-derivative of the analytic response curvature
    got1 <- distrib_cross2_y(d, y, th)
    for (q in d@params) {
      ref <- vapply(seq_along(y), function(i) {
        numDeriv::grad(function(v) {
          t2 <- th
          t2[[q]] <- v
          distrib_hess_y(d, y[i], t2)
        }, th[[q]])
      }, numeric(1))
      t2_close(got1[[q]], ref)
    }
    for (order in 1:2) {
      got <- if (order == 1L) distrib_grad_y_hess(d, y, th) else
        distrib_hess_y_hess(d, y, th)
      prs <- distributions7:::hess_pairs(d@params)
      for (nm in names(prs)) {
        t2_close(got[[nm]], t2_ref(d, y, th, order, prs[[nm]]))
      }
    }
  }
})

test_that("gaussian2 agrees with gaussian1 at the same law", {
  # the two describe the same distribution in different coordinates, so the
  # mapped components must equal the parent's carried by hand -- a check that
  # owes nothing to numDeriv
  y <- c(-1.1, 0.7, 2.0)
  s2 <- 1.7
  s <- sqrt(s2)
  d1 <- gaussian1_distrib()
  d2 <- gaussian2_distrib()
  th1 <- list(mu = 0.4, sigma = s)
  th2 <- list(mu = 0.4, sigma2 = s2)
  ds <- 0.5 / s
  d2s <- -0.25 / s2^1.5

  g1 <- distrib_grad_y_hess(d1, y, th1)
  c1 <- distrib_cross_y(d1, y, th1)
  g2 <- distrib_grad_y_hess(d2, y, th2)
  expect_equal(g2$mu_mu, g1$mu_mu, tolerance = 1e-12)
  expect_equal(g2$mu_sigma2, g1$mu_sigma * ds, tolerance = 1e-12)
  expect_equal(g2$sigma2_sigma2, g1$sigma_sigma * ds^2 + c1$sigma * d2s,
               tolerance = 1e-9)
})

test_that("a family with no location closes its scale", {
  skip_if_not_installed("numDeriv")
  y <- c(0.4, 1.3, 2.8)
  cases <- list(
    list(d = exponential_distrib(), th = list(mu = 1.6), at = "mu"),
    list(d = weibull1_distrib(), th = list(mu = 1.4, sigma = 1.8), at = "mu"),
    list(d = gpd_distrib(), th = list(sigma = 1.5, xi = 0.3), at = "sigma")
  )
  for (cs in cases) {
    d <- cs$d
    th <- cs$th
    nm <- paste(cs$at, cs$at, sep = "_")
    for (order in 1:2) {
      got <- if (order == 1L) distrib_grad_y_hess(d, y, th) else
        distrib_hess_y_hess(d, y, th)
      expect_named(got, hess_names(d@params))
      t2_close(got[[nm]], t2_ref(d, y, th, order, c(cs$at, cs$at)))
    }
    # and the first-order one in the scale
    ref <- vapply(seq_along(y), function(i) {
      numDeriv::grad(function(v) {
        t2 <- th
        t2[[cs$at]] <- v
        distrib_hess_y(d, y[i], t2)
      }, th[[cs$at]])
    }, numeric(1))
    t2_close(distrib_cross2_y(d, y, th)[[cs$at]], ref)
  }
})

test_that("the exponential is closed outright, having only a scale", {
  # one parameter and it is the scale, so there is nothing left to difference
  skip_if_not_installed("numDeriv")
  y <- c(0.3, 1.1, 2.5)
  d <- exponential_distrib()
  th <- list(mu = 1.7)
  expect_named(distrib_grad_y_hess(d, y, th), "mu_mu")
  expect_named(distrib_cross2_y(d, y, th), "mu")
  # l = -log(mu) - y/mu, so l^(y) = -1/mu and l^(yy) = 0: the curvature does
  # not move at all, and its derivatives are exactly zero
  expect_equal(distrib_cross2_y(d, y, th)$mu, rep(0, 3), tolerance = 1e-14)
  expect_equal(distrib_hess_y_hess(d, y, th)$mu_mu, rep(0, 3),
               tolerance = 1e-14)
  # while the gradient's do not vanish: d/dmu of -1/mu is 1/mu^2 and the
  # second is -2/mu^3
  expect_equal(distrib_grad_y_hess(d, y, th)$mu_mu, rep(-2 / 1.7^3, 3),
               tolerance = 1e-12)
})

test_that("the lognormal is the gaussian at log y, carried by the Jacobian", {
  skip_if_not_installed("numDeriv")
  y <- c(0.4, 1.2, 2.7)
  d <- lognormal1_distrib()
  th <- list(mu = 0.3, sigma2 = 0.8)

  got1 <- distrib_cross2_y(d, y, th)
  for (q in d@params) {
    ref <- vapply(seq_along(y), function(i) {
      numDeriv::grad(function(v) {
        t2 <- th
        t2[[q]] <- v
        distrib_hess_y(d, y[i], t2)
      }, th[[q]])
    }, numeric(1))
    t2_close(got1[[q]], ref)
  }
  prs <- distributions7:::hess_pairs(d@params)
  for (order in 1:2) {
    got <- if (order == 1L) distrib_grad_y_hess(d, y, th) else
      distrib_hess_y_hess(d, y, th)
    for (nm in names(prs)) t2_close(got[[nm]], t2_ref(d, y, th, order,
                                                      prs[[nm]]))
  }
})

test_that("the lognormal's chain is the gaussian's carried by hand", {
  # a check that owes nothing to numDeriv: the transformation carries no
  # parameter, so every theta-derivative is the gaussian's at t = log y and
  # only the Jacobian of the response derivatives enters
  y <- c(0.4, 1.2, 2.7)
  t <- log(y)
  th <- list(mu = 0.3, sigma2 = 0.8)
  d <- lognormal1_distrib()
  g2 <- gaussian2_distrib()

  expect_equal(distrib_cross_y(d, y, th)$mu,
               distrib_cross_y(g2, t, th)$mu / y, tolerance = 1e-12)
  expect_equal(distrib_cross2_y(d, y, th)$mu,
               (distrib_cross2_y(g2, t, th)$mu -
                  distrib_cross_y(g2, t, th)$mu) / y^2, tolerance = 1e-12)
  expect_equal(distrib_grad_y_hess(d, y, th)$mu_sigma2,
               distrib_grad_y_hess(g2, t, th)$mu_sigma2 / y,
               tolerance = 1e-12)
  expect_equal(distrib_hess_y_hess(d, y, th)$sigma2_sigma2,
               (distrib_hess_y_hess(g2, t, th)$sigma2_sigma2 -
                  distrib_grad_y_hess(g2, t, th)$sigma2_sigma2) / y^2,
               tolerance = 1e-12)
})

test_that("lognormal2 inherits it through the map", {
  # lognormal2 is a reparametrize() wrapper on lognormal1, so closing the
  # parent is what makes the child exact rather than merely registered
  skip_if_not_installed("numDeriv")
  y <- c(0.5, 1.4, 3.1)
  d <- lognormal2_distrib()
  th <- list(mean = 1.5, var = 0.9)
  prs <- distributions7:::hess_pairs(d@params)
  for (order in 1:2) {
    got <- if (order == 1L) distrib_grad_y_hess(d, y, th) else
      distrib_hess_y_hess(d, y, th)
    for (nm in names(prs)) t2_close(got[[nm]], t2_ref(d, y, th, order,
                                                      prs[[nm]]))
  }
})
