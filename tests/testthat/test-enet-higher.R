# The elastic net's third and fourth derivatives.
#
# They are assembled in the two rates (mu, a, c), where the log-density is
# quadratic in the response and linear in each rate, and carried onto
# (mu, lambda, alpha) by the bilinear map. Two things are checked: that the
# same assembly reproduces the hand-written score and Hessian, which is what
# licenses it at the orders where nothing else exists, and that orders three
# and four agree with ONE stencil applied to the analytic log-density --
# never a difference of a difference.

# one product stencil of the requested order, applied to the exact
# log-density; a repeated index uses the matching higher one-dimensional
# factor, distinct indices each contribute a central two-point factor
stencil_ld <- function(d, y, theta, idx, h) {
  fac <- list(
    "1" = list(o = c(-1, 1), w = c(-0.5, 0.5)),
    "2" = list(o = c(-1, 0, 1), w = c(1, -2, 1)),
    "3" = list(o = c(-2, -1, 1, 2), w = c(-0.5, 1, -1, 0.5)),
    "4" = list(o = c(-2, -1, 0, 1, 2), w = c(1, -4, 6, -4, 1))
  )
  tb <- table(idx)
  ks <- names(tb)
  fs <- lapply(as.integer(tb), function(m) fac[[as.character(m)]])
  grid <- expand.grid(lapply(fs, function(f) seq_along(f$o)))
  acc <- 0
  for (r in seq_len(nrow(grid))) {
    th <- theta
    w <- 1
    for (j in seq_along(ks)) {
      pick <- grid[r, j]
      th[[ks[j]]] <- th[[ks[j]]] + fs[[j]]$o[pick] * h
      w <- w * fs[[j]]$w[pick]
    }
    acc <- acc + w * distrib_pdf(d, y, th, log = TRUE)
  }
  acc / h^length(idx)
}

test_that("the third and fourth derivatives are the family's own", {
  d <- enet_distrib()
  for (g in list(distrib_deriv3, distrib_deriv4)) {
    m <- S7::method(g, S7::S7_class(d))
    owner <- attr(attr(m, "signature")[[1L]], "name")
    # the base class would mean the numerical fallback, which is what this
    # family had until the derivatives were written
    expect_identical(owner, "EnetDistrib")
  }
})

test_that("the assembly reproduces the hand-written score and Hessian", {
  # The same code path produces orders one to four, so agreement here at
  # machine precision is the licence for the orders that have nothing to be
  # compared against.
  d <- enet_distrib()
  y <- c(-2.3, -0.4, 0.15, 1.1, 3.7)
  for (th in list(list(mu = 0.2, lambda = 1.5, alpha = 0.6),
                  list(mu = -1.0, lambda = 0.4, alpha = 0.05),
                  list(mu = 2.0, lambda = 8.0, alpha = 0.97))) {
    g <- distributions7:::.enet_chain(y, th, 1L)
    h <- distributions7:::.enet_chain(y, th, 2L)
    expect_equal(g, distrib_gradient(d, y, th)[names(g)], tolerance = 1e-13)
    expect_equal(h, distrib_hessian(d, y, th)[names(h)], tolerance = 1e-12)
  }
})

test_that("orders three and four agree with one stencil on the log-density", {
  d <- enet_distrib()
  y <- 1.35
  nms <- c("mu", "lambda", "alpha")
  for (th in list(list(mu = 0.2, lambda = 1.5, alpha = 0.6),
                  list(mu = -1.0, lambda = 0.4, alpha = 0.2))) {
    for (ord in 3:4) {
      got <- distrib_deriv3(d, y, th)
      if (ord == 4L) got <- distrib_deriv4(d, y, th)
      idx <- deriv_indices(nms, ord)
      nm <- deriv_names(nms, ord)
      # The step is measured, not chosen. Swept over 1e-3 to 1e-1, the
      # third-order stencil is best at 1e-3 (1.4e-6) and the fourth at 3e-3
      # (8.3e-5): below that the rounding, which grows as h^-4, dominates,
      # and above it the h^2 truncation does. The step must also keep
      # alpha inside (0, 1), which 2h of 3e-3 does at these values.
      h <- if (ord == 3L) 1e-3 else 3e-3
      tol <- if (ord == 3L) 1e-4 else 1e-3
      ref <- vapply(idx, function(I) stencil_ld(d, y, th, nms[I], h),
                    numeric(1))
      sc <- max(abs(ref), 1)
      expect_lt(max(abs(unlist(got[nm]) - ref)) / sc, tol)
    }
  }
})

test_that("a derivative in mu alone vanishes past the second", {
  # the log-density is quadratic in y - mu away from the kink, so the third
  # and fourth derivatives in the location are exactly zero rather than
  # small, and the Hessian is minus the Gaussian rate
  d <- enet_distrib()
  th <- list(mu = 0.4, lambda = 2, alpha = 0.35)
  y <- c(-1, 0.9, 2.2)
  expect_equal(distrib_deriv3(d, y, th)$mu_mu_mu, rep(0, 3))
  expect_equal(distrib_deriv4(d, y, th)$mu_mu_mu_mu, rep(0, 3))
  expect_equal(distrib_hessian(d, y, th)$mu_mu,
               rep(-th$lambda * (1 - th$alpha), 3))
})

test_that("the expected derivatives go through the approximation", {
  d <- enet_distrib()
  th <- list(mu = 0, lambda = 1.5, alpha = 0.6)
  e <- distrib_deriv3(d, 0.5, th, expected = TRUE, approx = "bartlett")
  expect_type(e, "list")
  expect_named(e, deriv_names(c("mu", "lambda", "alpha"), 3L))
})
