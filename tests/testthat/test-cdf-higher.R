# Third and fourth derivatives of the distribution function.
#
# Two routes, as at the orders below: a discrete family sums the identity
# exactly and a continuous one differences its analytic cdf. Both are checked
# against something that shares no code with them -- the discrete against
# differences of the cdf, since checking it against the partial-expectation
# sum would be the same sum twice, and the continuous against the
# partial-expectation integral.

# one product stencil of the requested order on log F, the independent
# reference for the discrete route
fd_log_cdf <- function(d, q, theta, I, h) {
  nms <- d@params
  tb <- table(nms[I])
  ks <- names(tb)
  fac <- list(list(o = c(-1, 1), w = c(-0.5, 0.5)),
              list(o = c(-1, 0, 1), w = c(1, -2, 1)),
              list(o = c(-2, -1, 1, 2), w = c(-0.5, 1, -1, 0.5)),
              list(o = c(-2, -1, 0, 1, 2), w = c(1, -4, 6, -4, 1)))
  fs <- fac[as.integer(tb)]
  grid <- expand.grid(lapply(fs, function(f) seq_along(f$o)))
  acc <- 0
  for (r in seq_len(nrow(grid))) {
    t2 <- theta
    w <- 1
    for (j in seq_along(ks)) {
      pick <- grid[r, j]
      t2[[ks[j]]] <- theta[[ks[j]]] + fs[[j]]$o[pick] * h
      w <- w * fs[[j]]$w[pick]
    }
    acc <- acc + w * distrib_cdf(d, q, t2, log.p = TRUE)
  }
  acc / h^length(I)
}

test_that("the general conversion reproduces the written-out orders 1 and 2", {
  # cdf_scale_k() is the moment-to-cumulant relation at any order and
  # cdf_tail_scale() is the second order written out. They must agree
  # exactly, on both tails and both scales, or the new orders rest on a
  # different convention from the old ones.
  for (cs in list(list(gaussian1_distrib(), list(mu = 0.3, sigma = 1.4), 0.8),
                  list(poisson_distrib(), list(mu = 3), 4),
                  list(gamma1_distrib(), list(mu = 2, phi = 0.6), 1.7))) {
    d <- cs[[1]]; th <- cs[[2]]; q <- cs[[3]]
    Fq <- distrib_cdf(d, q, th)
    tb <- distributions7:::cdf_tables(d, q, th, 2L)
    for (lt in c(TRUE, FALSE)) {
      for (lg in c(TRUE, FALSE)) {
        g1 <- distributions7:::cdf_tail_scale(d, Fq, tb[[1]], NULL, lt, lg)
        g2 <- distributions7:::cdf_scale_k(d, Fq, tb, 1L, lt, lg)
        expect_equal(g2[names(g1)], g1)
        h1 <- distributions7:::cdf_tail_scale(d, Fq, tb[[1]], tb[[2]], lt, lg)
        h2 <- distributions7:::cdf_scale_k(d, Fq, tb, 2L, lt, lg)
        expect_equal(h2[names(h1)], h1)
      }
    }
  }
})

test_that("the general discrete sum reproduces the written-out one", {
  for (cs in list(list(poisson_distrib(), list(mu = 3), 4),
                  list(negbin2_distrib(), list(mu = 3, theta = 2), 5))) {
    d <- cs[[1]]; th <- cs[[2]]; q <- cs[[3]]
    for (k in 1:2) {
      a <- distributions7:::discrete_cdf_deriv(d, q, th, k)
      b <- distributions7:::discrete_cdf_deriv_k(d, q, th, k)
      expect_equal(b[names(a)], a)
    }
  }
})

test_that("the discrete exact sum agrees with differences of the cdf", {
  for (cs in list(list(poisson_distrib(), list(mu = 3), 4),
                  list(negbin2_distrib(), list(mu = 3, theta = 2), 5),
                  list(binomial_distrib(size = 10), list(mu = 0.4), 4))) {
    d <- cs[[1]]; th <- cs[[2]]; q <- cs[[3]]
    for (ord in 3:4) {
      got <- if (ord == 3L) distrib_deriv3_cdf(d, q, th) else
                            distrib_deriv4_cdf(d, q, th)
      idx <- deriv_indices(d@params, ord)
      # the step is the one the order can carry: the rounding of a fourth
      # difference grows as h^-4
      h <- if (ord == 3L) 1e-3 else 5e-3
      ref <- vapply(idx, function(I) fd_log_cdf(d, q, th, I, h), numeric(1))
      expect_lt(max(abs(unlist(got) - ref)) / max(abs(ref), 1), 1e-3)
    }
  }
})

test_that("the continuous stencil agrees with the partial expectation", {
  # d^I F(q) = integral over y <= q of f(y) (d^I f / f)(y), which shares no
  # code with differencing the cdf
  part_exp <- function(d, q, th, I) {
    nms <- d@params
    f <- function(y) {
      ell <- distributions7:::parent_ell(d, y, th, length(I), nms)
      distrib_pdf(d, y, th) * distributions7:::bell_f_ratio(nms[I], ell)
    }
    stats::integrate(f, max(d@bounds[1], q - 60), q,
                     rel.tol = 1e-10, subdivisions = 500L)$value
  }
  for (cs in list(list(gaussian1_distrib(), list(mu = 0.3, sigma = 1.4), 0.8),
                  list(logistic_distrib(), list(mu = 0.2, sigma = 1.1), 0.9))) {
    d <- cs[[1]]; th <- cs[[2]]; q <- cs[[3]]
    for (ord in 3:4) {
      got <- distributions7:::cdf_tables(d, q, th, ord)[[ord]]
      idx <- deriv_indices(d@params, ord)
      ref <- vapply(idx, function(I) part_exp(d, q, th, I), numeric(1))
      expect_lt(max(abs(unlist(got) - ref)) / max(abs(ref), 1), 1e-4)
    }
  }
})

test_that("both tails and both scales are served", {
  d <- gaussian1_distrib()
  th <- list(mu = 0, sigma = 1)
  q <- 0.6
  nm <- deriv_names(d@params, 3L)
  for (ord in 3:4) {
    f <- if (ord == 3L) distrib_deriv3_cdf else distrib_deriv4_cdf
    for (lt in c(TRUE, FALSE)) {
      for (lg in c(TRUE, FALSE)) {
        v <- f(d, q, th, lower.tail = lt, log = lg)
        expect_named(v, deriv_names(d@params, ord))
        expect_true(all(vapply(v, is.finite, logical(1))))
      }
    }
    # on the natural scale the upper tail is minus the lower, S = 1 - F
    lo <- f(d, q, th, lower.tail = TRUE, log = FALSE)
    up <- f(d, q, th, lower.tail = FALSE, log = FALSE)
    expect_equal(unlist(up), -unlist(lo))
  }
})
