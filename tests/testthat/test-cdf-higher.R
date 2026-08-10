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

test_that("the two routes to the truncation constant agree", {
  # d^B Z / Z can be had as F(u) - F(l^-) differentiated, or as a truncated
  # expectation of the same Bell quantity. They share no code, so their
  # agreement cross-validates both: the cdf derivatives of the parent on one
  # side and expectation()'s quadrature on the other.
  d <- gaussian1_distrib()
  th <- list(mu = 0, sigma = 1)
  tr <- truncated(d, lower = 0.2, upper = 4)
  params <- tr@params
  parent <- tr@parent_distrib
  Z <- distributions7:::trunc_constants(tr, th)$Z

  for (k in 1:4) {
    dZ <- distributions7:::trunc_mass_derivs(tr, th, k)
    expect_false(is.null(dZ))
    quad <- vapply(deriv_indices(params, k), function(I) {
      b <- params[I]
      expectation(tr, function(y, theta) {
        distributions7:::bell_f_ratio(
          b, distributions7:::parent_ell(parent, y, theta, length(b), params))
      }, th)
    }, numeric(1))
    expect_lt(max(abs(unlist(dZ) / Z - quad)) / max(abs(quad), 1), 1e-12)
  }
})

test_that("the gate refuses a stencil dressed as a closed form", {
  # the third and fourth defaults sit on `distrib` rather than on
  # `continuous_distrib`, so a gate that excluded only the latter would take
  # a stencil for a closed form and feed its noise into the truncation
  expect_true(all(vapply(1:4, function(k)
    distributions7:::has_exact_cdf_deriv(gaussian1_distrib(), k), logical(1))))
  expect_true(all(vapply(1:4, function(k)
    distributions7:::has_exact_cdf_deriv(poisson_distrib(), k), logical(1))))
  expect_false(any(vapply(1:4, function(k)
    distributions7:::has_exact_cdf_deriv(gamma1_distrib(), k), logical(1))))
})

test_that("the location-scale closed form reproduces the written-out orders", {
  for (cs in list(list(gaussian1_distrib(), list(mu = 0.3, sigma = 1.4), 0.9),
                  list(logistic_distrib(), list(mu = 0.2, sigma = 1.1), 0.7),
                  list(cauchy_distrib(), list(mu = 0.1, sigma = 1.3), 0.5),
                  list(laplace_distrib(), list(mu = 0.2, sigma = 1.1), 0.8))) {
    d <- cs[[1]]; th <- cs[[2]]; q <- cs[[3]]
    for (k in 1:2) {
      a <- distributions7:::loc_scale_cdf_deriv(d, q, th, k)
      b <- distributions7:::loc_scale_cdf_deriv_k(d, q, th, k)
      expect_equal(b[names(a)], a, tolerance = 1e-12)
    }
    # and at the new orders it agrees with the stencil route
    for (k in 3:4) {
      b <- unlist(distributions7:::loc_scale_cdf_deriv_k(d, q, th, k))
      r <- unlist(distributions7:::numerical_cdf_deriv_k(d, q, th, k))
      expect_lt(max(abs(b - r)) / max(abs(r), 1), 1e-4)
    }
  }
})


# --- the mapped and partially location-scale families -----------------------

cdfk_cases <- function() {
  list(
    list("gaussian2", gaussian2_distrib(), list(mu = 0.3, sigma2 = 1.4), 0.8, TRUE),
    list("gaussian3", gaussian3_distrib(), list(mu = 0.3, tau = 0.7), 0.8, TRUE),
    list("gumbel", gumbel_distrib(), list(mu = 0.2, sigma = 1.1), 0.9, TRUE),
    list("lognormal1", lognormal1_distrib(), list(mu = 0.3, sigma2 = 0.7), 1.4, TRUE),
    list("lognormal2", lognormal2_distrib(), list(mean = 1.5, var = 1.2), 1.4, TRUE),
    list("student_t1", student_t1_distrib(),
         list(mu = 0.2, sigma = 1.2, nu = 6), 0.7, FALSE),
    list("pseudohuber", pseudohuber_distrib(),
         list(mu = 0.2, sigma = 1.1, nu = 2), 0.7, FALSE),
    list("skewnormal1", skewnormal1_distrib(),
         list(mu = 0.2, sigma = 1.2, alpha = 1.3), 0.6, FALSE)
  )
}

test_that("the new routes agree with the partial expectation", {
  # the reference shares no code with the chain rule over a map, nor with the
  # location-scale construction, nor with the stencil
  part_exp <- function(d, q, th, I) {
    nms <- d@params
    f <- function(y) {
      ell <- distributions7:::parent_ell(d, y, th, length(I), nms)
      distrib_pdf(d, y, th) * distributions7:::bell_f_ratio(nms[I], ell)
    }
    stats::integrate(f, max(d@bounds[1] + 1e-10, q - 60), q,
                     rel.tol = 1e-11, subdivisions = 800L)$value
  }
  for (cs in cdfk_cases()) {
    nm <- cs[[1]]; d <- cs[[2]]; th <- cs[[3]]; q <- cs[[4]]; full <- cs[[5]]
    for (o in 3:4) {
      got <- if (o == 3L) distrib_deriv3_cdf(d, q, th, log = FALSE) else
                          distrib_deriv4_cdf(d, q, th, log = FALSE)
      ref <- vapply(deriv_indices(d@params, o),
                    function(I) part_exp(d, q, th, I), numeric(1))
      # a family closed throughout is held to machine precision; one whose
      # shape components still come from a stencil is held to the stencil's
      tol <- if (full) 1e-11 else 1e-4
      expect_lt(max(abs(unlist(got) - ref)) / max(1, max(abs(ref))), tol,
                label = sprintf("%s order %d", nm, o))
    }
  }
})

test_that("the gate refuses to carry a differenced parent", {
  # the gamma differences its own cdf at every order, the derivative of the
  # incomplete gamma in its shape having no elementary form, so a chain rule
  # over it must not report a closed form
  expect_false(any(vapply(1:4, function(k)
    distributions7:::has_exact_cdf_deriv(gamma1_distrib(), k), logical(1))))
  # while the gaussian and the Laplace are exact at every order, which is what
  # lets the mapped route close the lognormal and the second Laplace
  for (d in list(gaussian1_distrib(), laplace_distrib())) {
    expect_true(all(vapply(1:4, function(k)
      distributions7:::has_exact_cdf_deriv(d, k), logical(1))))
  }
  # and the weibull is exact now that the survival route serves it, which is
  # why weibull3 follows through the wrapper
  expect_true(all(vapply(1:4, function(k)
    distributions7:::has_exact_cdf_deriv(weibull1_distrib(), k), logical(1))))
})

test_that("both tails and both scales are served by the new routes", {
  for (cs in cdfk_cases()) {
    d <- cs[[2]]; th <- cs[[3]]; q <- cs[[4]]
    for (o in 3:4) {
      f <- if (o == 3L) distrib_deriv3_cdf else distrib_deriv4_cdf
      for (lt in c(TRUE, FALSE)) {
        for (lg in c(TRUE, FALSE)) {
          v <- f(d, q, th, lower.tail = lt, log = lg)
          expect_named(v, deriv_names(d@params, o))
          expect_true(all(vapply(v, is.finite, logical(1))), info = cs[[1]])
        }
      }
      lo <- f(d, q, th, lower.tail = TRUE, log = FALSE)
      up <- f(d, q, th, lower.tail = FALSE, log = FALSE)
      expect_equal(unlist(up), -unlist(lo), info = cs[[1]])
    }
  }
})


# --- the exponential-survival families, and the Laplace pair ----------------

# one stencil of the requested order on the ANALYTIC distribution function.
# This is the reference a non-regular family needs: the partial-expectation
# identity integrates the log-density's parameter derivatives, and the
# Laplace's second derivative in the location carries a point mass at the kink
# that no integral over the density can see. Measured, the shipped Laplace --
# validated long before this route existed -- disagrees with that integral by
# 1.3 relative at orders two to four while agreeing with the stencil below.
fd_on_cdf <- function(d, q, th, I, h) {
  nms <- d@params
  used <- sort(unique(I))
  mult <- tabulate(match(I, used), length(used))
  fac <- list(list(o = c(-1, 1), w = c(-0.5, 0.5)),
              list(o = c(-1, 0, 1), w = c(1, -2, 1)),
              list(o = c(-2, -1, 1, 2), w = c(-0.5, 1, -1, 0.5)),
              list(o = c(-2, -1, 0, 1, 2), w = c(1, -4, 6, -4, 1)))
  fs <- fac[mult]
  grid <- expand.grid(lapply(fs, function(x) seq_along(x$o)))
  acc <- 0
  for (r in seq_len(nrow(grid))) {
    t2 <- th
    w <- 1
    for (j in seq_along(used)) {
      pick <- grid[r, j]
      t2[[nms[used[j]]]] <- th[[nms[used[j]]]] + fs[[j]]$o[pick] * h
      w <- w * fs[[j]]$w[pick]
    }
    acc <- acc + w * distrib_cdf(d, q, t2)
  }
  acc / h^length(I)
}

test_that("the survival route and the Laplace map agree with a stencil on F", {
  hs <- c(1e-4, 1e-3, 3e-3, 1e-2)
  tols <- c(1e-6, 1e-5, 1e-3, 1e-2)
  for (cs in list(
        list("exponential", exponential_distrib(), list(mu = 2), 1.3),
        list("weibull1", weibull1_distrib(), list(mu = 2, sigma = 1.5), 1.7),
        list("weibull3", weibull3_distrib(), list(mean = 2, sigma = 1.5), 1.7),
        list("laplace2", laplace2_distrib(), list(mu = 0.3, lambda = 1.4), 0.9))) {
    nm <- cs[[1]]; d <- cs[[2]]; th <- cs[[3]]; q <- cs[[4]]
    for (o in 1:4) {
      got <- switch(o,
        distrib_grad_cdf(d, q, th, log = FALSE),
        distrib_hess_cdf(d, q, th, log = FALSE),
        distrib_deriv3_cdf(d, q, th, log = FALSE),
        distrib_deriv4_cdf(d, q, th, log = FALSE))
      ref <- vapply(deriv_indices(d@params, o),
                    function(I) fd_on_cdf(d, q, th, I, hs[o]), numeric(1))
      names(ref) <- deriv_names(d@params, o)
      expect_lt(max(abs(unlist(got)[names(ref)] - ref)) / max(1, max(abs(ref))),
                tols[o], label = sprintf("%s order %d", nm, o))
    }
  }
})

test_that("the upper tail of an exponential survival stays finite", {
  # log S is L, so its derivatives are L's own and need no division by
  # 1 - F, which is exactly one in double precision past q/mu = 37
  d <- exponential_distrib()
  th <- list(mu = 1)
  for (q in c(20, 40, 700)) {
    v <- distrib_grad_cdf(d, q, th, lower.tail = FALSE, log = TRUE)
    w <- distrib_hess_cdf(d, q, th, lower.tail = FALSE, log = TRUE)
    expect_equal(v[[1]], q / th$mu^2, tolerance = 1e-12)
    expect_equal(w[[1]], -2 * q / th$mu^3, tolerance = 1e-12)
    expect_true(all(is.finite(unlist(
      distrib_deriv4_cdf(d, q, th, lower.tail = FALSE, log = TRUE)))))
  }
  # and on the natural scale the two tails are still negatives of each other
  for (q in c(0.5, 2)) {
    lo <- distrib_deriv3_cdf(d, q, th, lower.tail = TRUE, log = FALSE)
    up <- distrib_deriv3_cdf(d, q, th, lower.tail = FALSE, log = FALSE)
    expect_equal(unlist(up), -unlist(lo))
  }
})

test_that("the partial-expectation reference is invalid at a kink", {
  # recorded as a test because it cost time: the identity holds where the
  # log-density is twice differentiable in the parameters, and the Laplace's
  # is not. The gaussian agrees with both references, the Laplace with only
  # one, and the code is the same in both cases.
  part_exp <- function(d, q, th, I) {
    nms <- d@params
    f <- function(y) {
      ell <- distributions7:::parent_ell(d, y, th, length(I), nms)
      distrib_pdf(d, y, th) * distributions7:::bell_f_ratio(nms[I], ell)
    }
    stats::integrate(f, q - 60, q, rel.tol = 1e-11, subdivisions = 800L)$value
  }
  q <- 0.9
  thg <- list(mu = 0.3, sigma = 0.7)
  gg <- unlist(distrib_hess_cdf(gaussian1_distrib(), q, thg, log = FALSE))
  rg <- vapply(deriv_indices(c("mu", "sigma"), 2L),
               function(I) part_exp(gaussian1_distrib(), q, thg, I), numeric(1))
  names(rg) <- deriv_names(c("mu", "sigma"), 2L)
  expect_lt(max(abs(gg[names(rg)] - rg)), 1e-12)

  gl <- unlist(distrib_hess_cdf(laplace_distrib(), q, thg, log = FALSE))
  rl <- vapply(deriv_indices(c("mu", "sigma"), 2L),
               function(I) part_exp(laplace_distrib(), q, thg, I), numeric(1))
  names(rl) <- deriv_names(c("mu", "sigma"), 2L)
  expect_gt(max(abs(gl[names(rl)] - rl)), 0.1)
  # while the stencil on F confirms the shipped value
  rs <- vapply(deriv_indices(c("mu", "sigma"), 2L),
               function(I) fd_on_cdf(laplace_distrib(), q, thg, I, 1e-3),
               numeric(1))
  names(rs) <- deriv_names(c("mu", "sigma"), 2L)
  expect_lt(max(abs(gl[names(rs)] - rs)), 1e-5)
})


# --- the generalized Pareto -------------------------------------------------

test_that("the generalized Pareto agrees with a stencil away from zero", {
  # the stencil is not usable near xi = 0: the step is 1e-3 to 3e-3, so at a
  # shape of 1e-6 the reference differentiates over a range a thousand times
  # the value itself. The exponential limit below is the check there.
  d <- gpd_distrib()
  hs <- c(1e-5, 1e-4, 1e-3, 3e-3)
  tols <- c(1e-8, 1e-6, 1e-4, 1e-3)
  q <- 1.3
  for (xi in c(-0.45, -0.3, -0.1, 0.1, 0.3, 0.5, 0.54, 1, 2)) {
    th <- list(sigma = 1.4, xi = xi)
    for (o in 1:4) {
      got <- switch(o,
        distrib_grad_cdf(d, q, th, log = FALSE),
        distrib_hess_cdf(d, q, th, log = FALSE),
        distrib_deriv3_cdf(d, q, th, log = FALSE),
        distrib_deriv4_cdf(d, q, th, log = FALSE))
      ref <- vapply(deriv_indices(d@params, o),
                    function(I) fd_on_cdf(d, q, th, I, hs[o]), numeric(1))
      names(ref) <- deriv_names(d@params, o)
      expect_lt(max(abs(unlist(got)[names(ref)] - ref)) / max(1, max(abs(ref))),
                tols[o], label = sprintf("xi = %g order %d", xi, o))
    }
  }
})

test_that("the exponential limit of the generalized Pareto is exact", {
  # at xi = 0 the scale components must be the exponential family's, and that
  # family reaches them through L = -q/mu, which shares no arithmetic with the
  # series for log1p(u)/u. This is the reference near zero, where a stencil in
  # the shape is not one.
  d <- gpd_distrib()
  e <- exponential_distrib()
  qq <- c(0.2, 1.3, 5, 40)
  sig <- 1.4
  for (o in 1:4) {
    a <- switch(o,
      distrib_grad_cdf(d, qq, list(sigma = sig, xi = 0), log = FALSE),
      distrib_hess_cdf(d, qq, list(sigma = sig, xi = 0), log = FALSE),
      distrib_deriv3_cdf(d, qq, list(sigma = sig, xi = 0), log = FALSE),
      distrib_deriv4_cdf(d, qq, list(sigma = sig, xi = 0), log = FALSE))
    b <- switch(o,
      distrib_grad_cdf(e, qq, list(mu = sig), log = FALSE),
      distrib_hess_cdf(e, qq, list(mu = sig), log = FALSE),
      distrib_deriv3_cdf(e, qq, list(mu = sig), log = FALSE),
      distrib_deriv4_cdf(e, qq, list(mu = sig), log = FALSE))
    expect_equal(a[[paste(rep("sigma", o), collapse = "_")]],
                 b[[paste(rep("mu", o), collapse = "_")]],
                 tolerance = 1e-12, info = sprintf("order %d", o))
  }
})

test_that("the generalized Pareto is continuous through a zero shape", {
  # the shape enters only through Lambda(u) = log1p(u)/u, which is analytic at
  # the origin, so every component moves linearly in xi rather than blowing up
  d <- gpd_distrib()
  qq <- c(0.2, 1.3, 5)
  base <- distrib_deriv4_cdf(d, qq, list(sigma = 1.4, xi = 0), log = FALSE)
  prev <- Inf
  for (xi in c(1e-7, 1e-9)) {
    v <- distrib_deriv4_cdf(d, qq, list(sigma = 1.4, xi = xi), log = FALSE)
    g <- max(vapply(names(base), function(k)
      max(abs(v[[k]] - base[[k]]) / pmax(1, abs(base[[k]]))), numeric(1)))
    expect_lt(g, 1e-5)
    # and it shrinks with the shape, which a cancelling form would not do
    expect_lt(g, prev)
    prev <- g
  }
  expect_true(all(is.finite(unlist(
    distrib_deriv4_cdf(d, qq, list(sigma = 1.4, xi = -1e-12), log = FALSE)))))
})

test_that("the series and the recursion agree at the switch", {
  # the branch changes at |u| = 0.5; a seam would show as a step here
  d <- gpd_distrib()
  th <- list(sigma = 1.4, xi = 1)
  # an EVEN grid: second differences of an uneven one measure the spacing
  qs <- seq(0.697, 0.703, by = 0.0015)
  v <- vapply(qs, function(q)
    distrib_deriv4_cdf(d, q, th, log = FALSE)[["xi_xi_xi_xi"]], numeric(1))
  # on an even grid the second differences of a smooth function are small next
  # to the first ones; a step at the seam would make them comparable
  d1 <- diff(v)
  d2 <- diff(d1)
  expect_lt(max(abs(d2)) / max(abs(d1)), 0.01)
})

test_that("a negative shape bounds the support above", {
  d <- gpd_distrib()
  th <- list(sigma = 1.4, xi = -0.5)
  # the upper endpoint is sigma/|xi| = 2.8
  v <- distrib_deriv3_cdf(d, c(1, 2.7, 2.9, 10), th, log = FALSE)
  expect_true(all(vapply(v, function(x) all(is.finite(x)), logical(1))))
  expect_true(all(vapply(v, function(x) all(x[3:4] == 0), logical(1))))
  expect_true(any(vapply(v, function(x) abs(x[1]) > 1e-8, logical(1))))
})

test_that("the generalized Pareto upper tail stays finite far out", {
  d <- gpd_distrib()
  th <- list(sigma = 1, xi = 0.1)
  for (q in c(50, 1e6, 1e12)) {
    v <- distrib_grad_cdf(d, q, th, lower.tail = FALSE, log = TRUE)
    expect_true(all(is.finite(unlist(v))))
  }
  # d log S / d sigma tends to 1/xi, but only once xi q / sigma is large: at
  # q = 50 it is 8.33 against 10, which is the function and not an error
  for (q in c(1e6, 1e12)) {
    v <- distrib_grad_cdf(d, q, th, lower.tail = FALSE, log = TRUE)
    expect_lt(abs(v[["sigma"]] - 1 / th$xi), 1e-3)
  }
  expect_true(all(is.finite(unlist(
    distrib_deriv4_cdf(d, 1e12, th, lower.tail = FALSE, log = TRUE)))))
})


# --- the families whose cdf is a sum of normal tails ------------------------

test_that("the inverse gaussian and the elastic net agree with a stencil on F", {
  hs <- c(1e-5, 1e-4, 1e-3, 3e-3)
  tols <- c(1e-6, 1e-5, 1e-3, 5e-3)
  ig <- invgauss1_distrib()
  en <- enet_distrib()
  cases <- list(
    list(ig, list(mu = 2, phi = 0.5), 1.3),
    list(ig, list(mu = 2, phi = 0.5), 0.3),
    list(ig, list(mu = 1, phi = 2), 3),
    list(ig, list(mu = 0.5, phi = 4), 0.2),
    # both sides of the location, the branch of the elastic net's cdf
    list(en, list(mu = 0.3, lambda = 2, alpha = 0.5), -0.4),
    list(en, list(mu = 0.3, lambda = 2, alpha = 0.5), 1.2),
    list(en, list(mu = 0, lambda = 1, alpha = 0.2), 0.8),
    list(en, list(mu = -0.5, lambda = 5, alpha = 0.7), -1.5),
    list(en, list(mu = 0, lambda = 0.5, alpha = 0.05), 2)
  )
  for (cs in cases) {
    d <- cs[[1]]; th <- cs[[2]]; q <- cs[[3]]
    for (o in 1:4) {
      got <- switch(o,
        distrib_grad_cdf(d, q, th, log = FALSE),
        distrib_hess_cdf(d, q, th, log = FALSE),
        distrib_deriv3_cdf(d, q, th, log = FALSE),
        distrib_deriv4_cdf(d, q, th, log = FALSE))
      ref <- vapply(deriv_indices(d@params, o),
                    function(I) fd_on_cdf(d, q, th, I, hs[o]), numeric(1))
      names(ref) <- deriv_names(d@params, o)
      expect_lt(max(abs(unlist(got)[names(ref)] - ref)) / max(1, max(abs(ref))),
                tols[o],
                label = sprintf("%s order %d at q = %g", d@distrib_name, o, q))
    }
  }
})

test_that("the inverse gaussian survives a weight that overflows alone", {
  # exp(2/(phi mu)) is Inf on its own at these settings while Phi(b)
  # underflows; the two are combined on the log scale
  d <- invgauss1_distrib()
  for (th in list(list(mu = 0.05, phi = 0.5), list(mu = 0.01, phi = 0.1))) {
    expect_true(is.infinite(exp(2 / (th$phi * th$mu))) ||
                  exp(2 / (th$phi * th$mu)) > 1e30)
    for (o in 1:4) {
      v <- switch(o,
        distrib_grad_cdf(d, 0.02, th, log = FALSE),
        distrib_hess_cdf(d, 0.02, th, log = FALSE),
        distrib_deriv3_cdf(d, 0.02, th, log = FALSE),
        distrib_deriv4_cdf(d, 0.02, th, log = FALSE))
      expect_true(all(is.finite(unlist(v))), info = sprintf("order %d", o))
    }
  }
})

test_that("the elastic net stays finite at both ends of its mixing", {
  # alpha near one is the Laplace and near zero the gaussian; the fourth
  # derivative in alpha grows as the family degenerates, which is the function
  d <- enet_distrib()
  for (al in c(0.001, 0.01, 0.5, 0.99, 0.999)) {
    v <- distrib_deriv4_cdf(d, 0.7, list(mu = 0, lambda = 2, alpha = al),
                            log = FALSE)
    expect_true(all(is.finite(unlist(v))), info = sprintf("alpha = %g", al))
  }
})

test_that("neither family is on the stencil any more", {
  for (d in list(invgauss1_distrib(), enet_distrib(), gpd_distrib())) {
    for (gen in list(distrib_grad_cdf, distrib_hess_cdf,
                     distrib_deriv3_cdf, distrib_deriv4_cdf)) {
      m <- S7::method(gen, S7::S7_class(d))
      owner <- attr(attr(m, "signature")[[1L]], "name")
      expect_false(owner %in% c("distrib", "continuous_distrib"),
                   info = sprintf("%s: %s", d@distrib_name, owner))
    }
  }
})


# --- the skew normal, in the shape as well ----------------------------------

test_that("the skew normal agrees with a stencil on F at every order", {
  hs <- c(1e-5, 1e-4, 1e-3, 3e-3)
  tols <- c(1e-7, 1e-6, 1e-4, 5e-3)
  d1 <- skewnormal1_distrib()
  d2 <- skewnormal2_distrib()
  for (cs in list(
        list(d1, list(mu = 0.2, sigma = 1.2, alpha = 1.5), 0.6),
        list(d1, list(mu = 0.2, sigma = 1.2, alpha = 1.5), -1.4),
        list(d1, list(mu = 0, sigma = 1, alpha = 0), 0.5),
        list(d1, list(mu = 0, sigma = 1, alpha = -4), 1.1),
        list(d2, list(mu = 0.3, sigma = 1.1, gamma1 = 0.4), 0.9),
        list(d2, list(mu = 0, sigma = 1, gamma1 = -0.6), -0.5))) {
    d <- cs[[1]]; th <- cs[[2]]; q <- cs[[3]]
    for (o in 1:4) {
      got <- switch(o,
        distrib_grad_cdf(d, q, th, log = FALSE),
        distrib_hess_cdf(d, q, th, log = FALSE),
        distrib_deriv3_cdf(d, q, th, log = FALSE),
        distrib_deriv4_cdf(d, q, th, log = FALSE))
      ref <- vapply(deriv_indices(d@params, o),
                    function(I) fd_on_cdf(d, q, th, I, hs[o]), numeric(1))
      names(ref) <- deriv_names(d@params, o)
      expect_lt(max(abs(unlist(got)[names(ref)] - ref)) / max(1, max(abs(ref))),
                tols[o],
                label = sprintf("%s order %d at q = %g", d@distrib_name, o, q))
    }
  }
})

test_that("at a zero shape the skew normal is the gaussian, exactly", {
  # the reference shares no code: the gaussian reaches these through the
  # location-scale construction and the skew normal through Owen's T
  d <- skewnormal1_distrib()
  g <- gaussian1_distrib()
  qq <- c(-1, 0.4, 2)
  th <- list(mu = 0.3, sigma = 1.2, alpha = 0)
  tg <- list(mu = 0.3, sigma = 1.2)
  for (o in 1:4) {
    a <- switch(o,
      distrib_grad_cdf(d, qq, th, log = FALSE),
      distrib_hess_cdf(d, qq, th, log = FALSE),
      distrib_deriv3_cdf(d, qq, th, log = FALSE),
      distrib_deriv4_cdf(d, qq, th, log = FALSE))
    b <- switch(o,
      distrib_grad_cdf(g, qq, tg, log = FALSE),
      distrib_hess_cdf(g, qq, tg, log = FALSE),
      distrib_deriv3_cdf(g, qq, tg, log = FALSE),
      distrib_deriv4_cdf(g, qq, tg, log = FALSE))
    for (k in names(b)) {
      expect_equal(a[[k]], b[[k]], tolerance = 1e-12,
                   info = sprintf("order %d component %s", o, k))
    }
  }
})

test_that("the first derivative in the quantile is the density", {
  # dF/dz = 2 phi(z) Phi(alpha z), which is sigma times the density: the check
  # that the Owen's T partial in its first argument is the right one
  d <- skewnormal1_distrib()
  th <- list(mu = 0.2, sigma = 1.3, alpha = 2)
  q <- c(-0.5, 0.4, 1.8)
  z <- (q - th$mu) / th$sigma
  expect_equal(2 * stats::dnorm(z) * stats::pnorm(th$alpha * z) / th$sigma,
               distrib_pdf(d, q, th), tolerance = 1e-14)
})

test_that("only the mathematical obstructions are left on the cdf stencil", {
  # every one of these differences its cdf because the derivative of an
  # incomplete gamma or beta in its shape is hypergeometric, or because the
  # distribution function is itself a quadrature. The test fails if a family
  # joins them, which is what would happen to a new one added without a route.
  exported <- grep("_distrib$", getNamespaceExports("distributions7"), value = TRUE)
  skip_these <- c("check_distrib", "fit_distrib", "continuous_distrib",
                  "discrete_distrib", "multivariate_distrib")
  open <- character()
  for (ctor in setdiff(exported, skip_these)) {
    d <- tryCatch(get(ctor, asNamespace("distributions7"))(),
                  error = function(e) NULL)
    if (is.null(d) || !S7::S7_inherits(d, distributions7:::continuous_distrib)) next
    if (!distributions7:::has_exact_cdf_deriv(d, 4L)) {
      open <- c(open, sub("_distrib$", "", ctor))
    }
  }
  expect_setequal(open, c("beta1", "beta2", "chisq", "gamma1", "gamma2",
                          "gengamma1", "vonmises1", "vonmises2"))
})
