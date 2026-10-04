# A single point evaluated at parameters that vary by observation.
#
# A compiled kernel sizes its output by the length of the response, so
# distrib_pdf(d, 0, theta) with theta of length n returned ONE value, read at
# the first observation's parameters, on eight families. zero_inflated()
# reads f(0) exactly that way, and with parameters varying by observation its
# derivatives were out by a relative 0.5 to 38 on those parents while its
# density at the full response was right. The generics now recycle the point
# as the derivative generics recycle the response.

scalar_loop <- function(f, d, x, ths) {
  vapply(seq_along(ths), function(i) f(d, x, ths[[i]]), numeric(1))
}
by_obs <- function(d, ths) {
  th <- lapply(d@params, function(p) vapply(ths, function(t) t[[p]], numeric(1)))
  stats::setNames(th, d@params)
}

test_that("pdf, cdf and quantile recycle a single point against vector parameters", {
  fams <- list(negbin1_distrib(), pig1_distrib(), pig2_distrib(),
               betabinom1_distrib(size = 10), betabinom2_distrib(size = 10),
               gengamma1_distrib(), gengamma2_distrib(), gpd_distrib(),
               poisson_distrib(), gaussian1_distrib())
  for (d in fams) {
    set.seed(3)
    ths <- replicate(4, generate_random_theta(d), simplify = FALSE)
    th <- by_obs(d, ths)
    x <- distrib_rng(d, 1, ths[[1]])
    pv <- distrib_pdf(d, x, th)
    expect_length(pv, 4)
    expect_equal(pv, scalar_loop(function(d, x, t) distrib_pdf(d, x, t), d, x, ths),
                 tolerance = 1e-12, label = d@distrib_name)
    expect_equal(distrib_pdf(d, x, th, log = TRUE),
                 scalar_loop(function(d, x, t) distrib_pdf(d, x, t, log = TRUE), d, x, ths),
                 tolerance = 1e-12, label = d@distrib_name)
    expect_equal(distrib_cdf(d, x, th),
                 scalar_loop(function(d, x, t) distrib_cdf(d, x, t), d, x, ths),
                 tolerance = 1e-10, label = d@distrib_name)
    expect_equal(distrib_quantile(d, 0.4, th),
                 scalar_loop(function(d, x, t) distrib_quantile(d, 0.4, t), d, x, ths),
                 tolerance = 1e-8, label = d@distrib_name)
  }
  # nothing moves where the lengths already agree, or where theta is scalar
  d <- pig1_distrib()
  expect_identical(recycle_point(d, c(0, 1), list(mu = c(1, 2), sigma = 1)), c(0, 1))
  expect_identical(recycle_point(d, 0, list(mu = 1, sigma = 1)), 0)
  expect_identical(recycle_point(d, 0, list(mu = c(1, 2, 3), sigma = 1)), c(0, 0, 0))
})


test_that("zero_inflated() derivatives are right with parameters by observation", {
  skip_on_cran()
  skip_if_not_installed("numDeriv")
  for (par in list(negbin1_distrib(), pig1_distrib(), pig2_distrib(),
                   betabinom2_distrib(size = 10), poisson_distrib())) {
    d <- zero_inflated(par)
    P <- d@params
    set.seed(7)
    ths <- replicate(5, generate_random_theta(d), simplify = FALSE)
    th <- by_obs(d, ths)
    y <- c(0, 0, 1, 3, 0)
    g <- distrib_gradient(d, y, th)
    h <- distrib_hessian(d, y, th)
    d3 <- distrib_deriv3(d, y, th)
    eh <- distrib_expected_hessian(d, y, th, approx = "bartlett")
    idx3 <- deriv_indices(P, 3L)
    for (i in seq_along(y)) {
      f <- function(u) stats::setNames(as.list(u), P)
      u0 <- unlist(ths[[i]][P])
      # order k against ONE Richardson pass on the analytic order k - 1
      ng <- numDeriv::grad(function(u) distrib_pdf(d, y[i], f(u), log = TRUE), u0)
      expect_equal(vapply(P, function(p) g[[p]][i], 0), ng, tolerance = 1e-7,
                   ignore_attr = TRUE, label = d@distrib_name)
      J <- numDeriv::jacobian(function(u) unlist(distrib_gradient(d, y[i], f(u))), u0)
      for (nm in names(h)) {
        ab <- match(hess_pairs(P)[[nm]], P)
        expect_equal(h[[nm]][i], J[ab[1], ab[2]], tolerance = 1e-7, label = nm)
      }
      Jh <- numDeriv::jacobian(function(u) unlist(distrib_hessian(d, y[i], f(u))), u0)
      for (r in seq_along(idx3)) {
        ix <- idx3[[r]]
        kh <- which(vapply(hess_pairs(P), function(ab) identical(sort(match(ab, P)), sort(ix[1:2])), TRUE))
        expect_equal(d3[[r]][i], Jh[kh, ix[3]], tolerance = 1e-6)
      }
      # the expected information read with parameters by observation against
      # the same quantity read at this observation's parameters alone: the
      # recycling this block is about. For the Poisson-inverse Gaussians that
      # is the whole reference -- their exactness is test-pig-expected.R's,
      # and a sum over 0..20000 with a row costing O(y) was 340 s of this
      # suite's 646 -- and the other families are summed over the support too.
      one <- distrib_expected_hessian(d, y[i], ths[[i]])
      for (nm in names(eh)) {
        expect_equal(eh[[nm]][i], one[[nm]], tolerance = 1e-12, label = nm)
      }
      if (grepl("poisson-inverse", par@distrib_name)) next
      # The pig2 draw at alpha = 0.67 has a tail heavy enough that a sum
      # stopped at 2000 is 2e-5 out on mu_mu while its mass already prints as
      # one, which is why the sum runs to 20000 where the support is not
      # bounded
      up <- if (is.finite(d@bounds[2])) d@bounds[2] else 20000
      yy <- 0:up
      pp <- distrib_pdf(d, yy, ths[[i]])
      k <- pp > 1e-300
      H <- distrib_hessian(d, yy[k], ths[[i]])
      # against the scale of the matrix: a component near zero -- mu_theta at
      # -1.4e-4 on one draw -- is 3e-13 out in absolute terms, 2.5e-9 relative
      sc <- max(abs(vapply(eh, function(v) v[i], 0)))
      for (nm in names(eh)) {
        expect_lt(abs(eh[[nm]][i] - sum(pp[k] * H[[nm]])), 1e-9 * sc, label = nm)
      }
    }
  }
})


test_that("the zero-inflated log density stays finite far in the tail", {
  d <- zero_inflated(poisson_distrib())
  th <- list(mu = 1, zi = 0.3)
  y <- c(0, 5, 200, 390)
  ref <- c(log(0.3 + 0.7 * dpois(0, 1)), log1p(-0.3) + dpois(y[-1], 1, log = TRUE))
  expect_equal(distrib_pdf(d, y, th, log = TRUE), ref, tolerance = 1e-12)
  expect_true(all(is.finite(distrib_pdf(d, y, th, log = TRUE))))
  # and the linear scale is untouched
  expect_equal(distrib_pdf(d, y, th),
               (1 - 0.3) * dpois(y, 1) + 0.3 * (y == 0), tolerance = 1e-15)
})


test_that("enet's third and fourth derivatives follow parameters by observation", {
  d <- enet_distrib()
  set.seed(3)
  ths <- replicate(4, generate_random_theta(d), simplify = FALSE)
  th <- by_obs(d, ths)
  y <- vapply(ths, function(t) distrib_rng(d, 1, t), 0)
  for (ord in 3:4) {
    fn <- if (ord == 3) distrib_deriv3 else distrib_deriv4
    v <- fn(d, y, th)
    for (nm in names(v)) {
      expect_equal(v[[nm]],
                   vapply(1:4, function(i) fn(d, y[i], ths[[i]])[[nm]], 0),
                   tolerance = 1e-12, label = nm)
    }
  }
})


test_that("betabinom2 keeps a parameter per observation aligned off the support", {
  d <- betabinom2_distrib(size = 10)
  th <- list(alpha = c(1, 2, 3, 4), beta = c(2, 2, 5, 1))
  y <- c(1, -1, 3, 11)
  v <- distrib_pdf(d, y, th, log = TRUE)
  expect_equal(v[c(2, 4)], c(-Inf, -Inf))
  expect_equal(v[c(1, 3)],
               c(distrib_pdf(d, 1, list(alpha = 1, beta = 2), log = TRUE),
                 distrib_pdf(d, 3, list(alpha = 3, beta = 5), log = TRUE)),
               tolerance = 1e-12)
})


# The third and fourth response and cdf derivatives took the point as given,
# so with a single point and parameters by observation the response
# derivatives came back length 1, read at the first observation's parameters,
# and the cdf derivatives at the right length with wrong values (-20.06 where
# the Poisson's mu_mu_mu at q = 3, mu = 8 is 0.0069). The truncated discrete
# wrapper reads F at its lower end that way, so its own third to fifth
# derivatives were wrong with parameters by observation. The four generics now
# recycle the point as the orders below them do.

point_families <- function() {
  ctors <- sort(grep("_distrib$", getNamespaceExports("distributions7"), value = TRUE))
  out <- list()
  for (nm in ctors) {
    f <- get(nm)
    if ("n_dim" %in% names(formals(f))) next
    d <- tryCatch(if ("size" %in% names(formals(f))) f(size = 10) else f(),
                  error = function(e) NULL)
    if (!is.null(d) && !S7::S7_inherits(d, multivariate_distrib)) out[[nm]] <- d
  }
  out
}
# component by component, the vector call against the scalar calls
expect_by_obs <- function(vec, refs, label, tolerance = 1e-8) {
  for (nm in names(refs[[1]])) {
    expect_length(vec[[nm]], length(refs))
    expect_equal(vec[[nm]], vapply(refs, function(r) r[[nm]][1], numeric(1)),
                 tolerance = tolerance, label = paste(label, nm))
  }
}

test_that("orders three and four recycle a single point against vector parameters", {
  fams <- point_families()
  expect_gte(length(fams), 42L)
  gens <- list(
    deriv3_y = function(d, x, t) distrib_deriv3_y(d, x, t),
    deriv4_y = function(d, x, t) distrib_deriv4_y(d, x, t),
    deriv3_cdf = function(d, x, t) distrib_deriv3_cdf(d, x, t),
    deriv4_cdf = function(d, x, t) distrib_deriv4_cdf(d, x, t))
  for (nm in names(fams)) {
    d <- fams[[nm]]
    set.seed(3)
    ths <- replicate(3, generate_random_theta(d), simplify = FALSE)
    th <- by_obs(d, ths)
    x <- distrib_rng(d, 1, ths[[1]])
    for (g in names(gens)) {
      refs <- tryCatch(lapply(ths, function(t) gens[[g]](d, x, t)), error = function(e) NULL)
      if (is.null(refs)) next                  # not defined for this family
      vec <- gens[[g]](d, x, th)
      if (!is.list(vec)) {
        expect_equal(vec, unlist(refs), tolerance = 1e-8, label = paste(nm, g))
      } else {
        expect_by_obs(vec, refs, paste(nm, g))
      }
    }
  }
})

test_that("the truncated discrete wrapper's higher derivatives follow parameters by observation", {
  y <- c(1, 3, 7)
  for (d in list(truncated(poisson_distrib(), lower = 1),
                 truncated(negbin2_distrib(), lower = 1, upper = 30))) {
    set.seed(4)
    ths <- replicate(3, generate_random_theta(d), simplify = FALSE)
    th <- by_obs(d, ths)
    for (f in list(distrib_hessian, distrib_deriv3, distrib_deriv4, distrib_deriv5)) {
      refs <- lapply(1:3, function(i) f(d, y[i], ths[[i]]))
      expect_by_obs(f(d, y, th), refs, d@distrib_name)
    }
  }
})

test_that("enet's constant response components recycle against vector parameters", {
  d <- enet_distrib()
  y <- c(-1.2, 0.3, 2.5)
  ths <- list(list(mu = 0, lambda = 1, alpha = 0.3), list(mu = 0.5, lambda = 2, alpha = 0.6),
              list(mu = -0.4, lambda = 0.7, alpha = 0.8))
  th <- by_obs(d, ths)
  hv <- distrib_hess_y(d, y, th)
  expect_length(hv, 3L)
  expect_equal(hv, vapply(1:3, function(i) distrib_hess_y(d, y[i], ths[[i]]), 0))
  expect_by_obs(distrib_cross_y(d, y, th), lapply(1:3, function(i) distrib_cross_y(d, y[i], ths[[i]])),
                "enet cross_y")
  # the numerical orders built on those two now have one entry per observation
  for (f in list(distrib_cross2_y, distrib_grad_y_hess, distrib_hess_y_hess)) {
    expect_by_obs(f(d, y, th), lapply(1:3, function(i) f(d, y[i], ths[[i]])), "enet",
                  tolerance = 1e-5)
  }
})
