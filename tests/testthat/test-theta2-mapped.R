# The mixed grid carried through a reparametrization.

test_that("a reparametrized family's mixed grid matches Richardson", {
  skip_if_not_installed("numDeriv")
  yr <- c(-1.6, 0.3, 1.9)
  yp <- c(0.4, 1.2, 2.7)
  cases <- list(
    list(d = gaussian2_distrib(), th = list(mu = 0.4, sigma2 = 1.7), y = yr),
    list(d = gaussian3_distrib(), th = list(mu = 0.4, tau = 0.6), y = yr),
    list(d = student_t2_distrib(), th = list(mu = 0.2, sigma = 1.3, nu = 8),
         y = yr),
    list(d = lognormal2_distrib(), th = list(mean = 1.4, var = 0.7), y = yp)
  )
  for (cs in cases) {
    d <- cs$d
    th <- cs$th
    y <- cs$y
    p <- d@params
    # cross2_y: one theta-derivative of the analytic response curvature
    got1 <- distrib_cross2_y(d, y, th)
    for (q in p) {
      ref <- vapply(seq_along(y), function(i) {
        numDeriv::grad(function(v) {
          t2 <- th
          t2[[q]] <- v
          distrib_hess_y(d, y[i], t2)
        }, th[[q]])
      }, numeric(1))
      expect_lt(max(abs(got1[[q]] - ref)) / max(1e-3, max(abs(ref))), 1e-7)
    }
    # and the two second-order ones, against ONE Richardson pass on the
    # analytic first-order component
    for (order in 1:2) {
      got <- if (order == 1L) distrib_grad_y_hess(d, y, th) else
        distrib_hess_y_hess(d, y, th)
      inner <- if (order == 1L) distrib_cross_y else distrib_cross2_y
      prs <- distributions7:::hess_pairs(p)
      for (nm in names(prs)) {
        pr <- prs[[nm]]
        ref <- vapply(seq_along(y), function(i) {
          numDeriv::grad(function(v) {
            t2 <- th
            t2[[pr[1L]]] <- v
            inner(d, y[i], t2)[[pr[2L]]]
          }, th[[pr[1L]]])
        }, numeric(1))
        # a mapped family is as exact as its parent: student_t2's pairs
        # touching nu route through student_t1's differenced ones, so the
        # tolerance here is the parent's and not the chain rule's
        expect_lt(max(abs(got[[nm]] - ref)) / max(1e-3, max(abs(ref))), 1e-5)
      }
    }
  }
})

test_that("the map's second partials are what the chain rule needs", {
  # gaussian2 is the gaussian at sigma = sqrt(sigma2), so d sigma/d sigma2 =
  # 1/(2 sqrt(sigma2)) and d2 sigma/d sigma2^2 = -1/(4 sigma2^(3/2)). The
  # sigma2_sigma2 component therefore carries BOTH the parent's second
  # derivative and its first through that curvature, and a version that
  # dropped the second term would still pass a first-order check.
  y <- c(-1.1, 0.7)
  s2 <- 1.7
  s <- sqrt(s2)
  th2 <- list(mu = 0.4, sigma2 = s2)
  th1 <- list(mu = 0.4, sigma = s)
  d1 <- gaussian1_distrib()
  d2 <- gaussian2_distrib()

  ds <- 1 / (2 * s)
  d2s <- -1 / (4 * s2^1.5)
  got <- distrib_grad_y_hess(d2, y, th2)$sigma2_sigma2
  by_hand <- distrib_grad_y_hess(d1, y, th1)$sigma_sigma * ds^2 +
    distrib_cross_y(d1, y, th1)$sigma * d2s
  # the two orderings of the same arithmetic, so the gap is rounding
  expect_equal(got, by_hand, tolerance = 1e-9)
  # and the term that would be dropped is not negligible
  expect_gt(max(abs(distrib_cross_y(d1, y, th1)$sigma * d2s)), 1e-3)
})

test_that("the laplace is registered at both orders", {
  # it was registered for cross2_y and not for the two second-order generics,
  # which the census caught: a family half-registered answers correctly at one
  # order and falls back at the next without anything failing
  d <- laplace_distrib()
  expect_false(identical(
    attr(attr(S7::method(distrib_grad_y_hess, S7::S7_class(d)),
              "signature")[[1]], "name"),
    "continuous_distrib"))
  expect_false(identical(
    attr(attr(S7::method(distrib_hess_y_hess, S7::S7_class(d)),
              "signature")[[1]], "name"),
    "continuous_distrib"))
})
