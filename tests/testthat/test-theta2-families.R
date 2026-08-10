# The location-scale identity at second order.

# ONE Richardson pass on the ANALYTIC first-order component: cross_y is
# d2l/dy dtheta and cross2_y is d3l/dy2 dtheta, both written out, so
# differentiating either once more in another parameter is a single pass on an
# analytic quantity. Nesting numDeriv to reach the second order instead makes
# the reference the weak side -- measured, it disagreed with correct code by
# up to 7e-4 relative, which is Richardson's own accuracy at that depth and
# not an error in what it was checking.
ref_theta2 <- function(d, y, th, order, pr) {
  inner <- if (order == 1L) distrib_cross_y else distrib_cross2_y
  vapply(seq_along(y), function(i) {
    numDeriv::grad(function(v) {
      t2 <- th
      t2[[pr[1L]]] <- v
      inner(d, y[i], t2)[[pr[2L]]]
    }, th[[pr[1L]]])
  }, numeric(1))
}

# a comparison that does not let a near-zero entry dominate: these components
# range over several orders of magnitude within one vector
expect_close <- function(got, ref, tol = 1e-8, info = NULL) {
  expect_lt(max(abs(got - ref)) / max(1e-3, max(abs(ref))), tol)
}

test_that("the identity closes the location-scale families", {
  skip_if_not_installed("numDeriv")
  y <- c(-2.4, -0.6, 0.7, 1.9)
  cases <- list(
    list(d = logistic_distrib(), th = list(mu = 0.3, sigma = 1.2)),
    list(d = cauchy_distrib(), th = list(mu = -0.4, sigma = 0.9)),
    list(d = gumbel_distrib(), th = list(mu = 0.1, sigma = 1.5))
  )
  for (cs in cases) {
    d <- cs$d
    th <- cs$th
    p <- d@params
    prs <- list(c(p[1], p[1]), c(p[2], p[2]), c(p[1], p[2]))
    for (order in 1:2) {
      got <- if (order == 1L) distrib_grad_y_hess(d, y, th) else
        distrib_hess_y_hess(d, y, th)
      for (pr in prs) {
        nm <- paste(pr, collapse = "_")
        expect_close(got[[nm]], ref_theta2(d, y, th, order, pr))
      }
    }
  }
})

test_that("the gaussian's written-out method agrees with the identity", {
  # the gaussian has its own method, written out; the identity reaches the
  # same numbers from the family's response derivatives, so the two are an
  # independent check of each other
  y <- c(-1.1, 0.4, 2.2)
  d <- gaussian1_distrib()
  th <- list(mu = 0.3, sigma = 1.4)
  expect_equal(distrib_grad_y_hess(d, y, th),
               distributions7:::loc_scale_grad_y_hess(d, y, th),
               tolerance = 1e-12)
  expect_equal(distrib_hess_y_hess(d, y, th),
               distributions7:::loc_scale_hess_y_hess(d, y, th),
               tolerance = 1e-12)
})

test_that("a shape family keeps its location-scale pairs closed", {
  skip_if_not_installed("numDeriv")
  y <- c(-1.8, 0.2, 1.6)
  cases <- list(
    list(d = student_t1_distrib(), th = list(mu = 0.2, sigma = 1.1, nu = 7)),
    list(d = pseudohuber_distrib(), th = list(mu = 0.1, sigma = 1.3, nu = 2))
  )
  for (cs in cases) {
    d <- cs$d
    th <- cs$th
    p <- d@params
    for (order in 1:2) {
      got <- if (order == 1L) distrib_grad_y_hess(d, y, th) else
        distrib_hess_y_hess(d, y, th)
      expect_named(got, hess_names(p))
      # the three location-scale pairs, closed
      for (pr in list(c(p[1], p[1]), c(p[2], p[2]), c(p[1], p[2]))) {
        nm <- paste(pr, collapse = "_")
        expect_close(got[[nm]], ref_theta2(d, y, th, order, pr))
      }
      # and the pairs touching the shape, which fall back and still answer
      for (nm in c(paste(p[1], p[3], sep = "_"),
                   paste(p[2], p[3], sep = "_"),
                   paste(p[3], p[3], sep = "_"))) {
        expect_true(all(is.finite(got[[nm]])), info = nm)
      }
    }
  }
})

test_that("the closed block is not merely the fallback under another name", {
  # the identity reads distrib_deriv3_y and distrib_deriv4_y and scales them;
  # the fallback differences distrib_cross_y. They agree to the fallback's own
  # accuracy, which is what says the identity is right, and they are not the
  # same arithmetic.
  skip_if_not_installed("numDeriv")
  y <- c(-1.5, 0.3, 2.0)
  d <- logistic_distrib()
  th <- list(mu = 0.2, sigma = 1.3)
  closed <- distrib_grad_y_hess(d, y, th)
  differenced <- numerical_theta2_y(d, y, th,
                                    function(t) distrib_cross_y(d, y, t))
  for (nm in names(closed)) expect_close(closed[[nm]], differenced[[nm]])
  ch <- distrib_hess_y_hess(d, y, th)
  dh <- numerical_theta2_y(d, y, th, function(t) distrib_cross2_y(d, y, t))
  for (nm in names(ch)) expect_close(ch[[nm]], dh[[nm]])
})

test_that("closing cross2_y is what makes the reference clean", {
  # Before the identity was registered on distrib_cross2_y, the numerical
  # routes at order 2 were differencing a DIFFERENCED quantity, and the three
  # routes disagreed among themselves by up to 3e-5 -- including the two
  # numerical ones with each other, which is what said the code was not the
  # weak side. With cross2_y closed they agree at 1e-10.
  skip_if_not_installed("numDeriv")
  for (d in list(logistic_distrib(), cauchy_distrib(), gumbel_distrib())) {
    p <- d@params
    th <- stats::setNames(list(0.2, 1.3), p)
    y <- c(-1.5, 0.3, 2.0)
    got <- distrib_cross2_y(d, y, th)
    for (q in p) {
      ref <- vapply(seq_along(y), function(i) {
        numDeriv::grad(function(v) {
          t2 <- th
          t2[[q]] <- v
          distrib_hess_y(d, y[i], t2)
        }, th[[q]])
      }, numeric(1))
      expect_close(got[[q]], ref, info = paste(d@distrib_name, q))
    }
  }
})
