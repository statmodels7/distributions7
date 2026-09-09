## The skew t's derivatives in (mu, sigma, alpha) are closed form at orders
## three and four.  Four things are asserted, and the first two are what make
## the other two mean anything.
##
##  (i)   the same assembly run at orders ONE and TWO reproduces the family's
##        independently hand-written score and Hessian.  Those were derived
##        separately and are already under check_distrib(), so agreement where
##        a closed form exists is what licenses the orders where none does;
##  (ii)  order three against Richardson applied to the family's own
##        hand-written Hessian, and order four against one second difference of
##        it -- references that do not run through the closed forms, so that an
##        error common to two orders cannot cancel;
##  (iii) the methods really return the closed forms, asserted by IDENTITY,
##        which a stencil cannot satisfy however accurate it is;
##  (iv)  the components carrying exactly one nu are closer to Richardson than
##        the route they replace, by a wide margin.

d <- skewt_distrib()
P <- d@params
y <- c(-1.7, -0.6, 0.2, 0.9, 2.3)
th <- list(mu = 0.2, sigma = 1.3, alpha = 0.7, nu = 8)

rel <- function(got, ref) max(abs(got - ref)) / max(max(abs(ref)), 1e-10)

## The Hessian names the diagonal first, so a pair's key comes from
## hess_names() rather than from pasting the two parameter names.
hess_key <- function(params, ij) {
  want <- sort(params[ij])
  pr <- hess_pairs(params)
  for (k in seq_along(pr)) {
    if (identical(sort(pr[[k]]), want)) {
      return(names(pr)[k])
    }
  }
  stop("no Hessian key for that index pair")
}

test_that("the closed block reproduces the hand-written score and Hessian", {
  tw <- skewt_msa_tower(y, th$mu, th$sigma, th$alpha, th$nu)
  cmp <- function(a, b, c) skewt_msa_component(tw, th$sigma, a, b, c)

  g <- distrib_gradient(d, y, th)
  expect_lt(rel(cmp(1L, 0L, 0L), g$mu), 1e-14)
  expect_lt(rel(cmp(0L, 1L, 0L), g$sigma), 1e-14)
  expect_lt(rel(cmp(0L, 0L, 1L), g$alpha), 1e-14)

  h <- distrib_hessian(d, y, th)
  expect_lt(rel(cmp(2L, 0L, 0L), h$mu_mu), 1e-14)
  expect_lt(rel(cmp(0L, 2L, 0L), h$sigma_sigma), 1e-14)
  expect_lt(rel(cmp(0L, 0L, 2L), h$alpha_alpha), 1e-14)
  expect_lt(rel(cmp(1L, 1L, 0L), h$mu_sigma), 1e-14)
  expect_lt(rel(cmp(1L, 0L, 1L), h$mu_alpha), 1e-14)
  expect_lt(rel(cmp(0L, 1L, 1L), h$sigma_alpha), 1e-14)
})

test_that("order three agrees with Richardson on the hand-written Hessian", {
  skip_if_not_installed("numDeriv")
  # The reference differentiates distrib_hessian(), which was derived
  # separately and shares no arithmetic with the closed forms.  Differentiating
  # the closed form instead would let an error common to two orders cancel:
  # measured, a uniform 5 per cent inflation passes such a check.
  x0 <- c(th$mu, th$sigma, th$alpha)
  hess_at <- function(x, nm) {
    t2 <- list(mu = x[1], sigma = x[2], alpha = x[3], nu = th$nu)
    distrib_hessian(d, y, t2)[[nm]]
  }
  compared <- 0L
  worst <- 0
  idx <- deriv_indices(P, 3L)
  for (t in seq_along(idx)) {
    ii <- idx[[t]]
    if (any(ii == 4L)) next
    m <- ii[3]
    nm2 <- hess_key(P, ii[1:2])
    ref <- as.vector(numDeriv::jacobian(
      function(v) {
        x <- x0
        x[m] <- v
        hess_at(x, nm2)
      }, x0[m]))
    tw <- skewt_msa_tower(y, x0[1], x0[2], x0[3], th$nu)
    got <- skewt_msa_component(tw, x0[2], sum(ii == 1L), sum(ii == 2L), sum(ii == 3L))
    worst <- max(worst, rel(got, ref))
    compared <- compared + 1L
  }
  # A skip inside the loop would abandon the block silently, so the count is
  # asserted rather than trusted.
  expect_identical(compared, 10L)
  expect_lt(worst, 1e-7)
})

test_that("order four agrees with a second difference of the Hessian", {
  # numerical_deriv4() is the route this release replaces: one second
  # difference of the analytic Hessian, whose own error was measured at 1e-5.
  # The tolerance sits well above that and well below any defect worth
  # catching, so this block fails on its own rather than leaning on the one
  # above.
  was <- numerical_deriv4(d, y, th)
  m4 <- skewt_msa_derivs(d, y, th, 4L)
  expect_identical(length(m4), 15L)
  worst <- 0
  for (nm in names(m4)) worst <- max(worst, rel(m4[[nm]], was[[nm]]))
  expect_lt(worst, 1e-3)
})

test_that("the Richardson check catches a corrupted component", {
  skip_if_not_installed("numDeriv")
  x0 <- c(th$mu, th$sigma, th$alpha)
  tw <- skewt_msa_tower(y, x0[1], x0[2], x0[3], th$nu)
  # (mu, mu, sigma): differentiate the analytic (mu, mu) in sigma
  ref <- as.vector(numDeriv::jacobian(
    function(v) {
      t2 <- skewt_msa_tower(y, x0[1], v, x0[3], th$nu)
      skewt_msa_component(t2, v, 2L, 0L, 0L)
    }, x0[2]))
  good <- skewt_msa_component(tw, x0[2], 2L, 1L, 0L)
  expect_lt(rel(good, ref), 1e-7)
  expect_gt(rel(1.05 * good, ref), 1e-3)
})

test_that("deriv3 and deriv4 return the closed forms themselves", {
  d3 <- distrib_deriv3(d, y, th)
  d4 <- distrib_deriv4(d, y, th)
  m3 <- skewt_msa_derivs(d, y, th, 3L)
  m4 <- skewt_msa_derivs(d, y, th, 4L)

  expect_identical(length(m3), 10L)
  expect_identical(length(m4), 15L)
  # Identity, not a tolerance: an accurate stencil would pass a tolerance.
  for (nm in names(m3)) expect_identical(d3[[nm]], m3[[nm]])
  for (nm in names(m4)) expect_identical(d4[[nm]], m4[[nm]])

  # No component is lost to the skipping, at either order.
  expect_identical(names(d3), deriv_names(P, 3L))
  expect_identical(names(d4), deriv_names(P, 4L))
  expect_false(any(vapply(d3, is.null, logical(1))))
  expect_false(any(vapply(d4, is.null, logical(1))))
})

test_that("the components carrying one nu beat the route they replace", {
  skip_if_not_installed("numDeriv")
  now <- distrib_deriv4(d, y, th)
  was <- numerical_deriv4(d, y, th)
  i4 <- deriv_indices(P, 4L)
  n4 <- deriv_names(P, 4L)
  compared <- 0L
  worst_now <- 0
  worst_was <- 0
  for (t in seq_along(n4)) {
    ii <- i4[[t]]
    if (sum(ii == 4L) != 1L) next
    j <- ii[ii != 4L]
    a <- sum(j == 1L)
    b <- sum(j == 2L)
    c <- sum(j == 3L)
    ref <- as.vector(numDeriv::jacobian(
      function(v) {
        tw <- skewt_msa_tower(y, th$mu, th$sigma, th$alpha, v)
        skewt_msa_component(tw, th$sigma, a, b, c)
      }, th$nu))
    worst_now <- max(worst_now, rel(now[[n4[t]]], ref))
    worst_was <- max(worst_was, rel(was[[n4[t]]], ref))
    compared <- compared + 1L
  }
  expect_identical(compared, 10L)
  expect_lt(worst_now, 1e-7)
  # The margin measured over fifteen settings of (nu, alpha) is 20x to 203x;
  # 5x is asked so the assertion is about the route and not about one cell.
  expect_gt(worst_was / worst_now, 5)
})

test_that("skip = NULL leaves the generic construction untouched", {
  expect_identical(numerical_deriv3(d, y, th), numerical_deriv3(d, y, th, skip = NULL))
  expect_identical(numerical_deriv4(d, y, th), numerical_deriv4(d, y, th, skip = NULL))
  # And a skipped component is left in place rather than dropped from the list.
  s <- numerical_deriv3(d, y, th, skip = "mu_mu_mu")
  expect_identical(names(s), deriv_names(P, 3L))
  expect_null(s[["mu_mu_mu"]])
  expect_false(is.null(s[["mu_mu_sigma"]]))
})

test_that("the closed block holds across the parameter space", {
  skip_if_not_installed("numDeriv")
  # Third order only, and against the hand-written Hessian, for the reason the
  # block above gives.  Fourth order is covered there over the same grid by
  # check_distrib()'s own battery.
  compared <- 0L
  worst <- 0
  for (nu in c(3, 8, 50)) {
    for (alpha in c(-2, 3)) {
      x0 <- c(0.2, 1.3, alpha)
      for (ii in list(c(1L, 1L, 1L), c(2L, 2L, 3L), c(1L, 3L, 3L),
                      c(1L, 2L, 2L), c(3L, 3L, 3L))) {
        m <- ii[3]
        nm2 <- hess_key(P, ii[1:2])
        ref <- as.vector(numDeriv::jacobian(
          function(v) {
            x <- x0
            x[m] <- v
            t2 <- list(mu = x[1], sigma = x[2], alpha = x[3], nu = nu)
            distrib_hessian(d, y, t2)[[nm2]]
          }, x0[m]))
        tw <- skewt_msa_tower(y, x0[1], x0[2], x0[3], nu)
        got <- skewt_msa_component(tw, x0[2], sum(ii == 1L), sum(ii == 2L),
                                   sum(ii == 3L))
        worst <- max(worst, rel(got, ref))
        compared <- compared + 1L
      }
    }
  }
  expect_identical(compared, 30L)
  expect_lt(worst, 1e-6)
})
