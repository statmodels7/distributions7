# The link scale at the fifth order, by two routes that share no arithmetic
# beyond the parameter-scale fourth.
#
#   (A) differentiate in eta the order-4 component ALREADY on the link scale.
#       That is what distrib_deriv5() does, and it needs neither the fifth
#       derivative of a link nor an order-5 entry in bell_partial().
#   (B) differentiate in theta the order-4 component on the PARAMETER scale and
#       carry it over with Faa di Bruno at order 5, which is what
#       to_link_scale() can do now: it consumes d5linkinv() from
#       linkfunctions7 and the order-5 Bell row.
#
# Two implementations of one object is the idiom this toolkit prefers to a
# tolerance chosen by hand, and it is the reason the links were done first.

skip_if_not_installed("numDeriv")


rel5 <- function(a, b) {
  glob <- max(vapply(b, function(v) max(abs(v)), numeric(1)), 1e-300)
  max(vapply(names(b), function(k) {
    den <- max(abs(b[[k]]))
    max(abs(a[[k]] - b[[k]])) / if (den > glob * 1e-12) den else glob
  }, numeric(1)))
}

# route (B), assembled in the test rather than taken from the package, so the
# comparison is between two things and not between one thing and itself
route_B <- function(d, y, theta) {
  nat <- link_scale_lower_orders(d, y, theta, FALSE, 5L)
  res <- numerical_deriv5(d, y, theta, scale = "parameter")
  to_link_scale(d, theta, c(nat, list(res)), 5L)
}

ls5_cases <- function() {
  list(
    gaussian1  = list(d = gaussian1_distrib(),  th = list(mu = 0.3, sigma = 1.2),
                      y = c(-1.3, -0.2, 0.4, 1.1)),
    gamma1     = list(d = gamma1_distrib(),     th = list(mu = 2.0, phi = 0.7),
                      y = c(0.4, 0.9, 1.6, 2.7)),
    poisson    = list(d = poisson_distrib(),    th = list(mu = 2.5),
                      y = c(0, 1, 2, 4)),
    negbin2    = list(d = negbin2_distrib(),    th = list(mu = 3.0, theta = 2.2),
                      y = c(0, 1, 3, 5)),
    student_t1 = list(d = student_t1_distrib(), th = list(mu = 0.2, sigma = 1.1, nu = 6),
                      y = c(-1.3, -0.2, 0.4, 1.1)),
    beta1      = list(d = beta1_distrib(),      th = list(mu = 0.4, phi = 6),
                      y = c(0.15, 0.31, 0.5, 0.68))
  )
}


test_that("the two routes to the link-scale fifth agree", {
  for (nm in names(ls5_cases())) {
    cs <- ls5_cases()[[nm]]
    A <- distrib_deriv5(cs$d, cs$y, cs$th, scale = "link")
    B <- route_B(cs$d, cs$y, cs$th)
    expect_lt(rel5(A, B), 1e-6, label = sprintf("%s: route A against route B", nm))
  }
})


test_that("and both agree with Richardson on the analytic fourth", {
  # The third opinion. Where two of the three agree the odd one out is the
  # third, which is what separated a family from a bad probe when the fifth
  # order was first measured.
  for (nm in names(ls5_cases())) {
    cs <- ls5_cases()[[nm]]
    params <- cs$d@params
    idx <- deriv_indices(params, 5)
    nms <- deriv_names(params, 5)
    links <- cs$d@link_params
    x0 <- vapply(params, function(p)
      linkfunctions7::linkfun(links[[p]], cs$th[[p]]), numeric(1))

    ref <- vector("list", length(nms))
    names(ref) <- nms
    for (t in seq_along(nms)) {
      m <- idx[[t]][5]
      key4 <- paste(params[idx[[t]][1:4]], collapse = "_")
      f <- function(v) {
        x <- x0
        x[m] <- v
        th2 <- cs$th
        for (i in seq_along(params)) {
          th2[[params[i]]] <- linkfunctions7::linkinv(links[[params[i]]], x[i])
        }
        distrib_deriv4(cs$d, cs$y, th2, scale = "link")[[key4]]
      }
      ref[[t]] <- as.vector(numDeriv::jacobian(f, x0[m]))
    }

    expect_lt(rel5(distrib_deriv5(cs$d, cs$y, cs$th, scale = "link"), ref), 1e-6,
              label = sprintf("%s: route A against Richardson", nm))
    expect_lt(rel5(route_B(cs$d, cs$y, cs$th), ref), 1e-6,
              label = sprintf("%s: route B against Richardson", nm))
  }
})


test_that("an identity link leaves every order alone", {
  # The control that costs nothing and would catch a chain applied where it
  # should not be: on the identity link the two scales are the same scale.
  d <- gaussian1_distrib(link_sigma = linkfunctions7::identity_link())
  th <- list(mu = 0.3, sigma = 1.2)
  y <- c(-1.3, -0.2, 0.4, 1.1)

  P <- distrib_deriv5(d, y, th, scale = "parameter")
  expect_lt(rel5(distrib_deriv5(d, y, th, scale = "link"), P), 1e-12)
  expect_lt(rel5(route_B(d, y, th), P), 1e-12)
})


test_that("the order-5 Bell row is the one an independent construction gives", {
  # B_{m,j} is the coefficient of t^m/m! in (sum_k x_k t^k/k!)^j / j!. Building
  # it by convolution shares no arithmetic with the written-out table, and the
  # coefficients summing to the Bell number B_5 = 52 is a second check on the
  # same row.
  bell_by_series <- function(x, n = 5L) {
    a <- c(0, x / factorial(seq_len(n)))
    out <- numeric(n)
    for (j in seq_len(n)) {
      p <- c(1, rep(0, n))
      for (r in seq_len(j)) {
        q <- numeric(n + 1)
        for (u in 0:n) {
          for (v in 0:(n - u)) q[u + v + 1] <- q[u + v + 1] + p[u + 1] * a[v + 1]
        }
        p <- q
      }
      out[j] <- p[n + 1] * factorial(n) / factorial(j)
    }
    out
  }

  set.seed(5)
  x <- rnorm(5)
  h <- as.list(x)
  tab <- vapply(1:5, function(j) bell_partial(5L, j, h), numeric(1))
  expect_equal(tab, bell_by_series(x), tolerance = 1e-12)

  # the coefficients of the row, read off at x = 1, sum to B_5
  ones <- as.list(rep(1, 5))
  expect_identical(sum(vapply(1:5, function(j) bell_partial(5L, j, ones), numeric(1))), 52)

  # and the orders below are unchanged, which is what says the table was
  # extended rather than rewritten
  expect_equal(bell_partial(4L, 2L, h), 4 * x[1] * x[3] + 3 * x[2]^2)
  expect_equal(bell_partial(3L, 2L, h), 3 * x[1] * x[2])
})


test_that("a wrong Bell coefficient would fail the comparison", {
  # The negative control for the route-B assembly: without it the agreement
  # above could be satisfied by something other than these polynomials.
  d <- gaussian1_distrib()
  th <- list(mu = 0.3, sigma = 1.2)
  y <- c(-1.3, -0.2, 0.4, 1.1)

  nat <- link_scale_lower_orders(d, y, th, FALSE, 5L)
  res <- numerical_deriv5(d, y, th, scale = "parameter")
  good <- to_link_scale(d, th, c(nat, list(res)), 5L)

  # 5 -> 5.25 in B_{5,2}, everything else exact. The original has to be
  # captured BEFORE the mock is installed: `bell_partial` is looked up in the
  # namespace at call time, so a fallback naming it would call the mock and
  # recurse until the node stack goes.
  orig_bell <- bell_partial
  bad_bell <- function(m, j, h) {
    if (m == 5L && j == 2L) return(5.25 * h[[1]] * h[[4]] + 10 * h[[2]] * h[[3]])
    orig_bell(m, j, h)
  }
  local_mocked_bindings(bell_partial = bad_bell)
  bad <- to_link_scale(d, th, c(nat, list(res)), 5L)

  expect_gt(rel5(bad, good), 1e-3)
})


test_that("the sixth order is still refused, and the message says which", {
  h <- as.list(rep(1, 6))
  expect_error(bell_partial(6L, 1L, h), "up to order 5")
})
