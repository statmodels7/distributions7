# The generalized gamma in Stacy's form, whose parametrization is chosen so
# that the families it nests are read off the parameters rather than derived.

test_that("it nests the gamma, the Weibull and the exponential exactly", {
  d <- gengamma1_distrib()
  y <- c(0.5, 2, 5, 11)

  expect_equal(distrib_pdf(d, y, list(a = 2, d = 3, p = 1)),
               stats::dgamma(y, shape = 3, scale = 2))
  expect_equal(distrib_pdf(d, y, list(a = 2, d = 1.5, p = 1.5)),
               stats::dweibull(y, shape = 1.5, scale = 2))
  expect_equal(distrib_pdf(d, y, list(a = 2, d = 1, p = 1)),
               stats::dexp(y, rate = 1 / 2))

  # and so do the distribution functions, which is a separate claim
  expect_equal(distrib_cdf(d, y, list(a = 2, d = 3, p = 1)),
               stats::pgamma(y, shape = 3, scale = 2))
  expect_equal(distrib_cdf(d, y, list(a = 2, d = 1.5, p = 1.5)),
               stats::pweibull(y, shape = 1.5, scale = 2))
})


test_that("the density integrates to one and the quantile inverts the cdf", {
  d <- gengamma1_distrib()
  for (th in list(list(a = 2, d = 3, p = 1.5), list(a = 0.7, d = 0.8, p = 3))) {
    expect_equal(stats::integrate(function(z) distrib_pdf(d, z, th),
                                  0, Inf)$value, 1, tolerance = 1e-8)
    # The round trip is asked where the distribution function has not
    # saturated: at a = 0.7, d = 0.8, p = 3 it is exactly 1 in double
    # precision already at y = 4, and no quantile can come back from that.
    y <- distrib_quantile(d, c(0.05, 0.5, 0.95), th)
    expect_equal(distrib_quantile(d, distrib_cdf(d, y, th), th), y,
                 tolerance = 1e-8)
  }
})


test_that("the moments are the ratio of gamma functions", {
  # E[Y^r] = a^r Gamma((d+r)/p) / Gamma(d/p), written out here rather than
  # taken from the package.
  d <- gengamma1_distrib()
  th <- list(a = 2, d = 3, p = 1.5)
  for (r in 1:3) {
    want <- th$a^r * gamma((th$d + r) / th$p) / gamma(th$d / th$p)
    expect_equal(moment(d, th, p = r), want, tolerance = 1e-5,
                 label = as.character(r))
  }
  expect_equal(mean(d, th),
               th$a * gamma((th$d + 1) / th$p) / gamma(th$d / th$p),
               tolerance = 1e-6)
})


test_that("the analytical derivatives match one Richardson differentiation", {
  skip_if_not_installed("numDeriv")
  d <- gengamma1_distrib()
  th <- list(a = 2, d = 3, p = 1.5)
  q0 <- c(2, 3, 1.5)
  at <- function(q) list(a = q[1], d = q[2], p = q[3])
  set.seed(1)
  y <- distrib_rng(d, 60, th)

  g <- distrib_gradient(d, y, th)
  expect_equal(vapply(g, sum, numeric(1)),
    numDeriv::grad(function(q) sum(distrib_pdf(d, y, at(q), log = TRUE)), q0),
    tolerance = 1e-6, ignore_attr = TRUE)

  J <- numDeriv::jacobian(function(q) {
    vapply(distrib_gradient(d, y, at(q)), sum, numeric(1))
  }, q0)
  h <- distrib_hessian(d, y, th)
  pairs <- list(a_a = c(1, 1), a_d = c(1, 2), a_p = c(1, 3),
                d_d = c(2, 2), d_p = c(2, 3), p_p = c(3, 3))
  for (nm in names(pairs)) {
    expect_equal(sum(h[[nm]]), J[pairs[[nm]][1], pairs[[nm]][2]],
                 tolerance = 1e-6, label = nm)
  }
})


test_that("the expected information is closed form", {
  # Every expectation is a moment of u = (Y/a)^p, which is Gamma(d/p, 1). The
  # three components free of the data must come out EXACTLY equal to a Monte
  # Carlo mean of the observed Hessian, since there is nothing to average.
  d <- gengamma1_distrib()
  th <- list(a = 2, d = 3, p = 1.5)
  eh <- distrib_expected_hessian(d, 0, th)
  set.seed(2)
  y <- distrib_rng(d, 5e5, th)
  h <- distrib_hessian(d, y, th)

  for (nm in c("a_d", "d_d", "d_p")) {
    expect_equal(eh[[nm]][1], mean(h[[nm]]), tolerance = 1e-12, label = nm)
  }
  for (nm in c("a_a", "a_p", "p_p")) {
    expect_equal(eh[[nm]][1], mean(h[[nm]]), tolerance = 5e-3, label = nm)
  }
})


test_that("the validator passes and a fit with one shape held recovers the rest", {
  d <- gengamma1_distrib()
  set.seed(3)
  res <- check_distrib(d, verbose = FALSE)
  expect_true(all(res$status == "OK"),
    label = paste(res$check[res$status != "OK"], collapse = ", "))

  # The three parameters are weakly identified together, d and p entering
  # largely through their ratio, so the fit that is checked here holds one.
  set.seed(5)
  th <- list(a = 2, d = 3, p = 1.5)
  y <- distrib_rng(d, 3000, th)
  f <- fit_distrib(fixed(gengamma1_distrib(), p = 1.5), y)
  expect_true(f@converged)
  expect_equal(unname(coef(f)), c(2, 3), tolerance = 0.2)
})


test_that("the generalized gamma's assembly reproduces the compiled kernel", {
  # The order 3 and 4 components are assembled from the five terms of the
  # log-density; running that assembly at orders 1 and 2 must return what the
  # independently written C++ kernel returns, which is what licenses trusting
  # it at the orders where there is no kernel to compare against.
  d <- gengamma1_distrib()
  y <- c(0.4, 1.1, 2.3)
  for (th in list(list(a = 1.3, d = 2.1, p = 1.6),
                  list(a = 0.7, d = 0.5, p = 3.2))) {
    g <- distributions7:::gengamma_components(y, th, 1L)
    r <- distrib_gradient(d, y, th)
    for (nm in names(r)) expect_equal(g[[nm]], r[[nm]], tolerance = 1e-12)
    h <- distributions7:::gengamma_components(y, th, 2L)
    r <- distrib_hessian(d, y, th)
    for (nm in names(r)) expect_equal(h[[nm]], r[[nm]], tolerance = 1e-12)
  }
})


test_that("the generalized gamma is closed at third and fourth order", {
  skip_if_not_installed("numDeriv")
  d <- gengamma1_distrib()
  y <- c(0.4, 1.1, 2.3)
  th <- list(a = 1.3, d = 2.1, p = 1.6)
  for (ord in 3:4) {
    lower <- if (ord == 3L) distrib_hessian else distrib_deriv3
    got <- if (ord == 3L) distrib_deriv3(d, y, th) else distrib_deriv4(d, y, th)
    for (nm in names(got)) {
      parts <- strsplit(nm, "_")[[1]]
      j <- match(parts[length(parts)], d@params)
      head_nm <- paste(parts[-length(parts)], collapse = "_")
      ref <- vapply(seq_along(y), function(i) {
        numDeriv::grad(function(z) {
          t2 <- th; t2[[j]] <- z
          lower(d, y[i], t2)[[head_nm]]
        }, th[[j]])
      }, numeric(1))
      expect_equal(got[[nm]], ref, tolerance = 1e-6, label = nm)
    }
  }
})
