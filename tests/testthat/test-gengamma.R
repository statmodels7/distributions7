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


test_that("the generalized gamma's assembly reproduces the compiled kernels", {
  # gengamma_components() assembles orders one to four from the five terms of
  # the log-density; the kernels are generated from the closed form in U and
  # e^U. The two routes were written independently. Below k = d/p of about
  # 0.01 the assembly is the weak side: it reads psi^(n)(k), whose poles the
  # kernels cancel in closed form (the next tests but one).
  d <- gengamma1_distrib()
  y <- c(0.4, 1.1, 2.3)
  gens <- list(distrib_gradient, distrib_hessian, distrib_deriv3, distrib_deriv4)
  for (th in list(list(a = 1.3, d = 2.1, p = 1.6),
                  list(a = 0.7, d = 0.5, p = 3.2),
                  list(a = 2, d = 0.4, p = 2))) {
    for (ord in 1:4) {
      g <- distributions7:::gengamma_components(y, th, ord)
      r <- gens[[ord]](d, y, th)
      expect_setequal(names(r), names(g))
      for (nm in names(r)) {
        expect_equal(r[[nm]], g[[nm]], tolerance = 1e-12, label = nm)
      }
    }
  }
})


test_that("the generalized gamma's fifth order is the derivative of its fourth", {
  skip_if_not_installed("numDeriv")
  d <- gengamma1_distrib()
  y <- c(0.4, 1.1, 2.3)
  th <- list(a = 1.3, d = 2.1, p = 1.6)
  got <- distrib_deriv5(d, y, th)
  expect_length(got, 21L)
  for (nm in names(got)) {
    parts <- strsplit(nm, "_")[[1]]
    j <- match(parts[5], d@params)
    head_nm <- paste(parts[1:4], collapse = "_")
    ref <- vapply(seq_along(y), function(i) {
      numDeriv::grad(function(z) {
        t2 <- th; t2[[j]] <- z
        distrib_deriv4(d, y[i], t2)[[head_nm]]
      }, th[[j]])
    }, numeric(1))
    expect_equal(got[[nm]], ref, tolerance = 1e-7, label = nm)
  }
})


test_that("the generalized gamma's expected higher derivatives are the integrals", {
  # the closed forms against quadrature of the observed derivatives over the
  # density, on the scale of u = (y/a)^p where the integrand is smooth
  d <- gengamma1_distrib()
  th <- list(a = 1.3, d = 2.1, p = 1.6)
  k <- th$d / th$p
  for (ord in 3:4) {
    got <- if (ord == 3L) distrib_deriv3(d, 1, th, expected = TRUE)
           else distrib_deriv4(d, 1, th, expected = TRUE)
    obs <- if (ord == 3L) distrib_deriv3 else distrib_deriv4
    for (nm in names(got)) {
      f <- function(u) {
        yy <- th$a * u^(1 / th$p)
        obs(d, yy, th)[[nm]] * stats::dgamma(u, shape = k)
      }
      ref <- stats::integrate(f, 0, Inf, rel.tol = 1e-12)$value
      expect_equal(got[[nm]], ref, tolerance = 1e-9, label = nm)
    }
  }
})


test_that("the generalized gamma's small-shape coefficients do not cancel", {
  # at small k = d/p, psi(k) + k psi'(k) and 2 psi'(k) + k psi''(k) are
  # differences of terms of size 1/k and 1/k^2; written at k + 1 they are not
  d <- gengamma1_distrib()
  for (k in c(1e-3, 1e-5)) {
    th <- list(a = 1, d = k, p = 1)
    e <- distrib_expected_hessian(d, 1, th)
    expect_equal(e$d_p, digamma(k + 1) + k * trigamma(k + 1), tolerance = 1e-14)
    de <- distrib_dexpected_hessian(d, 1, th)
    expect_equal(de$d_p_d, 2 * trigamma(k + 1) + k * psigamma(k + 1, 2),
                 tolerance = 1e-13)
  }
})


test_that("the generalized gamma's response derivatives are closed", {
  skip_if_not_installed("numDeriv")
  d <- gengamma1_distrib()
  y <- c(0.4, 1.1, 2.3)
  th <- list(a = 1.3, d = 2.1, p = 1.6)
  dy <- function(f, i, t2 = th) numDeriv::grad(function(z) f(d, z, t2), y[i])
  for (i in seq_along(y)) {
    expect_equal(distrib_deriv3_y(d, y, th)[i], dy(distrib_hess_y, i), tolerance = 1e-8)
    expect_equal(distrib_deriv4_y(d, y, th)[i], dy(distrib_deriv3_y, i), tolerance = 1e-8)
    c2 <- distrib_cross2_y(d, y, th)
    for (nm in names(c2)) {
      ref <- numDeriv::grad(function(z) distrib_cross_y(d, z, th)[[nm]], y[i])
      expect_equal(c2[[nm]][i], ref, tolerance = 1e-8, label = nm)
    }
    gh <- distrib_grad_y_hess(d, y, th)
    hh <- distrib_hess_y_hess(d, y, th)
    for (nm in names(gh)) {
      ref <- numDeriv::grad(function(z) distrib_hessian(d, z, th)[[nm]], y[i])
      expect_equal(gh[[nm]][i], ref, tolerance = 1e-8, label = nm)
      ref <- numDeriv::grad(function(z) distrib_grad_y_hess(d, z, th)[[nm]], y[i])
      expect_equal(hh[[nm]][i], ref, tolerance = 1e-8, label = nm)
    }
  }
})


test_that("a power varying by observation reaches the response derivatives", {
  # dy_pow() once reduced the falling factorial of a vector exponent to one
  # number; gengamma1 and weibull1 returned a wrong third and fourth
  # derivative in the response whenever the power varied
  y <- c(0.7, 1.3)
  cases <- list(list(gengamma1_distrib(), list(a = c(1.2, 1.2), d = c(2, 2), p = c(1.5, 3))),
                list(weibull1_distrib(), list(mu = c(1.2, 1.2), sigma = c(1.5, 3))))
  for (cs in cases) {
    for (f in list(distrib_deriv3_y, distrib_deriv4_y)) {
      v <- f(cs[[1]], y, cs[[2]])
      one <- vapply(1:2, function(i) f(cs[[1]], y[i], lapply(cs[[2]], `[`, i)), 0)
      expect_equal(v, one, tolerance = 1e-14)
    }
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


test_that("the generalized gamma's moments keep their digits towards the lognormal", {
  # exact values (300 digits, stabilita/gen_gg_moments_ref.py); formed from the
  # raw moments the variance lost 2.4e-4 at k = 1e6 and the kurtosis all of
  # its digits
  d1 <- gengamma1_distrib(); d2 <- gengamma2_distrib()
  ref <- list(
    list(d = 1e7, p = 10, mean1 = 7.962143052773451806075, var1 = 6.33957476681011702834e-7,
         var2 = 4.00000162000032939987e-8, skew = -0.0007000000607498916877316,
         kurt = 9.60000024299482976703e-7),
    list(d = 1e6, p = 1, mean1 = 2e6, var1 = 4e6, var2 = 4e-6, skew = 0.002, kurt = 6e-6),
    list(d = 2e8, p = 2, mean1 = 19999.99997500000001563, var1 = 0.999999998749999996875,
         var2 = 1.000000001249999996875e-8, skew = 0.00005000000015624999999023,
         kurt = 1.875000014062499956055e-17),
    list(d = 3, p = 10, mean1 = 1.482940592636823050565, var1 = 0.1708193088506930502525,
         var2 = 0.3107058605642900601734, skew = -0.54854248328472726573,
         kurt = -0.1832669087786699989039))
  for (r in ref) {
    t1 <- list(a = 2, d = r$d, p = r$p); t2 <- list(mean = 2, d = r$d, p = r$p)
    expect_equal(mean(d1, t1), r$mean1, tolerance = 1e-14)
    expect_equal(variance(d1, t1), r$var1, tolerance = 1e-14)
    expect_equal(variance(d2, t2), r$var2, tolerance = 1e-14)
    expect_equal(skewness(d1, t1), r$skew, tolerance = 1e-13)
    expect_equal(skewness(d2, t2), r$skew, tolerance = 1e-13)
    # at p = 2 the excess is of order 1/k^2 and ill conditioned in p itself
    expect_equal(kurtosis(d1, t1), r$kurt, tolerance = if (r$p == 2) 1e-5 else 1e-12)
  }
  # vectors of parameters, one row each
  v <- variance(d1, list(a = 2, d = c(1e7, 1e6), p = c(10, 1)))
  expect_equal(v, c(ref[[1]]$var1, ref[[2]]$var1), tolerance = 1e-14)
})
