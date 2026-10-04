# The generalized Pareto. Two things distinguish it from every other family
# here: it is not written in its mean, the mean existing only for xi < 1, and
# for xi < 0 its support depends on the parameters, which is what the license
# to differentiate under the integral sign rests on.

test_that("the density and the distribution function are the written formulas", {
  d <- gpd_distrib()
  s <- 1.4
  y <- c(0.1, 0.5, 1.2)
  for (x in c(-0.3, 0.25, 0.8)) {
    th <- list(sigma = s, xi = x)
    expect_equal(distrib_pdf(d, y, th),
                 (1 / s) * (1 + x * y / s)^(-1 / x - 1), label = as.character(x))
    expect_equal(distrib_cdf(d, y, th),
                 1 - (1 + x * y / s)^(-1 / x), label = as.character(x))
  }

  # the round trip through the quantile function
  for (x in c(-0.3, 0, 0.6)) {
    th <- list(sigma = s, xi = x)
    p <- c(0.1, 0.5, 0.9)
    expect_equal(distrib_cdf(d, distrib_quantile(d, p, th), th), p,
                 tolerance = 1e-10)
  }
})


test_that("the exponential is the limit at zero, reached by a series", {
  d <- gpd_distrib()
  s <- 1.4
  y <- c(0.1, 0.5, 1.2, 4)
  expect_equal(distrib_pdf(d, y, list(sigma = s, xi = 0)),
               stats::dexp(y, rate = 1 / s))
  # and the approach is continuous rather than a special case bolted on
  for (x in c(1e-3, 1e-6, 1e-10)) {
    expect_equal(distrib_pdf(d, y, list(sigma = s, xi = x)),
                 stats::dexp(y, rate = 1 / s), tolerance = 5 * x + 1e-9)
  }
})


test_that("the density integrates to one on both sides of zero", {
  d <- gpd_distrib()
  s <- 1.4
  for (x in c(-0.4, 0, 0.3, 0.9)) {
    up <- if (x < 0) -s / x else Inf
    v <- stats::integrate(function(z) distrib_pdf(d, z, list(sigma = s, xi = x)),
                          0, up, rel.tol = 1e-10)$value
    expect_equal(v, 1, tolerance = 1e-8, label = as.character(x))
  }
})


test_that("the analytical derivatives match one Richardson differentiation", {
  skip_if_not_installed("numDeriv")
  d <- gpd_distrib()
  s <- 1.4
  set.seed(1)
  for (x in c(-0.3, 0.25, 0.8)) {
    th <- list(sigma = s, xi = x)
    y <- distrib_rng(d, 60, th)
    at <- function(p) list(sigma = p[1], xi = p[2])

    g <- distrib_gradient(d, y, th)
    expect_equal(vapply(g, sum, numeric(1)),
      numDeriv::grad(function(p) sum(distrib_pdf(d, y, at(p), log = TRUE)),
                     c(s, x)),
      tolerance = 1e-6, ignore_attr = TRUE, label = as.character(x))

    J <- numDeriv::jacobian(function(p) {
      vapply(distrib_gradient(d, y, at(p)), sum, numeric(1))
    }, c(s, x))
    h <- distrib_hessian(d, y, th)
    expect_equal(sum(h$sigma_sigma), J[1, 1], tolerance = 1e-5)
    expect_equal(sum(h$sigma_xi), J[1, 2], tolerance = 1e-5)
    expect_equal(sum(h$xi_xi), J[2, 2], tolerance = 1e-5)
  }
})


test_that("the expected information is the closed form, checked by quadrature", {
  # Integrating on the PROBABILITY scale rather than on y: for xi < 0 the
  # second derivative blows up at the upper endpoint, and the transform turns
  # that endpoint into an ordinary point. A Monte Carlo average would agree
  # only to a few per cent there, converging slowly for the same reason.
  d <- gpd_distrib()
  s <- 1.4
  for (x in c(-0.3, -0.1, 0.25, 0.8, 1.5)) {
    th <- list(sigma = s, xi = x)
    q <- vapply(c("sigma_sigma", "sigma_xi", "xi_xi"), function(k) {
      stats::integrate(function(u) {
        distrib_hessian(d, distrib_quantile(d, u, th), th)[[k]]
      }, 0, 1, rel.tol = 1e-8, subdivisions = 8000L)$value
    }, numeric(1))
    got <- vapply(distrib_expected_hessian(d, 0, th)[names(q)], function(v) v[1],
                  numeric(1))
    expect_equal(got, q, tolerance = 1e-7, label = as.character(x))
  }
})


test_that("the information does not exist below xi = -1/2 and says so", {
  # The condition is exactly integrability: on the probability scale the
  # second derivative grows like (1-u)^(-2|xi|), integrable if and only if
  # |xi| < 1/2. Returning NA is the honest answer.
  d <- gpd_distrib()
  eh <- distrib_expected_hessian(d, 0, list(sigma = 1.4, xi = -0.7))
  expect_true(all(is.na(unlist(eh))))

  eh2 <- distrib_expected_hessian(d, 0, list(sigma = 1.4, xi = -0.5))
  expect_true(all(is.na(unlist(eh2))))

  # just above the boundary it exists again
  eh3 <- distrib_expected_hessian(d, 0, list(sigma = 1.4, xi = -0.49))
  expect_true(all(is.finite(unlist(eh3))))
})


test_that("the support moves with the parameters when the shape is negative", {
  d <- gpd_distrib()
  s <- 1.4
  x <- -0.4
  endpoint <- -s / x
  expect_equal(endpoint, 3.5)

  expect_gt(distrib_pdf(d, endpoint - 1e-6, list(sigma = s, xi = x)), 0)
  expect_identical(distrib_pdf(d, endpoint + 1e-6, list(sigma = s, xi = x)), 0)
  expect_equal(distrib_cdf(d, endpoint, list(sigma = s, xi = x)), 1)

  # the declared bounds cannot express it, being fixed at construction
  expect_identical(d@bounds, c(0, Inf))
})


test_that("a fit recovers the parameters in the regular regime", {
  d <- gpd_distrib()
  set.seed(9)
  th <- list(sigma = 1.5, xi = 0.3)
  y <- distrib_rng(d, 4000, th)
  f <- fit_distrib(d, y, start = list(sigma = 1, xi = 0.1))
  expect_true(f@converged, info = fit_report(f, d, th))
  expect_equal(unname(coef(f)), c(1.5, 0.3), tolerance = 0.2)
})


test_that("the generalized Pareto's assembly reproduces the compiled kernel", {
  # Running the order 3 and 4 assembly at orders 1 and 2, where a kernel
  # exists, is what licenses it at the orders where none does. Both routes
  # through the shape direction are exercised: the Leibniz form where xi*z is
  # ordinary, the series where it is small.
  d <- gpd_distrib()
  y <- c(0.2, 1.1, 3.5)
  for (th in list(list(sigma = 1.2, xi = 0.35), list(sigma = 2.0, xi = -0.25),
                  list(sigma = 0.7, xi = 1.8), list(sigma = 20, xi = 0.05),
                  list(sigma = 1.0, xi = 1e-7))) {
    g <- distributions7:::gpd_components(y, th, 1L)
    r <- distrib_gradient(d, y, th)
    for (nm in names(r)) expect_equal(g[[nm]], r[[nm]], tolerance = 1e-11)
    h <- distributions7:::gpd_components(y, th, 2L)
    r <- distrib_hessian(d, y, th)
    for (nm in names(r)) expect_equal(h[[nm]], r[[nm]], tolerance = 1e-12)
    for (ord in 3:4) {
      a <- distributions7:::gpd_components(y, th, ord)
      r <- if (ord == 3L) distrib_deriv3(d, y, th) else distrib_deriv4(d, y, th)
      expect_setequal(names(r), names(a))
      for (nm in names(r)) expect_equal(r[[nm]], a[[nm]], tolerance = 1e-10, label = nm)
    }
  }
})


test_that("the two routes through the shape agree where both are accurate", {
  # This is what pins the switch. At xi*z in [0.1, 0.35] the Leibniz form has
  # not yet lost digits to cancellation and the series has long converged, so
  # forcing each route must give the same numbers.
  th <- list(sigma = 1.0, xi = 0.3)
  y <- c(0.4, 0.8, 1.15)
  for (ord in 1:4) {
    a <- distributions7:::gpd_components(y, th, ord, cut = 0)    # Leibniz
    b <- distributions7:::gpd_components(y, th, ord, cut = 10)   # series
    for (nm in names(a)) {
      expect_equal(a[[nm]], b[[nm]], tolerance = 1e-9,
                   label = sprintf("order %d component %s", ord, nm))
    }
  }
})


test_that("the generalized Pareto is closed at third and fourth order", {
  skip_if_not_installed("numDeriv")
  d <- gpd_distrib()
  y <- c(0.2, 1.1, 3.5)
  # Richardson on the analytic order below, at shapes where the kernel it
  # differentiates is itself well conditioned
  for (th in list(list(sigma = 1.2, xi = 0.35), list(sigma = 2.0, xi = -0.25),
                  list(sigma = 0.7, xi = 1.8))) {
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
  }
})


test_that("the generalized Pareto's fifth order is the derivative of its fourth", {
  skip_if_not_installed("numDeriv")
  d <- gpd_distrib()
  y <- c(0.2, 1.1, 3.5)
  for (th in list(list(sigma = 1.2, xi = 0.35), list(sigma = 2.0, xi = -0.25))) {
    got <- distrib_deriv5(d, y, th)
    expect_length(got, 6L)
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
  }
})


test_that("the generalized Pareto's expected higher derivatives are the integrals", {
  # on the probability scale, as for the information; they exist above
  # xi = -1/3 and -1/4 and are NA at and below
  d <- gpd_distrib()
  for (x in c(-0.1, 0.25, 1.5)) {
    th <- list(sigma = 1.4, xi = x)
    for (ord in 3:4) {
      got <- if (ord == 3L) distrib_deriv3(d, 1, th, expected = TRUE)
             else distrib_deriv4(d, 1, th, expected = TRUE)
      obs <- if (ord == 3L) distrib_deriv3 else distrib_deriv4
      for (nm in names(got)) {
        ref <- stats::integrate(function(u) obs(d, distrib_quantile(d, u, th), th)[[nm]],
                                0, 1, rel.tol = 1e-10, subdivisions = 8000L)$value
        expect_equal(got[[nm]], ref, tolerance = 1e-7, label = paste(x, nm))
      }
    }
  }
  e3 <- distrib_deriv3(d, 1, list(sigma = 1.4, xi = c(-0.34, -1 / 3, -0.32)), expected = TRUE)
  expect_equal(is.na(e3$xi_xi_xi), c(TRUE, TRUE, FALSE))
  e4 <- distrib_deriv4(d, 1, list(sigma = 1.4, xi = c(-0.26, -0.25, -0.24)), expected = TRUE)
  expect_equal(is.na(e4$xi_xi_xi_xi), c(TRUE, TRUE, FALSE))
})


test_that("the generalized Pareto keeps its accuracy at the end of the support", {
  # y within a relative 1e-12 of sigma/|xi|: formed from the rounded product
  # xi*y, t lost 1e-4; exact values from stabilita/gen_gpd_ref.py
  d <- gpd_distrib()
  th <- list(sigma = 1.5, xi = -0.7)
  y <- 2.142857142855
  g <- distrib_gradient(d, y, th)
  expect_equal(g$sigma, 285704595153.8072843052, tolerance = 1e-13)
  expect_equal(g$xi, 612224132417.4830211302, tolerance = 1e-13)
  h <- distrib_hessian(d, y, th)
  expect_equal(h$xi_xi, -8.745762395554034366918e+23, tolerance = 1e-13)
  # outside the support every derivative of the log-density is NaN
  expect_true(all(is.nan(unlist(distrib_deriv3(d, c(-1, 3), th)))))
})


test_that("the generalized Pareto's response derivatives are closed", {
  skip_if_not_installed("numDeriv")
  d <- gpd_distrib()
  y <- c(0.2, 1.1, 3.5)
  for (th in list(list(sigma = 1.2, xi = 0.35), list(sigma = 2.0, xi = -0.25),
                  list(sigma = 1.5, xi = 0))) {
    for (i in seq_along(y)) {
      dy <- function(f) numDeriv::grad(function(z) f(d, z, th), y[i])
      expect_equal(distrib_deriv3_y(d, y, th)[i], dy(distrib_hess_y), tolerance = 1e-8)
      expect_equal(distrib_deriv4_y(d, y, th)[i], dy(distrib_deriv3_y), tolerance = 1e-8)
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
  }
})


test_that("the generalized Pareto's four cdf surfaces agree with one another", {
  # log F, log S, S and F from their own kernels: d log S = dS / S, and
  # d log F = dF / F, at first order; the second orders through the product
  d <- gpd_distrib()
  q <- c(0.3, 2, 8)
  for (th in list(list(sigma = 1.5, xi = 0.3), list(sigma = 1.5, xi = -0.1),
                  list(sigma = 1.5, xi = 0))) {
    S <- distrib_cdf(d, q, th, lower.tail = FALSE)
    F <- distrib_cdf(d, q, th)
    gS <- distrib_grad_cdf(d, q, th, lower.tail = FALSE, log = FALSE)
    gF <- distrib_grad_cdf(d, q, th, log = FALSE)
    lS <- distrib_grad_cdf(d, q, th, lower.tail = FALSE)
    lF <- distrib_grad_cdf(d, q, th)
    hS <- distrib_hess_cdf(d, q, th, lower.tail = FALSE, log = FALSE)
    hlS <- distrib_hess_cdf(d, q, th, lower.tail = FALSE)
    for (nm in names(gS)) {
      expect_equal(gF[[nm]], -gS[[nm]], tolerance = 1e-15)
      expect_equal(lS[[nm]], gS[[nm]] / S, tolerance = 1e-13)
      expect_equal(lF[[nm]], gF[[nm]] / F, tolerance = 1e-13)
    }
    expect_equal(hlS$sigma_xi, hS$sigma_xi / S - gS$sigma * gS$xi / S^2,
                 tolerance = 1e-12)
  }
  # below zero log F is -Inf: NaN; past the upper end every derivative is 0
  th <- list(sigma = 1, xi = -0.5)
  expect_true(all(is.nan(distrib_grad_cdf(d, c(-1, 0), th)$sigma)))
  expect_equal(distrib_deriv4_cdf(d, c(2, 3), th)$xi_xi_xi_xi, c(0, 0))
  expect_equal(distrib_deriv4_cdf(d, c(-1, 2), th, lower.tail = FALSE)$sigma_sigma_sigma_sigma,
               c(0, 0))
})
