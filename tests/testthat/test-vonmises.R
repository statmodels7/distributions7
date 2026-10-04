# The von Mises: the first family here whose support is a circle, so the two
# ends of the interval are the same point and the density need not vanish at
# either.

test_that("the density is the formula and integrates to one", {
  d <- vonmises1_distrib()
  th <- list(mu = 0.5, kappa = 2)
  y <- c(-2, -0.3, 0.5, 1.8, 3)

  expect_equal(distrib_pdf(d, y, th),
               exp(2 * cos(y - 0.5)) / (2 * pi * besselI(2, 0)))
  expect_equal(stats::integrate(function(z) distrib_pdf(d, z, th),
                                -pi, pi)$value, 1, tolerance = 1e-10)

  # the density does not vanish at the ends: they are the same point
  expect_gt(distrib_pdf(d, -pi, th), 0)
  expect_identical(d@bounds, c(-pi, pi))
})


test_that("a large concentration does not overflow the Bessel constant", {
  # I_0(kappa) itself is infinite past about 700; the scaled form with the
  # exponent added back is finite for any concentration.
  d <- vonmises1_distrib()
  expect_true(is.infinite(besselI(900, 0)))
  expect_true(is.finite(distrib_pdf(d, 0.5, list(mu = 0.5, kappa = 900))))
  expect_true(is.finite(distrib_pdf(d, 0.5, list(mu = 0.5, kappa = 5000),
                                    log = TRUE)))
})


test_that("the analytical derivatives match one Richardson differentiation", {
  skip_if_not_installed("numDeriv")
  d <- vonmises1_distrib()
  th <- list(mu = 0.5, kappa = 2)
  p0 <- c(0.5, 2)
  at <- function(p) list(mu = p[1], kappa = p[2])
  set.seed(1)
  y <- distrib_rng(d, 60, th)

  g <- distrib_gradient(d, y, th)
  expect_equal(vapply(g, sum, numeric(1)),
    numDeriv::grad(function(p) sum(distrib_pdf(d, y, at(p), log = TRUE)), p0),
    tolerance = 1e-7, ignore_attr = TRUE)

  J <- numDeriv::jacobian(function(p) {
    vapply(distrib_gradient(d, y, at(p)), sum, numeric(1))
  }, p0)
  h <- distrib_hessian(d, y, th)
  expect_equal(sum(h$mu_mu), J[1, 1], tolerance = 1e-7)
  expect_equal(sum(h$mu_kappa), J[1, 2], tolerance = 1e-7)
  expect_equal(sum(h$kappa_kappa), J[2, 2], tolerance = 1e-7)
})


test_that("the direction and the concentration are orthogonal", {
  # E[sin(Y - mu)] = 0 by symmetry, so the expected information is diagonal
  # and Fisher scoring updates the two independently.
  d <- vonmises1_distrib()
  for (th in list(list(mu = 0.5, kappa = 2), list(mu = -1, kappa = 0.3),
                  list(mu = 2, kappa = 12))) {
    eh <- distrib_expected_hessian(d, 0, th)
    expect_identical(eh$mu_kappa[1], 0)

    # and the diagonal against quadrature, which shares no code with it
    for (nm in c("mu_mu", "kappa_kappa")) {
      q <- stats::integrate(function(z) {
        distrib_hessian(d, z, th)[[nm]] * distrib_pdf(d, z, th)
      }, -pi, pi, rel.tol = 1e-11)$value
      expect_equal(eh[[nm]][1], q, tolerance = 1e-9, label = nm)
    }

    # A'(kappa) is the variance of cos(Y - mu), hence positive
    expect_gt(numericals7::bessel_i_ratio_d1(th$kappa), 0)
  }
})


test_that("the generator reproduces the direction and the resultant length", {
  d <- vonmises1_distrib()
  th <- list(mu = 0.5, kappa = 2)
  set.seed(2)
  y <- distrib_rng(d, 2e5, th)

  expect_true(all(y >= -pi & y < pi + 1e-9))
  expect_equal(atan2(mean(sin(y)), mean(cos(y))), th$mu, tolerance = 0.02)
  expect_equal(sqrt(mean(cos(y))^2 + mean(sin(y))^2), numericals7::bessel_i_ratio(th$kappa),
               tolerance = 0.01)
})


test_that("per-observation parameters give each draw its own mu and kappa", {
  # two groups interleaved, so a draw matched to the wrong observation's
  # parameters would mix the groups' directions and resultant lengths
  n <- 4e4
  g <- rep_len(1:2, n)
  mu <- c(-1, 1.5)[g]
  kappa <- c(0.5, 8)[g]
  direction <- function(y) atan2(mean(sin(y)), mean(cos(y)))

  set.seed(6)
  y1 <- expect_no_warning(
    distrib_rng(vonmises1_distrib(), n, list(mu = mu, kappa = kappa)))
  set.seed(7)
  y2 <- expect_no_warning(
    distrib_rng(vonmises2_distrib(), n,
                list(mu = mu, rho = numericals7::bessel_i_ratio(kappa))))

  for (y in list(y1, y2)) {
    expect_length(y, n)
    expect_true(all(y >= -pi & y < pi))
    for (j in 1:2) {
      expect_equal(direction(y[g == j]), c(-1, 1.5)[j], tolerance = 0.05)
      # mean of cos(y - mu) is unbiased for A(kappa); compared within four
      # standard errors, since at kappa = 0.5 a relative tolerance on a
      # resultant near 0.24 is about one standard error
      cy <- cos(y[g == j] - c(-1, 1.5)[j])
      expect_lt(abs(mean(cy) - numericals7::bessel_i_ratio(c(0.5, 8)[j])),
                4 * stats::sd(cy) / sqrt(length(cy)))
    }
  }
})


test_that("the validator passes and a fit recovers the parameters", {
  d <- vonmises1_distrib()
  set.seed(3)
  res <- check_distrib(d, verbose = FALSE)
  expect_true(all(res$status == "OK"),
    label = paste(res$check[res$status != "OK"], collapse = ", "))

  set.seed(5)
  th <- list(mu = 0.5, kappa = 2)
  y <- distrib_rng(d, 3000, th)
  f <- fit_distrib(d, y)
  expect_true(f@converged, info = fit_report(f, d, th))
  expect_equal(unname(coef(f)), c(0.5, 2), tolerance = 0.15)
})


test_that("the third and fourth derivatives are closed for both parametrizations", {
  skip_if_not_installed("numDeriv")
  # each order against ONE Richardson differentiation of the analytic order
  # below it -- never a chain of plain differences, which by order 3 would be
  # noise
  y <- c(-2.0, -0.4, 0.7, 2.6)
  chk <- function(d, th, order, tol) {
    lower <- if (order == 3L) distrib_hessian else distrib_deriv3
    got <- if (order == 3L) distrib_deriv3(d, y, th) else distrib_deriv4(d, y, th)
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
      expect_equal(got[[nm]], ref, tolerance = tol,
                   label = sprintf("%s %s", d@distrib_name, nm))
    }
  }
  for (th in list(list(mu = 0.3, kappa = 1.7), list(mu = -1.2, kappa = 0.2),
                  list(mu = 2.0, kappa = 12))) {
    chk(vonmises1_distrib(), th, 3L, 1e-6)
    chk(vonmises1_distrib(), th, 4L, 1e-6)
  }
  for (th in list(list(mu = 0.3, rho = 0.6), list(mu = -1.2, rho = 0.15),
                  list(mu = 2.0, rho = 0.93))) {
    chk(vonmises2_distrib(), th, 3L, 1e-6)
    chk(vonmises2_distrib(), th, 4L, 1e-6)
  }
})


test_that("the components mixing mu with two concentrations vanish exactly", {
  # kappa enters kappa*cos(y - mu) linearly, so it is a structural zero and
  # not a small number: asserted as such rather than to a tolerance
  d <- vonmises1_distrib()
  y <- c(-1, 0.5, 2)
  th <- list(mu = 0.4, kappa = 2.2)
  expect_identical(distrib_deriv3(d, y, th)$mu_kappa_kappa, rep(0, 3))
  d4 <- distrib_deriv4(d, y, th)
  expect_identical(d4$mu_mu_kappa_kappa, rep(0, 3))
  expect_identical(d4$mu_kappa_kappa_kappa, rep(0, 3))
})

test_that("the log-density is finite at concentrations past the scaled-Bessel underflow", {
  # log(besselI(k, 0, expon.scaled = TRUE)) is -Inf between kappa = 1e5 and
  # 1e6, where the scaled Bessel underflows to an exact zero;
  # numericals7::log_bessel_i is finite there, and the log-density must be too.
  d <- vonmises1_distrib()
  for (k in c(1e5, 5e5, 1e6)) {
    lp <- distrib_pdf(d, c(-0.2, 0, 0.3), list(mu = 0, kappa = k), log = TRUE)
    expect_true(all(is.finite(lp)))
  }
  # against the direct scaled route where that route still works
  k <- 50
  expect_equal(
    distrib_pdf(d, 0.4, list(mu = 0.1, kappa = k), log = TRUE),
    k * cos(0.4 - 0.1) - log(2 * pi) - (log(besselI(k, 0, expon.scaled = TRUE)) + k),
    tolerance = 1e-12
  )
})

test_that("the series distribution function agrees with the quadrature it replaces", {
  # the base class's numerical integration of the density, which is what the
  # method replaces and which shares no arithmetic with the Bessel series
  quad <- S7::method(distrib_cdf, continuous_distrib)
  d <- vonmises1_distrib()
  x <- seq(-pi, pi - 1e-9, length.out = 25)
  for (kap in c(0.05, 0.5, 2, 10, 50, 200)) {
    for (mu in c(0, 1.2, -2.9)) {
      th <- list(mu = mu, kappa = kap)
      expect_equal(distrib_cdf(d, x, th), quad(d, x, th), tolerance = 1e-11,
                   info = sprintf("kappa %g, mu %g", kap, mu))
    }
  }
})

test_that("the series has the properties a distribution function has", {
  d <- vonmises1_distrib()
  th <- list(mu = 0.9, kappa = 7)
  expect_equal(distrib_cdf(d, -pi, th), 0)
  expect_equal(distrib_cdf(d, pi - 1e-12, th), 1, tolerance = 1e-9)
  x <- seq(-pi, pi - 1e-12, length.out = 400)
  f <- distrib_cdf(d, x, th)
  expect_true(all(diff(f) >= -1e-14))
  expect_true(all(f >= 0 & f <= 1))
  # outside the support it is the constant it must be, not an extrapolation
  expect_equal(distrib_cdf(d, -4, th), 0)
  expect_equal(distrib_cdf(d, 4, th), 1)
  # and its derivative is the density
  h <- 1e-5
  inn <- x > -pi + 0.01 & x < pi - 0.01
  num <- (distrib_cdf(d, x[inn] + h, th) - distrib_cdf(d, x[inn] - h, th)) /
    (2 * h)
  expect_equal(num, distrib_pdf(d, x[inn], th), tolerance = 1e-8)
})

test_that("the term count is a safe upper bound wherever it is used", {
  # the rule is 8.5 sqrt(kappa) + 10; the check is that four times as many
  # terms change nothing, which is what says the series has converged
  d <- vonmises1_distrib()
  x <- seq(-pi + 1e-9, pi - 1e-9, length.out = 21)
  long <- function(mu, kappa, m) {
    R <- numericals7::bessel_i_ratios(rep(kappa, length(x)), m)
    j <- seq_len(m)
    s <- rowSums(R * (sin(outer(x - mu, j)) +
                        sin(outer(rep(pi + mu, length(x)), j))) /
                   rep(j, each = length(x)))
    (x + pi) / (2 * pi) + s / pi
  }
  for (kap in c(0.5, 5, 50, 300)) {
    m <- max(20L, as.integer(ceiling(8.5 * sqrt(kap))) + 10L)
    expect_equal(distrib_cdf(d, x, list(mu = 0.7, kappa = kap)),
                 long(0.7, kap, 4L * m), tolerance = 1e-13,
                 info = sprintf("kappa %g", kap))
  }
})

test_that("the resultant-length family is the concentration one at its kappa", {
  # the map touches the second parameter only and the response not at all
  d1 <- vonmises1_distrib()
  d2 <- vonmises2_distrib()
  x <- seq(-pi, pi - 1e-9, length.out = 21)
  for (rho in c(0.05, 0.3, 0.7, 0.95)) {
    k <- numericals7::bessel_i_ratio_inverse(rho)
    expect_identical(distrib_cdf(d2, x, list(mu = 0.5, rho = rho)),
                     distrib_cdf(d1, x, list(mu = 0.5, kappa = k)))
  }
  quad <- S7::method(distrib_cdf, continuous_distrib)
  th <- list(mu = -2.9, rho = 0.7)
  expect_equal(distrib_cdf(d2, x, th), quad(d2, x, th), tolerance = 1e-11)
})

test_that("a concentration varying by observation is read per row", {
  d <- vonmises1_distrib()
  x <- seq(-pi + 0.1, pi - 0.1, length.out = 6)
  kk <- c(0.2, 1, 3, 8, 20, 60)
  got <- distrib_cdf(d, x, list(mu = 0.3, kappa = kk))
  one <- vapply(seq_along(x), function(i)
    distrib_cdf(d, x[i], list(mu = 0.3, kappa = kk[i])), numeric(1))
  expect_equal(got, one, tolerance = 1e-12)
})

test_that("a scalar and a per-observation parameter give the same bits", {
  # the kernels evaluate the Bessel ratio once for a scalar parameter and
  # per observation otherwise; the two paths must agree exactly
  y <- c(-2.5, -1, 0, 0.5, 2, 3)
  n <- length(y)
  cases <- list(list(vonmises1_distrib(), list(mu = 0.4, kappa = 3.7)),
                list(vonmises1_distrib(), list(mu = 0.4, kappa = 1e-6)),
                list(vonmises2_distrib(), list(mu = 0.4, rho = 0.83)),
                list(vonmises2_distrib(), list(mu = 0.4, rho = 1e-4)))
  for (cs in cases) {
    d <- cs[[1]]; th <- cs[[2]]
    thv <- list(th[[1]], rep(th[[2]], n))
    names(thv) <- names(th)
    for (f in list(distrib_gradient, distrib_hessian, distrib_deriv3,
                   distrib_deriv4, distrib_expected_hessian)) {
      expect_identical(f(d, y, th), f(d, y, thv))
      expect_identical(f(d, y, thv, threads = 2L), f(d, y, thv))
    }
  }
})

test_that("the resultant-length derivatives match the chain rule written out", {
  # the kernels use d log I0(kappa(rho)) / d rho = rho kappa'(rho); the
  # Faa di Bruno expansion through A and the inverse's derivatives is the
  # same quantity by other arithmetic
  d2 <- vonmises2_distrib()
  y <- c(-2, 0.1, 1.3)
  for (rho in c(1e-4, 0.3, 0.9, 0.999)) {
    th <- list(mu = 0.2, rho = rho)
    k <- numericals7::bessel_i_ratio_inverse(rho)
    A <- numericals7::bessel_i_ratio(k)
    a1 <- numericals7::bessel_i_ratio_d1(k)
    a2 <- numericals7::bessel_i_ratio_d2(k)
    a3 <- numericals7::bessel_i_ratio_d3(k)
    k1 <- numericals7::bessel_i_ratio_inverse_d1(k)
    k2 <- numericals7::bessel_i_ratio_inverse_d2(k)
    k3 <- numericals7::bessel_i_ratio_inverse_d3(k)
    k4 <- numericals7::bessel_i_ratio_inverse_d4(k)
    cc <- cos(y - 0.2)
    phi3 <- a2 * k1^3 + 3 * a1 * k1 * k2 + A * k3
    phi4 <- a3 * k1^4 + 6 * a2 * k1^2 * k2 + a1 * (3 * k2^2 + 4 * k1 * k3) +
      A * k4
    scale3 <- max(abs(cc * k3), abs(phi3))
    scale4 <- max(abs(cc * k4), abs(phi4))
    expect_lt(max(abs(distrib_deriv3(d2, y, th)$rho_rho_rho -
                        (cc * k3 - phi3))) / scale3, 1e-13)
    expect_lt(max(abs(distrib_deriv4(d2, y, th)$rho_rho_rho_rho -
                        (cc * k4 - phi4))) / scale4, 1e-13)
  }
})
