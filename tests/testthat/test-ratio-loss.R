## log(1+w) - w where 1 + w is KNOWN EXACTLY and forming w loses it.
##
## `psi_Ew(z - 1)` represents 1 + w only to an absolute eps, so its error is
## about eps/(2z) and grows without bound as z -> 0; below z ~ 1.1e-16 the
## subtraction has lost the argument outright, z - 1 is exactly -1, and
## log1p(-1) is -Inf.  `psi_Ew2(z, z - 1)` carries the ratio as well and is
## bounded wherever the value is.  Two families reach the damaged region from
## opposite directions: a gamma of shape below about 0.4, and a negative
## binomial at a zero count with the dispersion at its Poisson limit.

## the form that was replaced, written out so the negative control is real
old_psi_Ew <- function(w) {
  ifelse(abs(w) < 1e-3,
         w * w * (-0.5 + w * (1 / 3 + w * (-0.25 + w * 0.2))),
         log1p(w) - w)
}

test_that("the probe points really are in the damaged region", {
  # Without this the tests below could pass by never reaching the regime they
  # are about.  The old form is EXACTLY -Inf at the smallest probe and finite
  # and wrong at the larger ones.
  expect_identical((1e-17 - 1), -1)                       # the subtraction
  expect_identical(old_psi_Ew(1e-17 - 1), -Inf)
  exact <- function(z) log(z) - z + 1
  expect_gt(abs(old_psi_Ew(1e-14 - 1) - exact(1e-14)) / abs(exact(1e-14)), 1e-6)
  # and a gamma of this shape puts observations there
  set.seed(4)
  y <- distrib_rng(gamma1_distrib(), 20000L, list(mu = 1.19, phi = 4.24))
  expect_gt(sum(y < 1e-15), 0L)
})

test_that("gamma1's dispersion derivatives are finite and right where the ratio vanishes", {
  d <- gamma1_distrib()
  th <- list(mu = 1, phi = 4)
  ld <- function(z, p) distrib_pdf(d, z, list(mu = 1, phi = p), log = TRUE)
  for (z in c(1e-6, 1e-10, 1e-14, 1e-16, 1e-17, 1e-20, 1e-40)) {
    g <- distrib_gradient(d, z, th)$phi
    h <- distrib_hessian(d, z, th)$phi_phi
    expect_true(is.finite(g), info = format(z))
    expect_true(is.finite(h), info = format(z))
    # against a central difference of the ANALYTIC log-density, which shares
    # no arithmetic with the derivative kernels
    s <- 1e-4 * th$phi
    n1 <- (ld(z, th$phi + s) - ld(z, th$phi - s)) / (2 * s)
    n2 <- (ld(z, th$phi + s) - 2 * ld(z, th$phi) + ld(z, th$phi - s)) / s^2
    expect_equal(g, n1, tolerance = 1e-6, info = format(z))
    expect_equal(h, n2, tolerance = 1e-6, info = format(z))
  }
})

test_that("no seed of generate_random_theta leaves gamma1 with a non-finite Hessian", {
  # 13 of these 40 seeds produced -Inf before the repair, always through one
  # observation of twenty thousand.
  d <- gamma1_distrib()
  bad <- 0L
  for (s in 1:12) {
    set.seed(s); th <- generate_random_theta(d)
    y <- distrib_rng(d, 5000L, th)
    H <- distrib_hessian(d, y, th)
    if (!all(is.finite(H$phi_phi))) bad <- bad + 1L
  }
  expect_identical(bad, 0L)
})

test_that("negbin2's dispersion score is exact at the Poisson limit with a zero count", {
  d <- negbin2_distrib()
  # At y = 0 the log-density is theta*log(theta/(theta+mu)), so the score is
  # log(theta) - log(theta+mu) + 1 - theta/(theta+mu), which at mu = 1 and a
  # tiny theta is log(theta) + 1 exactly.
  for (t0 in c(1e-6, 1e-12, 1e-16, 1e-20, 1e-40)) {
    g <- distrib_gradient(d, 0, list(mu = 1, theta = t0))$theta
    ex <- log(t0) - log(t0 + 1) + 1 - t0 / (t0 + 1)
    expect_true(is.finite(g), info = format(t0))
    expect_equal(g, ex, tolerance = 1e-9, info = format(t0))
  }
  # the form that was replaced gives -Inf from 1e-16 down
  expect_identical(old_psi_Ew((0 - 1) / (1e-16 + 1)), -Inf)
})

test_that("the repair costs nothing where the old form was accurate", {
  # Above the derived crossover the two spellings are the same expression, so
  # an ordinary sample is untouched.  The summed score is what a fit reads.
  d <- gamma1_distrib()
  for (shape in c(1, 4)) {
    set.seed(5)
    z <- stats::rgamma(4000L, shape = shape) / shape
    g <- distrib_gradient(d, z, list(mu = 1, phi = 1 / shape))$phi
    expect_true(all(is.finite(g)))
    # every observation of an ordinary sample sits above the crossover
    expect_gt(min(z), 1e-8)
  }
})
