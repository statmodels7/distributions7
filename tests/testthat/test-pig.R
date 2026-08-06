# The Poisson-inverse Gaussian, both parametrizations. References:
# gamlss.dist (an independent implementation, recursion-based where ours is
# the closed half-integer Bessel sum), one numDeriv pass per order on the
# analytic order below, the reparametrize() twin, and direct series.

test_that("pig1 matches gamlss's PIG", {
  skip_if_not_installed("gamlss.dist")
  d <- pig1_distrib()
  y <- c(0, 1, 2, 5, 17, 60, 200)
  for (th in list(c(1, 1), c(0.3, 4), c(10, 0.05), c(50, 2), c(2, 0.001))) {
    ref <- gamlss.dist::dPIG(y, mu = th[1], sigma = th[2], log = TRUE)
    got <- distrib_pdf(d, y, list(mu = th[1], sigma = th[2]), log = TRUE)
    expect_equal(got, ref, tolerance = 1e-12,
                 label = sprintf("pig1 at mu=%g sigma=%g", th[1], th[2]))
  }
})


test_that("pig2 matches gamlss's PIG2, whose sigma is the Bessel argument", {
  skip_if_not_installed("gamlss.dist")
  d <- pig2_distrib()
  y <- c(0, 1, 2, 5, 17, 60)
  for (th in list(c(3, 1.2), c(0.5, 0.1), c(20, 8))) {
    ref <- gamlss.dist::dPIG2(y, mu = th[1], sigma = th[2], log = TRUE)
    got <- distrib_pdf(d, y, list(mu = th[1], alpha = th[2]), log = TRUE)
    expect_equal(got, ref, tolerance = 1e-12,
                 label = sprintf("pig2 at mu=%g alpha=%g", th[1], th[2]))
  }
})


test_that("the mass is normalized and non-integer y has none", {
  d <- pig1_distrib()
  th <- list(mu = 4, sigma = 0.8)
  expect_equal(sum(distrib_pdf(d, 0:400, th)), 1, tolerance = 1e-12)
  expect_identical(distrib_pdf(d, c(1.5, -1), th), c(0, 0))
  expect_identical(distrib_pdf(d, c(1.5, -1), th, log = TRUE), c(-Inf, -Inf))
})


test_that("the kernel derivatives match one numerical pass per order", {
  skip_if_not_installed("numDeriv")
  d <- pig1_distrib()
  th <- list(mu = 3.2, sigma = 0.7)
  yv <- c(0, 1, 4, 11)

  g <- distrib_gradient(d, yv, th)
  for (p in c("mu", "sigma")) {
    nd <- vapply(yv, function(yy) numDeriv::grad(function(z) {
      tt <- th; tt[[p]] <- z
      distrib_pdf(d, yy, tt, log = TRUE)
    }, th[[p]]), numeric(1))
    expect_equal(g[[p]], nd, tolerance = 1e-6, label = paste("score", p))
  }

  h <- distrib_hessian(d, yv, th)
  nd <- vapply(yv, function(yy) numDeriv::grad(function(z)
    distrib_gradient(d, yy, list(mu = z, sigma = th$sigma))$mu, th$mu),
    numeric(1))
  expect_equal(h$mu_mu, nd, tolerance = 1e-6)
  nd <- vapply(yv, function(yy) numDeriv::grad(function(z)
    distrib_gradient(d, yy, list(mu = th$mu, sigma = z))$mu, th$sigma),
    numeric(1))
  expect_equal(h$mu_sigma, nd, tolerance = 1e-6)

  d3 <- distrib_deriv3(d, yv, th)
  nd <- vapply(yv, function(yy) numDeriv::grad(function(z)
    distrib_hessian(d, yy, list(mu = z, sigma = th$sigma))$mu_mu, th$mu),
    numeric(1))
  expect_equal(d3$mu_mu_mu, nd, tolerance = 1e-5)

  d4 <- distrib_deriv4(d, yv, th)
  nd <- vapply(yv, function(yy) numDeriv::grad(function(z)
    distrib_deriv3(d, yy, list(mu = th$mu, sigma = z))$mu_sigma_sigma,
    th$sigma), numeric(1))
  expect_equal(d4$mu_sigma_sigma_sigma, nd, tolerance = 1e-5)
})


test_that("pig2 agrees with reparametrize(pig1) at every order", {
  # two implementations of one object: the hand-written kernel with alpha a
  # seed variable, against the generic partition-sum wrapper over the map
  # sigma = (mu + sqrt(mu^2 + alpha^2)) / alpha^2
  d2 <- pig2_distrib()
  rp <- reparametrize(
    pig1_distrib(),
    map = function(psi) list(mu = psi$mu, sigma = pig2_sigma(psi$mu, psi$alpha)),
    params = c("mu", "alpha"),
    bounds = list(mu = c(0, Inf), alpha = c(0, Inf)),
    links = list(mu = log_link(), alpha = log_link())
  )
  th <- list(mu = 3.2, alpha = 1.4)
  yv <- c(0, 1, 4, 11)
  for (fn in list(distrib_gradient, distrib_hessian, distrib_deriv3,
                  distrib_deriv4)) {
    a <- fn(d2, yv, th)
    b <- fn(rp, yv, th)
    for (nm in names(b)) {
      expect_equal(a[[nm]], b[[nm]], tolerance = 1e-10, label = nm)
    }
  }
})


test_that("pig2's parameters are orthogonal and pig1's are not", {
  d2 <- pig2_distrib()
  eh2 <- distrib_expected_hessian(d2, 0, list(mu = 3.2, alpha = 1.4))
  expect_lt(abs(eh2$mu_alpha[1]), 1e-8)
  d1 <- pig1_distrib()
  eh1 <- distrib_expected_hessian(d1, 0, list(mu = 3.2, sigma = 0.7))
  expect_gt(abs(eh1$mu_sigma[1]), 1e-3)
})


test_that("the closed moments match direct series", {
  d <- pig1_distrib()
  th <- list(mu = 4, sigma = 0.8)
  pk <- distrib_pdf(d, 0:600, th)
  m <- sum((0:600) * pk)
  v <- sum(((0:600) - m)^2 * pk)
  expect_equal(mean(d, th), m, tolerance = 1e-10)
  expect_equal(variance(d, th), v, tolerance = 1e-9)
  expect_equal(skewness(d, th), sum(((0:600) - m)^3 * pk) / v^1.5,
               tolerance = 1e-8)
  expect_equal(kurtosis(d, th), sum(((0:600) - m)^4 * pk) / v^2 - 3,
               tolerance = 1e-7)

  # pig2's moments are pig1's at the implied dispersion
  a <- sqrt(1 + 2 * 0.8 * 4) / 0.8
  expect_equal(variance(pig2_distrib(), list(mu = 4, alpha = a)),
               variance(d, th), tolerance = 1e-12)
})


test_that("the validator passes both parametrizations", {
  set.seed(31)
  res <- check_distrib(pig1_distrib(), list(mu = 3, sigma = 0.6),
                       verbose = FALSE)
  expect_true(all(res$status == "OK"),
              info = paste(res$check[res$status != "OK"], collapse = ", "))
  set.seed(32)
  res <- check_distrib(pig2_distrib(), list(mu = 3, alpha = 1.5),
                       verbose = FALSE)
  expect_true(all(res$status == "OK"),
              info = paste(res$check[res$status != "OK"], collapse = ", "))
})


test_that("a fit recovers the parameters", {
  set.seed(33)
  d <- pig1_distrib()
  y <- distrib_rng(d, 4000, list(mu = 3, sigma = 0.8))
  f <- fit_distrib(d, y)
  expect_true(f@converged)
  expect_equal(unname(coef(f)), c(3, 0.8), tolerance = 0.15)
})


test_that("the explicit kernels agree with the jet kernels", {
  # Two implementations sharing no algebra: the hand-written closed forms
  # the methods run, against the bivariate-jet transcription kept for this
  # comparison. Measured 2x to 36x faster, the explicit route is what
  # ships; agreement is machine precision away from the sigma -> 0 corner,
  # where terms of order sigma^-5 cancel in BOTH implementations and the
  # comparison measures conditioning, not correctness.
  set.seed(41)
  n <- 3000
  y <- distrib_rng(pig1_distrib(), n, list(mu = 3, sigma = 0.8))
  mu <- runif(n, 0.2, 30)
  sg <- runif(n, 0.2, 4)
  A <- pig1_hd_jet_cpp(y, mu, sg)
  B <- pig1_hd_cpp(y, mu, sg)
  expect_lt(max(abs(A - B) / pmax(1, abs(A))), 1e-10)

  al <- runif(n, 0.2, 8)
  A2 <- pig2_hd_jet_cpp(y, mu, al)
  B2 <- pig2_hd_cpp(y, mu, al)
  expect_lt(max(abs(A2 - B2) / pmax(1, abs(A2))), 1e-10)

  # and the small-sigma corner stays within its conditioning floor
  sg2 <- runif(500, 0.02, 0.2)
  A3 <- pig1_hd_jet_cpp(y[1:500], mu[1:500], sg2)
  B3 <- pig1_hd_cpp(y[1:500], mu[1:500], sg2)
  expect_lt(max(abs(A3 - B3) / pmax(1, abs(A3))), 1e-6)
})
