# The derivatives of Sigma^-1 near a singular covariance. Where the
# log-Cholesky chart puts log L22 at -11.5, Sigma^-1 has entries of order
# 1e10 and nearly rank one, and forming its derivatives as the sandwich
# -Sigma^-1 A Sigma^-1 cancels: the second derivative came out entirely
# wrong there. Each order is checked against ONE Richardson difference of the
# order below, the order below reading Sigma^-1 through param_solve() at
# order zero and through the chart's own derivative of the inverse above it.

near_singular_theta <- function(d) {
  v <- c(mu1 = 0.3, mu2 = -0.2, sigma_log_L1 = -0.5, sigma_log_L2 = -11.5,
         sigma_L2.1 = 0.07, nu = 5)
  expect_true(all(d@params %in% names(v)))
  as.list(v[d@params])
}

rel <- function(a, b) max(abs(a - b)) / max(abs(b))

check_family <- function(d, y) {
  th <- near_singular_theta(d)
  mat <- grep("^sigma_", d@params, value = TRUE)
  step <- function(k, f) {
    matrix(numDeriv::jacobian(function(z) {
      t2 <- th
      t2[[k]] <- z
      as.vector(f(t2))
    }, th[[k]]))
  }
  c2 <- distrib_cross2_y(d, y, th)
  for (k in mat) {
    ref <- step(k, function(t2) distrib_hess_y(d, y, t2))
    expect_lt(rel(as.vector(c2[[k]]), as.vector(ref)), 1e-7)
  }
  hh <- distrib_hess_y_hess(d, y, th)
  for (a in mat) for (b in mat) {
    key <- paste(a, b, sep = "_")
    if (!key %in% names(hh)) key <- paste(b, a, sep = "_")
    expect_true(key %in% names(hh))
    ref <- step(b, function(t2) distrib_cross2_y(d, y, t2)[[a]])
    expect_lt(rel(as.vector(hh[[key]]), as.vector(ref)), 1e-6)
  }
}

test_that("the multivariate gaussian's response derivatives hold near a singular covariance", {
  skip_if_not_installed("numDeriv")
  set.seed(1)
  check_family(mvgaussian1_distrib(2), matrix(stats::rnorm(6), ncol = 2))
})

# The multivariate t is not checked here. Its response Hessian carries
# w w' / (nu + q) with q = r' Sigma^-1 r of order 1e10 at this point, which
# cancels against -c Sigma^-1 in its own right, and a Richardson difference
# of distrib_cross2_y() there moves by a factor of two when the order below
# changes in its seventh digit: the reference is not a reference. The t's
# second-order response derivatives near a singular scale matrix are an open
# item, before and after the change that motivated this file.
