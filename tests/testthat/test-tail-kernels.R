# The derivatives of the pseudo-Huber families and of the skew t far in the
# tail, where the squared residual overflows. Values from sympy's derivatives
# evaluated by mpmath at 400 digits (the pseudo-Huber families) and from
# mpmath.diff at 150 digits on the skew t's density.

test_that("the pseudo-Huber derivatives stay finite and exact far in the tail", {
  d <- pseudohuber_distrib()
  th <- function(v) list(mu = 0.2, sigma = 1.3, nu = v)
  expect_equal(distrib_hessian(d, -1e40, th(0.7))$mu_mu, -9.1e-121,
               tolerance = 1e-14)
  expect_equal(distrib_deriv3(d, 1e100, th(4))$mu_sigma_nu, -5.0e-201,
               tolerance = 1e-14)
  expect_equal(distrib_deriv4(d, 1e10, th(0.7))$mu_mu_mu_mu,
               -1.0920000001092e-49, tolerance = 1e-14)
  expect_equal(distrib_deriv4(d, 1e200, th(4))$sigma_sigma_sigma_sigma,
               -6.4638977842297053809e+200, tolerance = 1e-14)
  d2 <- pseudohuber2_distrib()
  expect_equal(distrib_gradient(d2, -1e40, th(0.7))$sigma,
               9.4606804505686612028e+39, tolerance = 1e-14)
  # 1 - R zD^2 written as nu/D^2: as a difference this read 2.7e-26
  expect_equal(distrib_hessian(d2, 1e10, th(0.7))$mu_mu,
               -5.6915730460115095221e-31, tolerance = 1e-13)
  expect_equal(distrib_hessian(d2, 1e100, th(4))$mu_sigma,
               -1.1271553380198322748, tolerance = 1e-14)
  expect_equal(distrib_hessian(d2, 1e200, th(4))$sigma_sigma,
               -1.7340851354151265766e+200, tolerance = 1e-14)
})

test_that("the skew t's derivatives stay finite and exact far in the tail", {
  d <- skewt_distrib()
  th <- list(mu = 0.2, sigma = 1.3, alpha = 0.7, nu = 4)
  expect_equal(distrib_deriv3(d, -1e60, th)$sigma_sigma_sigma,
               3.6413290851160673646, tolerance = 1e-13)
  expect_equal(distrib_deriv4(d, 1e30, th)$sigma_sigma_sigma_sigma,
               -8.4030671194986169952, tolerance = 1e-13)
  expect_equal(distrib_deriv4(d, -1e10, th)$alpha_alpha_alpha_alpha,
               -7.7560358760602505738, tolerance = 1e-13)
})

test_that("the heavy- and wide-tailed families' derivatives are finite at any y", {
  ys <- c(-1e300, -1e150, -1e40, 1e40, 1e150, 1e300)
  cases <- list(
    list(pseudohuber_distrib(), list(mu = 0.2, sigma = 1.3, nu = 4)),
    list(pseudohuber2_distrib(), list(mu = 0.2, sigma = 1.3, nu = 4)),
    list(skewt_distrib(), list(mu = 0.2, sigma = 1.3, alpha = 0.5, nu = 4)))
  for (cs in cases) {
    for (g in list(distrib_gradient, distrib_hessian, distrib_deriv3,
                   distrib_deriv4, distrib_deriv5)) {
      out <- g(cs[[1]], ys, lapply(cs[[2]], rep_len, length(ys)))
      expect_true(all(vapply(out, function(v) all(is.finite(v)), logical(1))),
                  label = cs[[1]]@distrib_name)
    }
  }
})
