# The Student t's derivatives far in the tail, where q = z^2/k overflows and
# t = 1/(1 + q) underflows. Values from sympy's derivatives evaluated by
# mpmath at 400 digits. Written in q and z, the fourth derivatives were NaN
# from |y| of about 1e40.

test_that("the Student t's derivatives stay finite and exact far in the tail", {
  th <- list(mu = 0.2, sigma = 1.3, nu = 4)
  d1 <- student_t1_distrib()
  expect_equal(distrib_gradient(d1, 1e200, th)$mu, 5.0e-200, tolerance = 1e-14)
  expect_equal(distrib_hessian(d1, 1e100, th)$sigma_nu, 0.76923076923076923077,
               tolerance = 1e-14)
  expect_equal(distrib_deriv4(d1, -1e40, th)$sigma_sigma_sigma_sigma,
               -8.4030671194986169952, tolerance = 1e-14)
  d2 <- student_t2_distrib()
  expect_equal(distrib_hessian(d2, -1e300, list(mu = 0.2, sigma = 1.3, nu = 50))$sigma_sigma,
               -29.585798816568047337, tolerance = 1e-14)
  ys <- c(-1e300, -1e150, -1e40, 1e40, 1e150, 1e300)
  for (d in list(d1, d2)) {
    for (g in list(distrib_gradient, distrib_hessian, distrib_deriv3,
                   distrib_deriv4, distrib_deriv5)) {
      out <- g(d, ys, lapply(th, rep_len, length(ys)))
      expect_true(all(vapply(out, function(v) all(is.finite(v)), logical(1))))
    }
  }
})
