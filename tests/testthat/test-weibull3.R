# weibull3 written out, against the same law obtained through
# reparametrize() on lognormal1 (weibull3_by_reparam(), an independent
# implementation) and against numDeriv where the reference stops.

w3_cases <- list(list(mean = 4, sigma = 1.7), list(mean = 0.4, sigma = 0.6),
                 list(mean = 60, sigma = 6))
w3_y <- c(0.3, 1.1, 2.7, 8)

expect_same_list <- function(a, b, tol, label) {
  expect_setequal(names(a), names(b))
  for (nm in names(b)) {
    expect_equal(rep_len(a[[nm]], length(b[[nm]])), b[[nm]], tolerance = tol,
                 label = paste(label, nm))
  }
}

test_that("weibull3 is the reparametrized weibull1 on every surface", {
  d <- weibull3_distrib()
  r <- weibull3_by_reparam()
  for (th in w3_cases) {
    y <- w3_y * th$mean
    expect_equal(distrib_pdf(d, y, th, log = TRUE), distrib_pdf(r, y, th, log = TRUE),
                 tolerance = 1e-12)
    expect_equal(distrib_cdf(d, y, th), distrib_cdf(r, y, th), tolerance = 1e-12)
    expect_equal(distrib_quantile(d, c(0.1, 0.5, 0.9), th),
                 distrib_quantile(r, c(0.1, 0.5, 0.9), th), tolerance = 1e-10)
    expect_same_list(distrib_gradient(d, y, th), distrib_gradient(r, y, th), 1e-9, "gradient")
    expect_same_list(distrib_hessian(d, y, th), distrib_hessian(r, y, th), 1e-9, "hessian")
    expect_same_list(distrib_deriv3(d, y, th), distrib_deriv3(r, y, th), 1e-8, "deriv3")
    expect_same_list(distrib_deriv4(d, y, th), distrib_deriv4(r, y, th), 1e-7, "deriv4")
    expect_same_list(distrib_expected_hessian(d, y, th),
                     distrib_expected_hessian(r, y, th), 1e-9, "expected")
    expect_same_list(distrib_deriv3(d, y, th, expected = TRUE),
                     distrib_deriv3(r, y, th, expected = TRUE), 1e-8, "deriv3 E")
    expect_same_list(distrib_deriv4(d, y, th, expected = TRUE),
                     distrib_deriv4(r, y, th, expected = TRUE), 1e-7, "deriv4 E")
    expect_same_list(distrib_dexpected_hessian(d, y, th),
                     distrib_dexpected_hessian(r, y, th), 1e-6, "dexpected")
    expect_same_list(distrib_d2expected_hessian(d, y, th),
                     distrib_d2expected_hessian(r, y, th), 1e-5, "d2expected")
    for (f in c("distrib_grad_y", "distrib_hess_y", "distrib_deriv3_y",
                "distrib_deriv4_y")) {
      expect_equal(get(f)(d, y, th), get(f)(r, y, th), tolerance = 1e-8, label = f)
    }
    for (f in c("distrib_cross_y", "distrib_cross2_y", "distrib_grad_y_hess",
                "distrib_hess_y_hess")) {
      # the reference's hess_y_hess carries a stencil error of about 5e-7 in
      # sigma; the closed form is checked against numDeriv below
      expect_same_list(get(f)(d, y, th), get(f)(r, y, th),
                       if (f == "distrib_hess_y_hess") 1e-5 else 1e-8, f)
    }
    for (lt in c(TRUE, FALSE)) for (lg in c(TRUE, FALSE)) {
      expect_same_list(distrib_grad_cdf(d, y, th, lower.tail = lt, log = lg),
                       distrib_grad_cdf(r, y, th, lower.tail = lt, log = lg), 1e-8, "grad_cdf")
      expect_same_list(distrib_hess_cdf(d, y, th, lower.tail = lt, log = lg),
                       distrib_hess_cdf(r, y, th, lower.tail = lt, log = lg), 1e-8, "hess_cdf")
      expect_same_list(distrib_deriv3_cdf(d, y, th, lower.tail = lt, log = lg),
                       distrib_deriv3_cdf(r, y, th, lower.tail = lt, log = lg), 1e-7, "d3_cdf")
      expect_same_list(distrib_deriv4_cdf(d, y, th, lower.tail = lt, log = lg),
                       distrib_deriv4_cdf(r, y, th, lower.tail = lt, log = lg), 1e-6, "d4_cdf")
    }
  }
})

test_that("weibull3's fifth order is the derivative of its fourth", {
  d <- weibull3_distrib()
  nm <- c("mean", "sigma")
  for (th in w3_cases) {
    y <- w3_y * th$mean
    p0 <- c(th$mean, th$sigma)
    hi <- distrib_deriv5(d, y, th)
    expect_length(hi, 6)
    for (key in names(hi)) {
      parts <- strsplit(key, "_")[[1]]
      last <- match(parts[5], nm)
      lo_key <- paste(parts[1:4], collapse = "_")
      num <- sapply(seq_along(y), function(i) numDeriv::grad(function(p)
        distrib_deriv4(d, y[i], list(mean = p[1], sigma = p[2]))[[lo_key]], p0)[last])
      # the reference differences an analytic fourth order with numDeriv,
      # whose own noise was measured at about 6e-6 relative on these cases
      expect_equal(hi[[key]], num, tolerance = 1e-4, label = key)
    }
  }
})

test_that("weibull3's moments are its mean and the Weibull's shape", {
  d <- weibull3_distrib()
  th <- list(mean = 4, sigma = 1.7)
  b <- 4 / gamma(1 + 1 / 1.7)
  m <- function(k) integrate(function(t) (t - 4)^k * distrib_pdf(d, t, th), 0, Inf,
                             rel.tol = 1e-10)$value
  expect_equal(mean(d, th), integrate(function(t) t * distrib_pdf(d, t, th), 0, Inf,
                                      rel.tol = 1e-10)$value, tolerance = 1e-8)
  expect_equal(variance(d, th), m(2), tolerance = 1e-8)
  expect_equal(skewness(d, th), m(3) / m(2)^1.5, tolerance = 1e-7)
  expect_equal(kurtosis(d, th), m(4) / m(2)^2 - 3, tolerance = 1e-6)
  expect_equal(distrib_pdf(d, c(1, 3), th), dweibull(c(1, 3), shape = 1.7, scale = b))
})

test_that("weibull3's second mixed derivative is the derivative of hess_y", {
  d <- weibull3_distrib()
  th <- list(mean = 4, sigma = 1.7)
  y <- c(0.3, 1.1, 2.7, 8) * 4
  got <- distrib_hess_y_hess(d, y, th)
  for (i in seq_along(y)) {
    h <- numDeriv::hessian(function(p) distrib_hess_y(d, y[i],
                                                      list(mean = p[1], sigma = p[2])),
                           c(4, 1.7))
    expect_equal(c(got$mean_mean[i], got$sigma_sigma[i], got$mean_sigma[i]),
                 c(h[1, 1], h[2, 2], h[1, 2]), tolerance = 1e-8)
  }
})
