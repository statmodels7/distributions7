# student_t2 written out, against the same law obtained through
# reparametrize() on student_t1 (student_t2_by_reparam(), an independent
# implementation), against numDeriv where the reference stops, and against
# exact values at large nu, where both of those lose their digits.

t2_cases <- list(list(mu = 0.2, sigma = 1.3, nu = 8),
                 list(mu = -1, sigma = 0.4, nu = 3.5),
                 list(mu = 5, sigma = 4, nu = 40))
t2_z <- c(-3.1, -0.4, 0.7, 2.2)

expect_same_list <- function(a, b, tol, label) {
  expect_setequal(names(a), names(b))
  for (nm in names(b)) {
    expect_equal(rep_len(a[[nm]], length(b[[nm]])), b[[nm]], tolerance = tol,
                 label = paste(label, nm))
  }
}

test_that("student_t2 is the reparametrized student_t1 on every surface", {
  d <- student_t2_distrib()
  r <- student_t2_by_reparam()
  for (th in t2_cases) {
    y <- th$mu + t2_z * th$sigma
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
    # the reference's expected third and fourth orders are quadratures
    expect_same_list(distrib_deriv3(d, y, th, expected = TRUE),
                     distrib_deriv3(r, y, th, expected = TRUE), 1e-6, "deriv3 E")
    expect_same_list(distrib_deriv4(d, y, th, expected = TRUE),
                     distrib_deriv4(r, y, th, expected = TRUE), 1e-6, "deriv4 E")
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
      expect_same_list(get(f)(d, y, th), get(f)(r, y, th), 1e-6, f)
    }
    for (lt in c(TRUE, FALSE)) for (lg in c(TRUE, FALSE)) {
      expect_same_list(distrib_grad_cdf(d, y, th, lower.tail = lt, log = lg),
                       distrib_grad_cdf(r, y, th, lower.tail = lt, log = lg), 1e-6, "grad_cdf")
      expect_same_list(distrib_hess_cdf(d, y, th, lower.tail = lt, log = lg),
                       distrib_hess_cdf(r, y, th, lower.tail = lt, log = lg), 1e-5, "hess_cdf")
    }
  }
})

test_that("student_t2's expectations are integrals of its derivatives", {
  d <- student_t2_distrib()
  th <- list(mu = 0.2, sigma = 1.3, nu = 6.5)
  E <- function(fun, key) {
    integrate(function(t) fun(t)[[key]] * distrib_pdf(d, t, th), -Inf, Inf,
              rel.tol = 1e-11)$value
  }
  e3 <- distrib_deriv3(d, 0, th, expected = TRUE)
  e4 <- distrib_deriv4(d, 0, th, expected = TRUE)
  for (key in names(e3)) {
    expect_equal(e3[[key]], E(function(t) distrib_deriv3(d, t, th), key),
                 tolerance = 1e-8, label = key)
  }
  for (key in names(e4)) {
    expect_equal(e4[[key]], E(function(t) distrib_deriv4(d, t, th), key),
                 tolerance = 1e-8, label = key)
  }
})

test_that("student_t2's fifth order is the derivative of its fourth", {
  d <- student_t2_distrib()
  nm <- c("mu", "sigma", "nu")
  for (th in t2_cases[1:2]) {
    y <- th$mu + t2_z * th$sigma
    p0 <- c(th$mu, th$sigma, th$nu)
    hi <- distrib_deriv5(d, y, th)
    expect_length(hi, 21)
    for (key in names(hi)) {
      parts <- strsplit(key, "_")[[1]]
      last <- match(parts[5], nm)
      lo_key <- paste(parts[1:4], collapse = "_")
      num <- sapply(seq_along(y), function(i) numDeriv::grad(function(p)
        distrib_deriv4(d, y[i], list(mu = p[1], sigma = p[2], nu = p[3]))[[lo_key]],
        p0)[last])
      expect_equal(hi[[key]], num, tolerance = 1e-5, label = key)
    }
  }
})

test_that("student_t2's second mixed derivative is the derivative of hess_y", {
  d <- student_t2_distrib()
  th <- list(mu = 0.2, sigma = 1.3, nu = 8)
  y <- th$mu + t2_z * th$sigma
  got <- distrib_hess_y_hess(d, y, th)
  for (i in seq_along(y)) {
    h <- numDeriv::hessian(function(p) distrib_hess_y(d, y[i],
                                                      list(mu = p[1], sigma = p[2], nu = p[3])),
                           c(0.2, 1.3, 8))
    expect_equal(c(got$mu_mu[i], got$sigma_sigma[i], got$nu_nu[i], got$mu_sigma[i],
                   got$mu_nu[i], got$sigma_nu[i]),
                 c(h[1, 1], h[2, 2], h[3, 3], h[1, 2], h[1, 3], h[2, 3]), tolerance = 1e-7)
  }
})

test_that("student_t2's derivatives in nu keep their digits as nu grows", {
  # exact values from the closed forms at 250 digits (mpmath, written by
  # stabilita/gen_student_t2.py --reference), at y = 0.7, mu = 0.2, sigma = 1.3.
  # Every one vanishes as a power of 1/nu and is a difference of terms of a
  # lower order, so the direct forms lose every digit by nu = 1e8.
  d <- student_t2_distrib()
  ref <- list(
    list(1e6, "gradient", "nu", -5.3357842175450464011e-13),
    list(1e10, "gradient", "nu", -5.3357725581493539778e-21),
    list(1e6, "hessian", "nu_nu", 1.0671580095668976958e-18),
    list(1e10, "hessian", "nu_nu", 1.0671545117464762417e-30),
    list(1e10, "deriv3", "mu_nu_nu", 1.6876159807703140239e-30),
    list(1e6, "deriv3", "nu_nu_nu", -3.2014775268777857923e-24),
    list(1e10, "deriv4", "sigma_sigma_nu_nu", -1.4460899537695253821e-30),
    list(1e6, "deriv4", "nu_nu_nu_nu", 1.2805924100233223831e-29),
    list(1e10, "deriv5", "nu_nu_nu_nu_nu", -6.4029270725777554809e-59),
    list(1e6, "expected_hessian", "nu_nu", -1.499997000021499979e-24),
    list(1e10, "deriv4_expected", "nu_nu_nu_nu", -5.3999999967600000041e-59),
    list(1e6, "dexpected2", "nu_nu_nu_nu", -2.9999910000902998824e-35))
  fun <- list(
    gradient = function(th) distrib_gradient(d, 0.7, th),
    hessian = function(th) distrib_hessian(d, 0.7, th),
    deriv3 = function(th) distrib_deriv3(d, 0.7, th),
    deriv4 = function(th) distrib_deriv4(d, 0.7, th),
    deriv5 = function(th) distrib_deriv5(d, 0.7, th),
    expected_hessian = function(th) distrib_expected_hessian(d, 0.7, th),
    deriv4_expected = function(th) distrib_deriv4(d, 0.7, th, expected = TRUE),
    dexpected2 = function(th) distrib_d2expected_hessian(d, 0.7, th))
  for (r in ref) {
    th <- list(mu = 0.2, sigma = 1.3, nu = r[[1]])
    expect_equal(fun[[r[[2]]]](th)[[r[[3]]]], r[[4]], tolerance = 1e-12,
                 label = paste(r[[2]], r[[3]], r[[1]]))
  }
})

test_that("student_t2's surfaces stay finite where the link can take nu", {
  d <- student_t2_distrib()
  y <- c(-40, 0.7, 3)
  for (v in c(1e150, 1e300, .Machine$double.xmax)) {
    th <- list(mu = 0.2, sigma = 1.3, nu = v)
    outs <- c(distrib_gradient(d, y, th), distrib_hessian(d, y, th),
              distrib_deriv3(d, y, th), distrib_deriv4(d, y, th),
              distrib_deriv5(d, y, th), distrib_expected_hessian(d, y, th),
              distrib_deriv3(d, y, th, expected = TRUE),
              distrib_deriv4(d, y, th, expected = TRUE),
              distrib_dexpected_hessian(d, y, th),
              distrib_d2expected_hessian(d, y, th),
              distrib_cross_y(d, y, th), distrib_hess_y_hess(d, y, th))
    for (nm in names(outs)) {
      expect_true(all(is.finite(outs[[nm]])), label = paste(nm, v))
    }
    # the location and scale components reach the gaussian's
    expect_equal(distrib_gradient(d, y, th)$mu, (y - 0.2) / 1.3^2)
    expect_equal(distrib_expected_hessian(d, y, th)$sigma_sigma, rep(-2 / 1.3^2, 3))
  }
})

test_that("student_t2's moments are its location and standard deviation", {
  d <- student_t2_distrib()
  th <- list(mu = 1, sigma = 2, nu = 6)
  m <- function(k) integrate(function(t) (t - 1)^k * distrib_pdf(d, t, th), -Inf, Inf,
                             rel.tol = 1e-10)$value
  expect_equal(mean(d, th), 1)
  expect_equal(variance(d, th), m(2), tolerance = 1e-8)
  expect_equal(variance(d, th), 4)
  expect_equal(skewness(d, th), 0)
  expect_equal(kurtosis(d, th), 3)
  expect_equal(kurtosis(d, list(mu = 0, sigma = 1, nu = c(3, 10))), c(Inf, 1))
  expect_equal(skewness(d, list(mu = 0, sigma = 1, nu = c(2.5, 10))), c(NaN, 0))
})

test_that("student_t2 passes its own validator to order five", {
  d <- student_t2_distrib()
  set.seed(3)
  res <- check_distrib(d, list(mu = 0.2, sigma = 1.3, nu = 7), orders = 1:5,
                       verbose = FALSE)
  expect_true(all(res$status == "OK"),
              info = paste(res$check[res$status != "OK"], collapse = ", "))
})
