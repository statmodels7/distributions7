# gengamma2 written out, against the same law obtained through
# reparametrize() on gengamma1 (gengamma2_by_reparam(), an independent
# implementation), against numDeriv where the reference stops, and against
# exact values where d/p is large and the reference cannot be evaluated.

gg2_cases <- list(list(mean = 1.2, d = 2, p = 1.5), list(mean = 0.4, d = 0.6, p = 0.8),
                  list(mean = 5, d = 3, p = 4))
gg2_y <- c(0.2, 0.7, 1.3, 2.9)

expect_same_list <- function(a, b, tol, label) {
  expect_setequal(names(a), names(b))
  for (nm in names(b)) {
    expect_equal(rep_len(a[[nm]], length(b[[nm]])), b[[nm]], tolerance = tol,
                 label = paste(label, nm))
  }
}

test_that("gengamma2 is the reparametrized gengamma1 on every surface", {
  d <- gengamma2_distrib()
  r <- gengamma2_by_reparam()
  for (th in gg2_cases) {
    y <- gg2_y * th$mean
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
                     distrib_deriv3(r, y, th, expected = TRUE), 1e-6, "deriv3 E")
    expect_same_list(distrib_deriv4(d, y, th, expected = TRUE),
                     distrib_deriv4(r, y, th, expected = TRUE), 1e-5, "deriv4 E")
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
      # the reference's hess_y_hess carries a stencil error of up to 1.3e-4;
      # the closed form is checked against numDeriv below
      expect_same_list(get(f)(d, y, th), get(f)(r, y, th),
                       if (f == "distrib_hess_y_hess") 1e-3 else 1e-6, f)
    }
  }
})

test_that("gengamma2's expectations are integrals of its derivatives", {
  d <- gengamma2_distrib()
  th <- list(mean = 1.2, d = 2, p = 1.5)
  E <- function(fun, key) {
    integrate(function(t) fun(t)[[key]] * distrib_pdf(d, t, th), 0, Inf,
              rel.tol = 1e-11)$value
  }
  e2 <- distrib_expected_hessian(d, 1, th)
  e3 <- distrib_deriv3(d, 1, th, expected = TRUE)
  e4 <- distrib_deriv4(d, 1, th, expected = TRUE)
  for (key in names(e2)) {
    expect_equal(e2[[key]], E(function(t) distrib_hessian(d, t, th), key),
                 tolerance = 1e-8, label = key)
  }
  for (key in names(e3)) {
    expect_equal(e3[[key]], E(function(t) distrib_deriv3(d, t, th), key),
                 tolerance = 1e-8, label = key)
  }
  for (key in names(e4)) {
    expect_equal(e4[[key]], E(function(t) distrib_deriv4(d, t, th), key),
                 tolerance = 1e-8, label = key)
  }
})

test_that("gengamma2's fifth order is the derivative of its fourth", {
  d <- gengamma2_distrib()
  nm <- c("mean", "d", "p")
  for (th in gg2_cases[1:2]) {
    y <- gg2_y * th$mean
    p0 <- c(th$mean, th$d, th$p)
    hi <- distrib_deriv5(d, y, th)
    expect_length(hi, 21)
    for (key in names(hi)) {
      parts <- strsplit(key, "_")[[1]]
      last <- match(parts[5], nm)
      lo_key <- paste(parts[1:4], collapse = "_")
      num <- sapply(seq_along(y), function(i) numDeriv::grad(function(p)
        distrib_deriv4(d, y[i], list(mean = p[1], d = p[2], p = p[3]))[[lo_key]],
        p0)[last])
      expect_equal(hi[[key]], num, tolerance = 1e-5, label = key)
    }
  }
})

test_that("gengamma2's second mixed derivative is the derivative of hess_y", {
  d <- gengamma2_distrib()
  th <- list(mean = 1.2, d = 2, p = 1.5)
  y <- gg2_y * 1.2
  got <- distrib_hess_y_hess(d, y, th)
  for (i in seq_along(y)) {
    h <- numDeriv::hessian(function(p) distrib_hess_y(d, y[i],
                                                      list(mean = p[1], d = p[2], p = p[3])),
                           c(1.2, 2, 1.5))
    expect_equal(c(got$mean_mean[i], got$d_d[i], got$p_p[i], got$mean_d[i],
                   got$mean_p[i], got$d_p[i]),
                 c(h[1, 1], h[2, 2], h[3, 3], h[1, 2], h[1, 3], h[2, 3]), tolerance = 1e-7)
  }
})

test_that("gengamma2 agrees with exact values, including where d/p is large", {
  # exact values at 1200 bits (mpmath, written by stabilita/gen_gengamma2.py
  # --reference). Where d/p and d are large every derivative in d and p is a
  # difference of terms agreeing to several orders; the kernels carry the
  # cancellation on coefficients formed in double-double.
  d <- gengamma2_distrib()
  ref <- list(
    list(c(0.7, 1.2, 2, 1.5), "gradient", "p", -0.048648209007135585321),
    list(c(5, 4, 40, 0.5), "gradient", "p", 0.045166658592804386287),
    list(c(0.05, 0.4, 0.6, 0.8), "hessian", "d_p", -0.7954831938368368325),
    list(c(5, 4, 40, 0.5), "hessian", "d_p", -0.026999931702950628686),
    list(c(5, 4, 40, 0.5), "expected_hessian", "p_p", -1.995630499579616036),
    list(c(0.7, 1.2, 2, 1.5), "deriv3", "d_p_p", 0.032117938478213033688),
    list(c(5, 4, 40, 0.5), "deriv3", "d_p_p", -0.0027934719548109612014),
    list(c(0.05, 0.4, 0.6, 0.8), "deriv4", "d_d_p_p", -0.36235521790676203018),
    list(c(5, 4, 40, 0.5), "deriv4", "d_d_p_p", -0.000061003170043260093241),
    list(c(5, 4, 40, 0.5), "deriv4_expected", "d_d_p_p", -0.000060344634339442042899),
    list(c(0.7, 1.2, 2, 1.5), "deriv5", "mean_d_p_p_p", -0.049496580113041424474),
    list(c(5, 4, 40, 0.5), "deriv5", "mean_d_p_p_p", 0.0030936970011387884367),
    list(c(5, 4, 40, 0.5), "dexpected2", "d_p_d_p", -0.0013102154029237909718),
    # towards the lognormal, where the former forms lost 1e-3 to all digits
    list(c(4.04, 4, 1e4, 1), "gradient", "p", 0.0066499663357041036988),
    list(c(4.04, 4, 1e4, 1), "hessian", "d_p", -0.000049833320033099623408),
    list(c(4.01, 4, 1e5, 0.5), "hessian", "p_p", -2.000536255467123248),
    list(c(4.04, 4, 1e4, 1), "expected_hessian", "p_p", -0.49996666833329999667),
    list(c(4.04, 4, 1e4, 1), "deriv3", "d_p_p", -3.2834327746922792888e-7),
    list(c(4.01, 4, 1e5, 0.5), "deriv3", "d_p_p", -4.9937008792094294129e-9),
    list(c(4.04, 4, 1e4, 1), "deriv4", "d_d_p_p", -5.0163377681562863258e-13),
    list(c(4.01, 4, 1e5, 0.5), "deriv4", "d_p_p_p", -1.2097146005926567587e-9),
    list(c(4.04, 4, 1e4, 1), "deriv4_expected", "d_d_p_p", -4.9996666522228893429e-13),
    list(c(4.04, 4, 1e4, 1), "deriv5", "mean_d_p_p_p", 2.4874981087059559789e-7),
    list(c(4.04, 4, 1e4, 1), "deriv5", "d_p_p_p_p", 2.9977331310685413578e-8),
    list(c(4.01, 4, 1e5, 0.5), "deriv5", "d_d_p_p_p", 2.3999638424403137175e-14))
  fun <- list(
    gradient = function(y, th) distrib_gradient(d, y, th),
    hessian = function(y, th) distrib_hessian(d, y, th),
    deriv3 = function(y, th) distrib_deriv3(d, y, th),
    deriv4 = function(y, th) distrib_deriv4(d, y, th),
    deriv5 = function(y, th) distrib_deriv5(d, y, th),
    expected_hessian = function(y, th) distrib_expected_hessian(d, y, th),
    deriv4_expected = function(y, th) distrib_deriv4(d, y, th, expected = TRUE),
    dexpected2 = function(y, th) distrib_d2expected_hessian(d, y, th))
  for (r in ref) {
    pt <- r[[1]]
    th <- list(mean = pt[2], d = pt[3], p = pt[4])
    expect_equal(fun[[r[[2]]]](pt[1], th)[[r[[3]]]], r[[4]], tolerance = 1e-11,
                 label = paste(r[[2]], r[[3]], pt[3] / pt[4]))
  }
})

test_that("gengamma2's moments are its mean and the shape's ratios", {
  d <- gengamma2_distrib()
  th <- list(mean = 5, d = 3, p = 1.5)
  m <- function(k) integrate(function(t) (t - 5)^k * distrib_pdf(d, t, th), 0, Inf,
                             rel.tol = 1e-10)$value
  expect_equal(mean(d, th), integrate(function(t) t * distrib_pdf(d, t, th), 0, Inf,
                                      rel.tol = 1e-10)$value, tolerance = 1e-8)
  expect_equal(variance(d, th), m(2), tolerance = 1e-8)
  expect_equal(skewness(d, th), m(3) / m(2)^1.5, tolerance = 1e-7)
  expect_equal(kurtosis(d, th), m(4) / m(2)^2 - 3, tolerance = 1e-6)
  r <- gengamma2_by_reparam()
  expect_equal(variance(d, th), variance(r, th), tolerance = 1e-12)
})

test_that("gengamma2 passes its own validator to order five", {
  d <- gengamma2_distrib()
  set.seed(3)
  res <- check_distrib(d, list(mean = 1.2, d = 2, p = 1.5), orders = 1:5,
                       verbose = FALSE)
  expect_true(all(res$status == "OK"),
              info = paste(res$check[res$status != "OK"], collapse = ", "))
})
