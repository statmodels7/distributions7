# A parameter name containing an underscore.
#
# Component names are built by joining parameter names with "_", so recovering a
# multi-index by splitting one back apart is ambiguous the moment a parameter
# name contains an underscore of its own: "mu_log_scale_log_scale" splits into
# five pieces, none of which except the first matches a parameter. Nothing in the
# shipped catalog is affected -- none of the fourteen has such a name -- so the
# only way this surfaces is through a distribution someone defines themselves,
# which is precisely the case the package is built to support.

# A distribution defining nothing but its density, so that every derivative comes
# from the numerical fallbacks. The class and its method are created once:
# S7::method(gen, cls) <- fn mutates the generic in place, so re-registering it
# inside a helper called per test prints "Overwriting method" for each call.
PlainGauss <- S7::new_class("PlainGauss", parent = continuous_distrib)
S7::method(distrib_pdf, PlainGauss) <- function(distrib, y, theta, log = FALSE) {
  th <- align_theta(distrib, theta)
  stats::dnorm(y, th[[1]], th[[2]], log = log)
}

# The two instances differ only in what the scale parameter is called.
underscore_gauss <- function(params = c("mu", "log_scale")) {
  PlainGauss(
    distrib_name = "plain gauss", dimension = "univariate",
    bounds = c(-Inf, Inf), params = params, n_params = 2L,
    params_interpretation = stats::setNames(c("location", "scale"), params),
    params_bounds = stats::setNames(list(c(-Inf, Inf), c(0, Inf)), params),
    link_params = stats::setNames(
      list(linkfunctions7::identity_link(), linkfunctions7::log_link()), params)
  )
}

test_that("higher-order numerical derivatives accept an underscore in a name", {
  y <- c(-0.5, 0.2, 1.1)
  d_us <- underscore_gauss(c("mu", "log_scale"))
  d_ok <- underscore_gauss(c("mu", "sigma"))
  th_us <- list(mu = 0, log_scale = 1)
  th_ok <- list(mu = 0, sigma = 1)

  for (fn in list(numerical_deriv3, numerical_deriv4)) {
    a <- fn(d_us, y, th_us)
    b <- fn(d_ok, y, th_ok)
    # same numbers, only the keys differ
    expect_equal(unname(a), unname(b))
    expect_equal(length(a), length(b))
  }

  expect_equal(names(numerical_deriv3(d_us, y, th_us)),
               deriv_names(c("mu", "log_scale"), 3))
})

test_that("the Bartlett expected Hessian accepts an underscore in a name", {
  y <- c(-0.5, 0.2, 1.1)
  d_us <- underscore_gauss(c("mu", "log_scale"))
  d_ok <- underscore_gauss(c("mu", "sigma"))

  a <- distrib_expected_hessian(d_us, y, list(mu = 0, log_scale = 1),
                                approx = "bartlett")
  b <- distrib_expected_hessian(d_ok, y, list(mu = 0, sigma = 1),
                                approx = "bartlett")
  expect_equal(unname(a), unname(b))

  # and it is the right answer: for a standard Gaussian the information is
  # 1/sigma^2 in mu and 2/sigma^2 in sigma, so E[H] is -1 and -2
  expect_equal(a[["mu_mu"]], rep(-1, 3), tolerance = 1e-6)
  expect_equal(a[["log_scale_log_scale"]], rep(-2, 3), tolerance = 1e-6)
})

test_that("deriv_indices and deriv_names describe the same enumeration", {
  for (params in list(c("mu", "sigma"), c("mu", "log_scale"), c("a", "b", "c"))) {
    for (order in 1:4) {
      idx <- deriv_indices(params, order)
      nms <- deriv_names(params, order)
      expect_equal(length(idx), length(nms))
      expect_equal(vapply(idx, function(i) paste(params[i], collapse = "_"), ""),
                   nms)
      # non-decreasing, so each multi-index is listed once
      expect_true(all(vapply(idx, function(i) all(diff(i) >= 0), TRUE)))
      expect_true(all(vapply(idx, length, 0L) == order))
    }
  }
})

test_that("hess_pairs inverts hess_names for awkward parameter names", {
  params <- c("mu", "log_scale")
  pr <- hess_pairs(params)
  expect_equal(names(pr), hess_names(params))
  expect_equal(pr[["log_scale_log_scale"]], c("log_scale", "log_scale"))
  expect_equal(pr[["mu_log_scale"]], c("mu", "log_scale"))
})

test_that("order_indices agrees with deriv_names ordering", {
  for (params in list(c("mu", "sigma"), c("mu", "log_scale"), c("a", "b", "c"))) {
    for (order in 3:4) {
      expect_equal(
        vapply(order_indices(params, order), paste, "", collapse = "_"),
        deriv_names(params, order)
      )
    }
  }
})
