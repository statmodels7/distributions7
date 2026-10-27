# Values returned with method = NULL before `method` was introduced (measured
# on the previous version of the package), fixed here to guard the defaults:
# NULL resolves to the closed form where the family has one, to expect_opg()
# at order 2 and to expect_integrate() (identity = "hessian") above it.

dens_default <- function() {
  Dens <- S7::new_class("DensDefaults", parent = continuous_distrib)
  S7::method(distrib_pdf, Dens) <- function(distrib, y, theta, log = FALSE, ...) {
    stats::dnorm(y, theta[[1]], theta[[2]], log = log)
  }
  Dens(distrib_name = "dens", dimension = "univariate",
       bounds = c(-Inf, Inf), params = c("mu", "sigma"),
       params_interpretation = c(mu = "location", sigma = "scale"),
       n_params = 2L,
       params_bounds = list(mu = c(-Inf, Inf), sigma = c(0, Inf)),
       link_params = list(mu = linkfunctions7::identity_link(),
                          sigma = linkfunctions7::log_link()),
       params_smooth = c(mu = TRUE, sigma = TRUE))
}

test_that("the default of distrib_expected_hessian() is opg for a density-only family", {
  d <- dens_default()
  th <- list(mu = 0.2, sigma = 1.7)
  y <- c(-1, 0.5, 2)
  eh <- distrib_expected_hessian(d, y, th)
  expect_equal(eh$sigma_sigma,
               c(-0.0871048778710322297, -0.3248048716162006966, -0.0050750761251496072),
               tolerance = 1e-12)
  expect_identical(eh, distrib_expected_hessian(d, y, th, method = expect_opg()))
})

test_that("the default of the third and fourth derivatives is integration of the observed derivative", {
  d <- dens_default()
  th <- list(mu = 0.2, sigma = 1.7)
  y <- c(-1, 0.5, 2)
  # the family is a gaussian written as a density only, so the references are
  # the closed forms 2 / sigma^3 and -6 / sigma^4; the default route
  # integrates a stencil derivative, whose result differs by 1.3e-7 (relative)
  # on arm64 macOS and by less on Windows, hence the tolerance
  expect_equal(distrib_deriv3(d, y, th, expected = TRUE)$mu_mu_sigma,
               rep(2 / th$sigma^3, 3), tolerance = 1e-6)
  expect_equal(distrib_deriv4(d, y, th, expected = TRUE)$mu_mu_sigma_sigma,
               rep(-6 / th$sigma^4, 3), tolerance = 1e-6)
  # the pig1 third derivative is not a closed form
  expect_equal(distrib_deriv3(pig1_distrib(), c(0, 1, 3, 5), list(mu = 2, sigma = 0.8),
                              expected = TRUE)$mu_mu_mu,
               rep(0.28580977888728287, 4), tolerance = 1e-7)
})

test_that("the default fit of a density-only family is the Fisher scoring fit on opg", {
  d <- dens_default()
  set.seed(11)
  y <- stats::rnorm(80, 1, 2)
  f <- fit_distrib(d, y)
  expect_equal(unname(c(f@coefficients, f@loglik)),
               c(0.71690470154923591, 1.75380930599392282, -158.45831329624951422),
               tolerance = 1e-7)
  f2 <- fit_distrib(d, y, method = fisher_scoring(method = expect_opg()))
  expect_identical(f@coefficients, f2@coefficients)
})

test_that("the default of the derivatives of the expected information is the closed form", {
  g <- gaussian1_distrib()
  th <- list(mu = 0.2, sigma = 1.7)
  for (gen in list(distrib_dexpected_hessian, distrib_d2expected_hessian,
                   distrib_d3expected_hessian, distrib_d4expected_hessian)) {
    expect_identical(gen(g, 0, th), gen(g, 0, th, method = NULL))
  }
  # and a family without a closed form has no default
  expect_error(distrib_dexpected_hessian(dens_default(), 0, th),
               "approximates its expected information")
})
