# The information objects, and the choice of the matrix behind the standard
# errors of a fit after the fit.

test_that("the constructors return objects of the classes of the hierarchy", {
  expect_true(S7::S7_inherits(observed_info(), Information))
  expect_true(S7::S7_inherits(observed_info(), ObservedInfo))
  expect_true(S7::S7_inherits(expected_info(), Information))
  expect_true(S7::S7_inherits(expected_info(), ExpectedInfo))
  expect_false(S7::S7_inherits(observed_info(), ExpectedInfo))
  expect_error(Information())
})

test_that("expected_info() carries and validates its approximation", {
  expect_true(S7::S7_inherits(expected_info()@approx, ExpectOPG))
  a <- expected_info(approx = expect_mc(nsim = 300))
  expect_true(S7::S7_inherits(a@approx, ExpectMC))
  expect_identical(a@approx@nsim, 300)
  for (bad in list("opg", 1, NULL, observed_info())) {
    expect_error(expected_info(approx = bad), "must be an ExpectApprox object")
  }
  expect_false("approx" %in% names(S7::S7_class(observed_info())@properties))
})

test_that("the print methods name the matrix", {
  expect_output(print(observed_info()), "^Observed information$")
  expect_output(print(expected_info()), "Expected information \\(expectation by the outer product")
  expect_output(print(expected_info(approx = expect_mc(nsim = 500))),
                "Monte Carlo \\(nsim = 500\\)")
})

test_that("vcov() rejects an argument that is not an information object", {
  set.seed(1)
  d <- gaussian1_distrib()
  fit <- fit_distrib(d, distrib_rng(d, 100, list(mu = 1, sigma = 2)))
  msg <- "'information' must be observed_info\\(\\) or expected_info\\(\\)"
  expect_error(vcov(fit, information = "observed"), msg)
  expect_error(vcov(fit, information = 1), msg)
  expect_error(vcov(fit, information = expect_opg()), msg)
  expect_error(confint(fit, information = "expected"), msg)
  expect_error(print(fit, information = "expected"), msg)
})

test_that("the default is the observed information for every fit", {
  set.seed(2)
  d <- gamma1_distrib()
  y <- distrib_rng(d, 200, list(mu = 3, phi = 0.6))
  for (m in list(fisher_scoring(), optimizers7::newton())) {
    fit <- fit_distrib(d, y, method = m)
    H <- -distributions7:::fit_hess_matrix(d, y, as.list(fit@coefficients),
                                           expected = FALSE)
    expect_equal(unname(vcov(fit, scale = "link")), unname(solve(H)),
                 tolerance = 1e-12)
    expect_identical(vcov(fit), vcov(fit, information = observed_info()))
    expect_equal(fit@se, sqrt(diag(vcov(fit))))
    expect_equal(fit@se_eta, sqrt(diag(vcov(fit, scale = "link"))))
  }
})

test_that("a density-only family reads its expected information through approx", {
  Dens <- S7::new_class("DensInformation", parent = continuous_distrib)
  S7::method(distrib_pdf, Dens) <- function(distrib, y, theta, log = FALSE, ...) {
    stats::dnorm(y, theta[[1]], theta[[2]], log = log)
  }
  d <- Dens(distrib_name = "dens", dimension = "univariate",
            bounds = c(-Inf, Inf), params = c("mu", "sigma"),
            params_interpretation = c(mu = "location", sigma = "scale"),
            n_params = 2L,
            params_bounds = list(mu = c(-Inf, Inf), sigma = c(0, Inf)),
            link_params = list(mu = linkfunctions7::identity_link(),
                               sigma = linkfunctions7::log_link()),
            params_smooth = c(mu = TRUE, sigma = TRUE))
  set.seed(3)
  y <- stats::rnorm(150, 1, 2)
  fit <- fit_distrib(d, y)
  V_opg <- vcov(fit, information = expected_info())
  V_bar <- vcov(fit, information = expected_info(approx = expect_bartlett()))
  V_obs <- vcov(fit)
  expect_false(isTRUE(all.equal(V_opg, V_obs)))
  expect_equal(V_bar, V_obs, tolerance = 1e-3)
  expect_equal(V_obs, vcov(fit, information = observed_info()))
})

test_that("print, confint and mv_summary follow the chosen information", {
  set.seed(4)
  d <- negbin2_distrib()
  y <- distrib_rng(d, 300, list(mu = 5, theta = 2))
  fit <- fit_distrib(d, y)
  se_o <- fit@se
  se_e <- sqrt(diag(vcov(fit, information = expected_info())))
  expect_false(isTRUE(all.equal(unname(se_o), unname(se_e))))

  ci_o <- confint(fit, scale = "link")
  ci_e <- confint(fit, scale = "link", information = expected_info())
  z <- stats::qnorm(0.975)
  se_eta_e <- sqrt(diag(vcov(fit, scale = "link", information = expected_info())))
  expect_equal(unname(ci_e[, 2] - ci_e[, 1]), unname(2 * z * se_eta_e))
  expect_equal(unname(ci_o), unname(confint(fit, scale = "link",
                                            information = observed_info())))

  txt_o <- capture.output(print(fit))
  txt_e <- capture.output(print(fit, information = expected_info()))
  expect_false(any(grepl("expected information", txt_o)))
  expect_true(any(grepl("Standard errors from the expected information", txt_e)))
  expect_false(identical(txt_o, txt_e))

  set.seed(5)
  m <- mvgaussian1_distrib(2)
  th <- list(mu1 = 0, mu2 = 1, sigma_log_L1 = 0, sigma_log_L2 = 0, sigma_L2.1 = 0.7)
  fm <- fit_distrib(m, distrib_rng(m, 400, th))
  s_o <- mv_summary(fm)
  s_e <- mv_summary(fm, information = expected_info())
  expect_identical(s_o[["Estimate"]], s_e[["Estimate"]])
  expect_equal(s_o[["Std. Error"]], s_e[["Std. Error"]], tolerance = 1e-3)
})

test_that("a family with a kink has no observed standard error and an expected one", {
  set.seed(13)
  d <- laplace_distrib()
  y <- distrib_rng(d, 500, list(mu = 1, sigma = 2))
  fit <- fit_distrib(d, y, start = list(mu = 0, sigma = 1))
  expect_true(all(is.na(fit@se)))
  V <- vcov(fit, information = expected_info())
  expect_true(all(is.finite(V)))
})
