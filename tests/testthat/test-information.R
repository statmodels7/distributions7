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

test_that("expected_info() carries and validates its method", {
  expect_null(expected_info()@method)
  expect_null(expected_info(method = NULL)@method)
  a <- expected_info(method = expect_mc(nsim = 300))
  expect_true(S7::S7_inherits(a@method, ExpectMC))
  expect_identical(a@method@nsim, 300)
  expect_true(S7::S7_inherits(expected_info(method = expect_opg())@method, ExpectOPG))
  for (bad in list("opg", 1, observed_info())) {
    expect_error(expected_info(method = bad), "must be NULL or an ExpectMethod object")
  }
  expect_false("method" %in% names(S7::S7_class(observed_info())@properties))
})

test_that("the print methods name the matrix", {
  expect_output(print(observed_info()), "^Observed information$")
  expect_output(print(expected_info()),
                "Expected information \\(closed form where the family has one")
  expect_output(print(expected_info(method = expect_opg())),
                "Expected information \\(expectation by the outer product")
  expect_output(print(expected_info(method = expect_mc(nsim = 500))),
                "Monte Carlo \\(nsim = 500\\)")
})

test_that("vcov() rejects an argument that is not an information object", {
  set.seed(1)
  d <- gaussian1_distrib()
  fit <- fit_distrib(d, distrib_rng(d, 100, list(mu = 1, sigma = 2)))
  msg <- "'information' must be NULL, observed_info\\(\\) or expected_info\\(\\)"
  expect_error(vcov(fit, information = "observed"), msg)
  expect_error(vcov(fit, information = 1), msg)
  expect_error(vcov(fit, information = expect_opg()), msg)
  expect_error(confint(fit, information = "expected"), msg)
  expect_error(print(fit, information = "expected"), msg)
})

test_that("the default is the observed information for a smooth family", {
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
    expect_identical(vcov(fit), vcov(fit, information = NULL))
    expect_true(S7::S7_inherits(distributions7::default_information(d), ObservedInfo))
    expect_equal(fit@se, sqrt(diag(vcov(fit))))
    expect_equal(fit@se_eta, sqrt(diag(vcov(fit, scale = "link"))))
  }
})

test_that("a density-only family reads its expected information through method", {
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
  V_bar <- vcov(fit, information = expected_info(method = expect_integrate(identity = "score")))
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

test_that("a family with a kink takes the expected information by default", {
  set.seed(13)
  d <- laplace_distrib()
  y <- distrib_rng(d, 500, list(mu = 1, sigma = 2))
  fit <- fit_distrib(d, y, start = list(mu = 0, sigma = 1))
  expect_true(S7::S7_inherits(distributions7::default_information(d), ExpectedInfo))
  expect_null(distributions7::default_information(d)@method)
  expect_true(all(is.finite(fit@se)))
  V <- vcov(fit)
  expect_true(all(is.finite(V)))
  expect_identical(V, vcov(fit, information = expected_info()))
  expect_equal(fit@se, sqrt(diag(vcov(fit))))
  expect_true(all(is.na(vcov(fit, information = observed_info()))))
  # print states the information that was used, the confidence interval follows
  txt <- capture.output(print(fit))
  expect_true(any(grepl("Standard errors from the expected information", txt)))
  expect_true(all(is.finite(confint(fit))))
  expect_true(all(is.na(confint(fit, information = observed_info()))))
  # the other families with a kink in a parameter
  set.seed(14)
  d2 <- laplace2_distrib()
  y2 <- distrib_rng(d2, 500, list(mu = 1, lambda = 0.5))
  f2 <- fit_distrib(d2, y2, start = list(mu = 0, lambda = 1))
  expect_true(all(is.finite(f2@se)))
  expect_true(S7::S7_inherits(distributions7::default_information(enet_distrib()),
                              ExpectedInfo))
})

test_that("an explicit information reaches print, confint and vcov of a smooth fit", {
  set.seed(15)
  d <- gamma1_distrib()
  fit <- fit_distrib(d, distrib_rng(d, 200, list(mu = 3, phi = 0.6)))
  e <- expected_info()
  expect_equal(unname(sqrt(diag(vcov(fit, information = e)))),
               unname(sqrt(diag(vcov(fit, information = expected_info(method = NULL))))))
  txt <- capture.output(print(fit, information = expected_info(method = expect_opg())))
  expect_true(any(grepl("Standard errors from the expected information", txt)))
  expect_true(any(grepl("outer product of the scores", txt)))
})
