# The approximation of an expectation as an S7 object: construction, the
# errors for the strings of the former interface, the routing to the three
# computations, and the defaults of each consumer.

dens_only <- function() {
  Dens <- S7::new_class("DensExpectApprox", parent = continuous_distrib)
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

test_that("the constructors return objects of the classes of the hierarchy", {
  for (a in list(expect_opg(), expect_bartlett(), expect_integrate(), expect_mc())) {
    expect_true(S7::S7_inherits(a, ExpectApprox))
  }
  expect_true(S7::S7_inherits(expect_opg(), ExpectOPG))
  expect_true(S7::S7_inherits(expect_bartlett(), ExpectBartlett))
  expect_true(S7::S7_inherits(expect_integrate(), ExpectIntegrate))
  expect_true(S7::S7_inherits(expect_mc(), ExpectMC))
  expect_false(S7::S7_inherits(expect_opg(), ExpectMC))
  # the parent is abstract
  expect_error(ExpectApprox())
})

test_that("expect_mc() carries nsim and validates it", {
  expect_identical(expect_mc()@nsim, 10000)
  expect_identical(expect_mc(nsim = 250)@nsim, 250)
  expect_identical(expect_mc(nsim = 20L)@nsim, 20L)
  for (bad in list(0, -1, NA_real_, Inf, c(10, 20), "10", NULL, 0.5)) {
    expect_error(expect_mc(nsim = bad), "single positive number")
  }
})

test_that("the print methods name the approximation", {
  expect_output(print(expect_opg()), "outer product")
  expect_output(print(expect_bartlett()), "Bartlett")
  expect_output(print(expect_integrate()), "integration")
  expect_output(print(expect_mc(nsim = 1000)), "Monte Carlo \\(nsim = 1000\\)")
  expect_output(print(fisher_scoring(approx = expect_mc(nsim = 1000))),
                "Monte Carlo \\(nsim = 1000\\)")
})

test_that("a string signals an error that names the constructor", {
  d <- gaussian1_distrib()
  th <- list(mu = 0, sigma = 1)
  expect_error(distrib_expected_hessian(d, 0, th, approx = "opg"),
               "no longer accepted.*expect_opg\\(\\)")
  expect_error(distrib_expected_hessian(d, 0, th, approx = "bartlett"),
               "expect_bartlett\\(\\)")
  expect_error(distrib_deriv3(d, 0, th, expected = TRUE, approx = "integrate"),
               "expect_integrate\\(\\)")
  expect_error(distrib_deriv4(d, 0, th, expected = TRUE, approx = "mc"),
               "expect_mc\\(nsim = \\.\\.\\.\\)")
  expect_error(fisher_scoring(approx = "mc"), "expect_mc")
  expect_error(fisher_scoring(approx = "nonsense"), "expectation object")
  expect_error(fisher_scoring(approx = 1), "expectation object")
  # the check runs on a family with a closed form too, where the argument
  # would otherwise be dropped without a message
  expect_error(distrib_d2expected_hessian(d, 0, th, approx = "opg"),
               "expect_opg")
  fit <- fit_distrib(d, c(-1, 0.2, 0.7, 1.4, -0.3, 0.8, 2))
  expect_error(vcov(fit, information = "expected", approx = "bartlett"),
               "expect_bartlett")
})

test_that("nsim is no longer an argument of the consumers", {
  d <- gaussian1_distrib()
  th <- list(mu = 0, sigma = 1)
  expect_error(distrib_expected_hessian(d, 0, th, nsim = 10), "expect_mc")
  expect_error(distrib_deriv3(d, 0, th, nsim = 10), "expect_mc")
  expect_error(fisher_scoring(nsim = 10), "unused argument")
  expect_false("nsim" %in% names(formals(fisher_scoring)))
  for (g in list(distrib_expected_hessian, distrib_deriv3, distrib_deriv4)) {
    expect_false("nsim" %in% names(formals(g)))
  }
})

test_that("each object selects its route, and the routes agree with the helpers", {
  d <- dens_only()
  th <- list(mu = 0.3, sigma = 2)
  y <- c(-1, 0.5, 2)
  for (k in 2:4) {
    opg <- distributions7:::expected_derivative(d, y, th, order = k,
                                                approx = expect_opg())
    bar <- distributions7:::expected_derivative(d, y, th, order = k,
                                                approx = expect_bartlett())
    int <- distributions7:::expected_derivative(d, y, th, order = k,
                                                approx = expect_integrate())
    if (k == 2L) {
      expect_identical(opg, distributions7:::expected_by_opg(d, y, th))
    } else {
      # the outer product of the scores has no counterpart above order 2
      expect_identical(opg, bar)
    }
    expect_identical(bar, distributions7:::expected_by_bartlett(d, y, th, k))
    expect_identical(int, distributions7:::expected_by_integrate(d, y, th, k))
  }
  set.seed(4)
  mc <- distributions7:::expected_derivative(d, y, th, order = 2L,
                                             approx = expect_mc(nsim = 300))
  set.seed(4)
  ref <- distributions7:::expected_by_mc(d, y, th, 2L, 300)
  expect_identical(mc, ref)
})

test_that("the generics keep their own defaults", {
  d <- dens_only()
  th <- list(mu = 0.3, sigma = 2)
  y <- c(-1, 0.5, 2)
  expect_identical(distrib_expected_hessian(d, y, th),
                   distrib_expected_hessian(d, y, th, approx = expect_opg()))
  expect_identical(distrib_deriv3(d, y, th, expected = TRUE),
                   distrib_deriv3(d, y, th, expected = TRUE,
                                  approx = expect_integrate()))
  expect_identical(distrib_deriv4(d, y, th, expected = TRUE),
                   distrib_deriv4(d, y, th, expected = TRUE,
                                  approx = expect_integrate()))
  expect_s3_class(eval(formals(distrib_expected_hessian)$approx),
                  "distributions7::ExpectOPG")
  expect_s3_class(eval(formals(distrib_deriv3)$approx),
                  "distributions7::ExpectIntegrate")
  expect_s3_class(eval(formals(distrib_deriv4)$approx),
                  "distributions7::ExpectIntegrate")
  expect_s3_class(eval(formals(fisher_scoring)$approx),
                  "distributions7::ExpectOPG")
})

test_that("the sample size reaches the Monte Carlo draw", {
  d <- dens_only()
  th <- list(mu = 0.3, sigma = 2)
  set.seed(1)
  a <- distrib_expected_hessian(d, 0, th, approx = expect_mc(nsim = 200))
  set.seed(1)
  b <- distrib_expected_hessian(d, 0, th, approx = expect_mc(nsim = 200))
  set.seed(1)
  cc <- distrib_expected_hessian(d, 0, th, approx = expect_mc(nsim = 2000))
  expect_identical(a, b)
  expect_false(identical(a, cc))
})

test_that("a multivariate family reads the object", {
  d <- mvgaussian1_distrib(2)
  th <- distributions7:::align_theta(
    d, list(mu1 = 0, mu2 = 0, sigma_log_L1 = 0, sigma_log_L2 = 0, sigma_L2.1 = 0.5))
  base <- S7::method(distrib_expected_hessian, multivariate_distrib)
  y <- matrix(0, 3, 2)
  expect_error(base(d, y, th, approx = expect_integrate()), "expect_integrate")
  expect_error(base(d, y, th, approx = "mc"), "expect_mc")
  set.seed(2)
  o <- base(d, y, th, approx = expect_opg())
  set.seed(2)
  b <- base(d, y, th, approx = expect_bartlett())
  expect_identical(o, b)
  set.seed(3)
  m1 <- base(d, y, th, approx = expect_mc(nsim = 500))
  set.seed(3)
  m2 <- base(d, y, th, approx = expect_mc(nsim = 500))
  expect_identical(m1, m2)
})

test_that("the multivariate third derivative takes its draw count from the object", {
  d <- mvstudent_t1_distrib(2)
  th <- list(mu1 = 0, mu2 = 0, sigma_log_L1 = 0, sigma_log_L2 = 0,
             sigma_L2.1 = 0.3, nu = 5)
  set.seed(2)
  y <- distrib_rng(d, 4, th)
  set.seed(9)
  a <- distrib_deriv3(d, y, th, expected = TRUE, approx = expect_mc(nsim = 300))
  set.seed(9)
  b <- distrib_deriv3(d, y, th, expected = TRUE, approx = expect_mc(nsim = 300))
  set.seed(9)
  cc <- distrib_deriv3(d, y, th, expected = TRUE, approx = expect_mc(nsim = 3000))
  expect_identical(a, b)
  expect_false(identical(a, cc))
})

test_that("FisherScoring holds the object in one property", {
  fs <- fisher_scoring(approx = expect_mc(nsim = 123))
  expect_true(S7::S7_inherits(fs@approx, ExpectMC))
  expect_identical(fs@approx@nsim, 123)
  expect_false("nsim" %in% names(S7::S7_class(fs)@properties))
  expect_true(S7::S7_inherits(fisher_scoring()@approx, ExpectOPG))
  # the raw class rejects anything that is not an approximation
  expect_error(FisherScoring(approx = "opg", criterion = NULL, maxit = NULL))
})

test_that("fit_distrib() accepts objects and rejects the former strings", {
  set.seed(31)
  d <- gaussian1_distrib()
  y <- distrib_rng(d, 100, list(mu = 1, sigma = 2))
  for (s in c("fisher", "newton", "bfgs")) {
    expect_error(fit_distrib(d, y, method = s),
                 sprintf("string \"%s\" is no longer accepted", s), fixed = FALSE)
  }
  expect_error(fit_distrib(d, y, method = "fisher"), "fisher_scoring\\(")
  expect_error(fit_distrib(d, y, method = "newton"),
               "optimizers7::newton\\(criterion = optimizers7::crit_grad\\(\\)\\)")
  expect_error(fit_distrib(d, y, method = "bfgs"),
               "optimizers7::bfgs\\(criterion = optimizers7::crit_grad\\(\\)\\)")
  expect_error(fit_distrib(d, y, method = "lbfgs"), "optimizers7 optimizer object")
  expect_error(fit_distrib(d, y, method = 3), "optimizers7 optimizer object")
  f1 <- fit_distrib(d, y, method = fisher_scoring())
  f2 <- fit_distrib(d, y)
  expect_identical(f1@coefficients, f2@coefficients)
  expect_identical(f1@method, "Fisher scoring")
  fn <- fit_distrib(d, y, method = optimizers7::newton(
    criterion = optimizers7::crit_grad()))
  expect_equal(unname(fn@coefficients), unname(f1@coefficients), tolerance = 1e-6)
})
