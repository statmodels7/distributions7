# The method of an expectation as an S7 object: construction, the validation
# of the argument, the routing to the computations, the meaning of NULL and
# the use of an explicit method by a family that has a closed form.

dens_only <- function() {
  Dens <- S7::new_class("DensExpectMethod", parent = continuous_distrib)
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

kinked_only <- function() {
  Dens <- S7::new_class("DensKinkExpectMethod", parent = continuous_distrib)
  S7::method(distrib_pdf, Dens) <- function(distrib, y, theta, log = FALSE, ...) {
    lp <- -log(2 * theta[[2]]) - abs(y - theta[[1]]) / theta[[2]]
    if (log) lp else exp(lp)
  }
  Dens(distrib_name = "densk", dimension = "univariate",
       bounds = c(-Inf, Inf), params = c("mu", "sigma"),
       params_interpretation = c(mu = "location", sigma = "scale"),
       n_params = 2L,
       params_bounds = list(mu = c(-Inf, Inf), sigma = c(0, Inf)),
       link_params = list(mu = linkfunctions7::identity_link(),
                          sigma = linkfunctions7::log_link()),
       params_smooth = c(mu = FALSE, sigma = TRUE))
}

test_that("the constructors return objects of the classes of the hierarchy", {
  for (a in list(expect_opg(), expect_integrate(), expect_mc(),
                 expect_integrate(identity = "score"),
                 expect_mc(identity = "score"))) {
    expect_true(S7::S7_inherits(a, ExpectMethod))
  }
  expect_true(S7::S7_inherits(expect_opg(), ExpectOPG))
  expect_true(S7::S7_inherits(expect_integrate(), ExpectIntegrate))
  expect_true(S7::S7_inherits(expect_mc(), ExpectMC))
  expect_false(S7::S7_inherits(expect_opg(), ExpectMC))
  expect_false(S7::S7_inherits(expect_integrate(), ExpectOPG))
  # the parent is abstract
  expect_error(ExpectMethod())
})

test_that("the identity is a property of integrate and mc, defaulting to hessian", {
  expect_identical(expect_integrate()@identity, "hessian")
  expect_identical(expect_integrate(identity = "score")@identity, "score")
  expect_identical(expect_mc()@identity, "hessian")
  expect_identical(expect_mc(identity = "score")@identity, "score")
  expect_identical(expect_mc(nsim = 50, identity = "score")@nsim, 50)
  expect_error(expect_integrate(identity = "bartlett"), "should be one of")
  expect_error(expect_mc(identity = "opg"), "should be one of")
  expect_error(expect_integrate(identity = c("hessian", "score", "x")))
  expect_false("identity" %in% names(S7::S7_class(expect_opg())@properties))
})

test_that("expect_mc() carries nsim and validates it", {
  expect_identical(expect_mc()@nsim, 10000)
  expect_identical(expect_mc(nsim = 250)@nsim, 250)
  expect_identical(expect_mc(nsim = 20L)@nsim, 20L)
  for (bad in list(0, -1, NA_real_, Inf, c(10, 20), "10", NULL, 0.5)) {
    expect_error(expect_mc(nsim = bad), "single positive number")
  }
})

test_that("the print methods name the method and the identity", {
  expect_output(print(expect_opg()), "outer product")
  expect_output(print(expect_integrate()), "integration")
  expect_output(print(expect_integrate(identity = "score")), "identity = score")
  expect_output(print(expect_mc(nsim = 1000)), "Monte Carlo \\(nsim = 1000\\)")
  expect_output(print(expect_mc(nsim = 1000, identity = "score")),
                "Monte Carlo \\(nsim = 1000, identity = score\\)")
  expect_output(print(fisher_scoring(method = expect_mc(nsim = 1000))),
                "Monte Carlo \\(nsim = 1000\\)")
  expect_output(print(fisher_scoring()), "closed form where the family has one")
})

test_that("a value that is not NULL or an ExpectMethod object is rejected", {
  d <- gaussian1_distrib()
  th <- list(mu = 0, sigma = 1)
  msg <- "must be NULL or an ExpectMethod object such as expect_opg\\(\\)"
  expect_error(distrib_expected_hessian(d, 0, th, method = "opg"), msg)
  expect_error(distrib_deriv3(d, 0, th, expected = TRUE, method = "integrate"), msg)
  expect_error(distrib_deriv4(d, 0, th, expected = TRUE, method = "mc"), msg)
  expect_error(fisher_scoring(method = "nonsense"), msg)
  expect_error(fisher_scoring(method = 1), msg)
  expect_error(distrib_d2expected_hessian(d, 0, th, method = "opg"), msg)
  expect_error(expected_info(method = "opg"), msg)
})

test_that("the sample size is an argument of expect_mc() only", {
  expect_error(fisher_scoring(nsim = 10), "unused argument")
  expect_false("nsim" %in% names(formals(fisher_scoring)))
  for (g in list(distrib_expected_hessian, distrib_deriv3, distrib_deriv4)) {
    expect_false("nsim" %in% names(formals(g)))
  }
})

test_that("expect_integrate() carries the tolerances of the quadrature and the series", {
  a <- expect_integrate()
  qf <- formals(numericals7::quad_vec)
  sf <- formals(numericals7::series_vec)
  expect_identical(a@quad_atol, qf$atol)
  expect_identical(a@quad_rtol, qf$rtol)
  expect_identical(a@max_panels, qf$max_panels)
  expect_identical(a@series_atol, sf$atol)
  expect_identical(a@series_rtol, sf$rtol)
  expect_identical(a@max_terms, sf$max_terms)
  for (nm in c("quad_atol", "quad_rtol", "max_panels", "series_atol",
               "series_rtol", "max_terms")) {
    for (bad in list(0, -1, NA_real_, Inf, c(1, 2), "1", NULL)) {
      expect_error(do.call(expect_integrate, stats::setNames(list(bad), nm)),
                   "single positive finite number")
    }
  }
  expect_output(print(expect_integrate(quad_atol = 1e-12)), "quad_atol = 0.000000000001")
  expect_output(print(expect_integrate()), "density$")
})

test_that("the tolerances of expect_integrate() reach the quadrature and the series", {
  d <- dens_only()
  th <- list(mu = 0.3, sigma = 2)
  y <- c(-1, 0.5)
  def <- distrib_deriv3(d, y, th, expected = TRUE, method = expect_integrate())
  same <- distrib_deriv3(d, y, th, expected = TRUE,
                         method = expect_integrate(quad_atol = 1e-10, quad_rtol = 1e-8))
  expect_identical(def, same)
  loose <- distrib_deriv3(d, y, th, expected = TRUE,
                          method = expect_integrate(quad_atol = 1e-2, quad_rtol = 1e-2))
  expect_false(identical(def, loose))
  # the tolerances apply to the score identity too
  s1 <- distrib_deriv3(d, y, th, expected = TRUE,
                       method = expect_integrate(identity = "score"))
  s2 <- distrib_deriv3(d, y, th, expected = TRUE,
                       method = expect_integrate(identity = "score",
                                                 quad_atol = 1e-2, quad_rtol = 1e-2))
  expect_false(identical(s1, s2))
  # a series that cannot meet its tolerance in max_terms terms is an error
  pm <- function(a) {
    old <- distributions7:::set_integration_settings(a)
    on.exit(distributions7:::set_integration_settings(old))
    expectation(poisson_distrib(), function(y, theta) y, list(mu = 500))
  }
  expect_equal(pm(expect_integrate()), 500, tolerance = 1e-8)
  expect_error(pm(expect_integrate(max_terms = 20)), "did not converge")
  # the settings are restored after the call
  distrib_deriv3(d, y, th, expected = TRUE,
                 method = expect_integrate(quad_atol = 1e-2, quad_rtol = 1e-2))
  expect_identical(distrib_deriv3(d, y, th, expected = TRUE), def)
})

test_that("each object selects its route, and the routes agree with the helpers", {
  d <- dens_only()
  th <- list(mu = 0.3, sigma = 2)
  y <- c(-1, 0.5, 2)
  for (k in 2:4) {
    opg <- distributions7:::expected_derivative(d, y, th, order = k,
                                                method = expect_opg())
    sco <- distributions7:::expected_derivative(
      d, y, th, order = k, method = expect_integrate(identity = "score"))
    hes <- distributions7:::expected_derivative(d, y, th, order = k,
                                                method = expect_integrate())
    expect_identical(opg, distributions7:::expected_by_opg(d, y, th, k))
    expect_identical(sco, distributions7:::expected_by_bartlett(d, y, th, k))
    expect_identical(hes, distributions7:::expected_by_integrate(d, y, th, k))
    # the empirical sums are not the quadrature, at any order
    expect_false(isTRUE(all.equal(opg, sco)))
    set.seed(4)
    mc <- distributions7:::expected_derivative(d, y, th, order = k,
                                               method = expect_mc(nsim = 300))
    set.seed(4)
    expect_identical(mc, distributions7:::expected_by_mc(d, y, th, k, 300))
    set.seed(4)
    mcs <- distributions7:::expected_derivative(
      d, y, th, order = k, method = expect_mc(nsim = 300, identity = "score"))
    set.seed(4)
    expect_identical(mcs, distributions7:::expected_by_mc(d, y, th, k, 300,
                                                          identity = "score"))
  }
})

test_that("NULL resolves to opg at order 2 and to integration above it", {
  d <- dens_only()
  th <- list(mu = 0.3, sigma = 2)
  y <- c(-1, 0.5, 2)
  expect_identical(distributions7:::resolve_method(NULL, 2L), expect_opg())
  expect_identical(distributions7:::resolve_method(NULL, 3L), expect_integrate())
  expect_identical(distributions7:::resolve_method(NULL, 4L), expect_integrate())
  expect_identical(distributions7:::resolve_method(expect_mc(), 2L), expect_mc())
  expect_identical(distrib_expected_hessian(d, y, th),
                   distrib_expected_hessian(d, y, th, method = expect_opg()))
  expect_identical(distrib_deriv3(d, y, th, expected = TRUE),
                   distrib_deriv3(d, y, th, expected = TRUE,
                                  method = expect_integrate()))
  expect_identical(distrib_deriv4(d, y, th, expected = TRUE),
                   distrib_deriv4(d, y, th, expected = TRUE,
                                  method = expect_integrate()))
  expect_identical(distributions7:::expected_derivative(d, y, th, order = 3L),
                   distributions7:::expected_derivative(
                     d, y, th, order = 3L, method = expect_integrate()))
  for (g in list(distrib_expected_hessian, distrib_deriv3, distrib_deriv4,
                 fisher_scoring, expected_info)) {
    expect_null(formals(g)$method)
  }
})

test_that("expect_opg() above order 2 is the empirical sum of products, not a quadrature", {
  g <- gaussian1_distrib()
  d <- dens_only()
  th <- list(mu = 0.3, sigma = 2)
  y <- c(-1, 0.5, 2, 3.1)
  # order 3: -(l_ij l_k + l_ik l_j + l_jk l_i + l_i l_j l_k), by hand from the
  # observed derivatives of the Gaussian
  l1 <- distrib_gradient(g, y, th)
  l2 <- distrib_hessian(g, y, th)
  o3 <- distrib_deriv3(d, y, th, expected = TRUE, method = expect_opg())
  by_hand <- -(l2$mu_mu * l1$sigma + 2 * l2$mu_sigma * l1$mu + l1$mu^2 * l1$sigma)
  expect_equal(o3$mu_mu_sigma, by_hand, tolerance = 1e-5)
  # it varies with the observation, which an expectation does not
  expect_gt(length(unique(round(o3$mu_mu_sigma, 8))), 1L)
  # order 4: the partitions of four indices into at least two blocks
  l3 <- distrib_deriv3(g, y, th)
  o4 <- distrib_deriv4(d, y, th, expected = TRUE, method = expect_opg())
  # the partition sum evaluated independently from the observed derivatives
  parts <- Filter(function(pp) length(pp) > 1L, numericals7::set_partitions(4L))
  idx <- c(1L, 1L, 2L, 2L)              # mu mu sigma sigma
  obs <- list(l1, l2, l3)
  val <- 0
  for (pp in parts) {
    pr <- 1
    for (B in pp) {
      blk <- sort(idx[B])
      key <- paste(g@params[blk], collapse = "_")
      pr <- pr * obs[[length(blk)]][[key]]
    }
    val <- val + pr
  }
  expect_equal(o4$mu_mu_sigma_sigma, -val, tolerance = 1e-5)
  # the sum over many observations estimates the expectation
  set.seed(1)
  yy <- distrib_rng(g, 40000, th)
  est <- mean(distrib_deriv3(d, yy, th, expected = TRUE, method = expect_opg())$sigma_sigma_sigma)
  exact <- distrib_deriv3(g, 0, th, expected = TRUE)$sigma_sigma_sigma
  expect_equal(est, exact, tolerance = 0.05)
  # and no quadrature is involved: a density that cannot be integrated is read
  set.seed(2)
  expect_no_error(distrib_deriv3(d, y, th, expected = TRUE, method = expect_opg()))
})

test_that("an explicit method is used also by a family with a closed form", {
  g <- gaussian1_distrib()
  th <- list(mu = 0.3, sigma = 2)
  y <- c(-1, 0.5, 2)
  exact <- distrib_expected_hessian(g, y, th)
  # opg is the outer product of the scores at the data, not the closed form
  o <- distrib_expected_hessian(g, y, th, method = expect_opg())
  expect_identical(o, distributions7:::expected_by_opg(g, y, th))
  expect_false(isTRUE(all.equal(o, exact)))
  # integration reproduces the closed form to its tolerance
  i <- distrib_expected_hessian(g, y, th, method = expect_integrate())
  expect_equal(i, exact, tolerance = 1e-7)
  expect_false(identical(i, exact))
  s <- distrib_expected_hessian(g, y, th, method = expect_integrate(identity = "score"))
  expect_equal(s, exact, tolerance = 1e-7)
  # simulation carries sampling error and depends on the seed
  set.seed(1)
  m1 <- distrib_expected_hessian(g, y, th, method = expect_mc(nsim = 500))
  set.seed(2)
  m2 <- distrib_expected_hessian(g, y, th, method = expect_mc(nsim = 500))
  expect_false(identical(m1, m2))
  expect_equal(m1$sigma_sigma, exact$sigma_sigma, tolerance = 0.3)
  # the same holds at the higher orders and on the link scale
  d3 <- distrib_deriv3(g, y, th, expected = TRUE)
  expect_equal(distrib_deriv3(g, y, th, expected = TRUE,
                              method = expect_integrate()), d3, tolerance = 1e-7)
  expect_false(identical(
    distrib_deriv3(g, y, th, expected = TRUE, method = expect_opg()), d3))
  d4l <- distrib_deriv4(g, y, th, expected = TRUE, scale = "link")
  expect_equal(distrib_deriv4(g, y, th, expected = TRUE, scale = "link",
                              method = expect_integrate()), d4l, tolerance = 1e-6)
  # with expected = FALSE the method is not read
  expect_identical(distrib_deriv3(g, y, th, method = expect_opg()),
                   distrib_deriv3(g, y, th))
})

test_that("an explicit method is used by a wrapper and by Fisher scoring", {
  g <- gaussian1_distrib()
  set.seed(31)
  y <- distrib_rng(g, 200, list(mu = 1, sigma = 2))
  f0 <- fit_distrib(g, y)
  f1 <- fit_distrib(g, y, method = fisher_scoring(method = expect_integrate()))
  expect_equal(f1@coefficients, f0@coefficients, tolerance = 1e-6)
  set.seed(5)
  f2 <- fit_distrib(g, y, method = fisher_scoring(method = expect_mc(nsim = 200)))
  expect_equal(f2@coefficients, f0@coefficients, tolerance = 1e-3)
  f3 <- fit_distrib(g, y, method = fisher_scoring(method = expect_opg()))
  expect_equal(f3@coefficients, f0@coefficients, tolerance = 1e-6)
})

test_that("the Hessian identity at a kink signals an error that names the score identity", {
  k <- kinked_only()
  th <- list(mu = 0.2, sigma = 1.5)
  y <- c(-1, 0.5, 2)
  expect_error(distrib_expected_hessian(k, y, th, method = expect_integrate()),
               'identity = "score"')
  expect_error(distrib_expected_hessian(k, y, th, method = expect_mc(nsim = 100)),
               'identity = "score"')
  expect_error(distrib_dexpected_hessian(k, y, th, method = expect_integrate()),
               'identity = "score"')
  # the score identity and opg give the information 1 / sigma^2 for the location
  s <- distrib_expected_hessian(k, y, th, method = expect_integrate(identity = "score"))
  expect_equal(s$mu_mu, rep(-1 / th$sigma^2, 3), tolerance = 1e-6)
  o <- distrib_expected_hessian(k, y, th, method = expect_opg())
  expect_equal(o$mu_mu, rep(-1 / th$sigma^2, 3), tolerance = 1e-6)
  set.seed(1)
  m <- distrib_expected_hessian(k, y, th, method = expect_mc(nsim = 500, identity = "score"))
  expect_equal(m$mu_mu, rep(-1 / th$sigma^2, 3), tolerance = 1e-6)
  # the default of a family without a closed form is opg and does not error
  expect_identical(distrib_expected_hessian(k, y, th), o)
  # a smooth family is unaffected, and the check is at order 2 only
  expect_no_error(distrib_deriv3(k, y, th, expected = TRUE, method = expect_integrate()))
})

test_that("the sample size reaches the Monte Carlo draw", {
  d <- dens_only()
  th <- list(mu = 0.3, sigma = 2)
  set.seed(1)
  a <- distrib_expected_hessian(d, 0, th, method = expect_mc(nsim = 200))
  set.seed(1)
  b <- distrib_expected_hessian(d, 0, th, method = expect_mc(nsim = 200))
  set.seed(1)
  cc <- distrib_expected_hessian(d, 0, th, method = expect_mc(nsim = 2000))
  expect_identical(a, b)
  expect_false(identical(a, cc))
})

test_that("a multivariate family reads the object", {
  d <- mvgaussian1_distrib(2)
  th <- distributions7:::align_theta(
    d, list(mu1 = 0, mu2 = 0, sigma_log_L1 = 0, sigma_log_L2 = 0, sigma_L2.1 = 0.5))
  base <- S7::method(distrib_expected_hessian, multivariate_distrib)
  y <- matrix(0, 3, 2)
  expect_error(base(d, y, th, method = expect_integrate()), "expect_integrate")
  expect_error(base(d, y, th, method = "mc"), "must be NULL or an ExpectMethod object")
  set.seed(2)
  o <- base(d, y, th, method = expect_opg())
  set.seed(2)
  b <- base(d, y, th, method = expect_mc(nsim = 10000, identity = "score"))
  expect_identical(o, b)
  set.seed(2)
  expect_identical(base(d, y, th), o)
  set.seed(3)
  m1 <- base(d, y, th, method = expect_mc(nsim = 500))
  set.seed(3)
  m2 <- base(d, y, th, method = expect_mc(nsim = 500))
  expect_identical(m1, m2)
  # an explicit method reaches a multivariate family that has a closed form
  set.seed(4)
  e <- distrib_expected_hessian(d, y, th, method = expect_mc(nsim = 2000))
  ex <- distrib_expected_hessian(d, y, th)
  expect_false(identical(e, ex))
  expect_equal(e$mu1_mu1, ex$mu1_mu1, tolerance = 0.1)
})

test_that("the multivariate third derivative takes its draw count from the object", {
  d <- mvstudent_t1_distrib(2)
  th <- list(mu1 = 0, mu2 = 0, sigma_log_L1 = 0, sigma_log_L2 = 0,
             sigma_L2.1 = 0.3, nu = 5)
  set.seed(2)
  y <- distrib_rng(d, 4, th)
  set.seed(9)
  a <- distrib_deriv3(d, y, th, expected = TRUE, method = expect_mc(nsim = 300))
  set.seed(9)
  b <- distrib_deriv3(d, y, th, expected = TRUE, method = expect_mc(nsim = 300))
  set.seed(9)
  cc <- distrib_deriv3(d, y, th, expected = TRUE, method = expect_mc(nsim = 3000))
  expect_identical(a, b)
  expect_false(identical(a, cc))
})

test_that("FisherScoring holds NULL or the object in one property", {
  fs <- fisher_scoring(method = expect_mc(nsim = 123))
  expect_true(S7::S7_inherits(fs@method, ExpectMC))
  expect_identical(fs@method@nsim, 123)
  expect_false("nsim" %in% names(S7::S7_class(fs)@properties))
  expect_null(fisher_scoring()@method)
  # the raw class rejects anything that is not NULL or a method
  expect_error(FisherScoring(method = "opg", criterion = NULL, maxit = NULL))
})

test_that("fit_distrib() accepts objects and rejects other values of method", {
  set.seed(31)
  d <- gaussian1_distrib()
  y <- distrib_rng(d, 100, list(mu = 1, sigma = 2))
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
