# What a multivariate fit reports. The free values of a structure are
# coordinates, not quantities anybody reads, so the fit carries its variance
# matrix onto the standard deviations and correlations by the delta method.
# Every closed-form Jacobian here is checked against Richardson extrapolation,
# which shares no code with it, and every derived value against arithmetic
# written out by hand.

test_that("the parameter names say which matrix the structure describes", {
  ds <- mvgaussian_distrib(2)
  do <- mvgaussian_distrib(2, omega = parameters7::log_cholesky(2))

  expect_identical(
    ds@params, c("mu1", "mu2", "sigma_log_L1", "sigma_log_L2", "sigma_L2.1")
  )
  expect_identical(
    do@params, c("mu1", "mu2", "omega_log_L1", "omega_log_L2", "omega_L2.1")
  )
  # The two are different models, so identical names would be ambiguous about
  # the one thing the reader has to know.
  expect_false(identical(ds@params, do@params))

  # the scale matrix of a t is written Sigma, so it takes the same prefix
  expect_identical(
    mvstudent_t_distrib(2)@params,
    c("mu1", "mu2", "sigma_log_L1", "sigma_log_L2", "sigma_L2.1", "nu")
  )

  # a diagonal structure carries its own labels through the same prefix
  expect_identical(
    mvgaussian_distrib(3, sigma = parameters7::diagonal_matrix(3))@params,
    c("mu1", "mu2", "mu3", "sigma_log_d1", "sigma_log_d2", "sigma_log_d3")
  )
})

test_that("the derived quantities are the decomposition, written out by hand", {
  d <- mvgaussian_distrib(2)
  th <- list(mu1 = 0.5, mu2 = -0.3, sigma_log_L1 = 0.1, sigma_log_L2 = -0.2,
             sigma_L2.1 = 0.4)
  s <- mv_sigma(d, th)
  der <- mv_derived(d, th)

  expect_identical(names(der$value), c("sd_v1", "sd_v2", "cor_v1_v2"))
  expect_equal(unname(der$value[1:2]), unname(sqrt(diag(s))))
  expect_equal(
    unname(der$value[3]), unname(s[1, 2] / sqrt(s[1, 1] * s[2, 2]))
  )
  # the scales the intervals are built on
  expect_identical(unname(der$transform), c("log", "log", "atanh"))
})

test_that("the derived Jacobian is closed form and agrees with an independent one", {
  for (d in list(mvgaussian_distrib(2),
                 mvgaussian_distrib(3),
                 mvgaussian_distrib(3, sigma = parameters7::diagonal_matrix(3)),
                 mvgaussian_distrib(3, omega = parameters7::log_cholesky(3)),
                 mvstudent_t_distrib(2))) {
    set.seed(61)
    th <- generate_random_theta(d)
    v0 <- unlist(th)
    gfun <- function(v) {
      mv_derived(d, as.list(stats::setNames(v, d@params)))$value
    }
    j <- numDeriv::jacobian(gfun, v0)
    expect_equal(unname(mv_derived(d, th)$jacobian), j,
      tolerance = 1e-6, label = d@distrib_name
    )
  }
})

test_that("a precision parametrisation reports what it describes directly", {
  ds <- mvgaussian_distrib(3)
  do <- mvgaussian_distrib(3, omega = parameters7::log_cholesky(3))

  set.seed(62)
  ths <- generate_random_theta(ds)
  sigma <- mv_sigma(ds, ths)
  eta <- parameters7::param_free(do@param, solve(unname(sigma)))
  tho <- as.list(stats::setNames(
    c(unlist(ths)[1:3], unname(eta)), do@params
  ))

  ds_der <- mv_derived(ds, ths)
  do_der <- mv_derived(do, tho)

  # the standard deviations and correlations of the response are the same
  # whichever side is parametrised: they are properties of the law
  expect_equal(do_der$value[names(ds_der$value)], ds_der$value)

  # and the precision form adds the readings that are its own
  omega <- solve(unname(sigma))
  expect_equal(
    unname(do_der$value[paste0("cvar_v", 1:3)]), 1 / diag(omega)
  )
  expect_equal(
    unname(do_der$value[["pcor_v1_v2"]]),
    -omega[1, 2] / sqrt(omega[1, 1] * omega[2, 2])
  )

  # in two dimensions there is nothing to condition on, so the partial
  # correlation would repeat the correlation and is not printed
  d2 <- mvgaussian_distrib(2, omega = parameters7::log_cholesky(2))
  set.seed(63)
  v2 <- mv_derived(d2, generate_random_theta(d2))$value
  expect_false(any(grepl("^pcor_", names(v2))))
  expect_true(any(grepl("^cvar_", names(v2))))
})

test_that("the base class reports the matrix on its own scale", {
  # A family that says nothing more specific still gets its matrix back, named
  # after the coordinates, rather than the structure's coordinates.
  Odd <- S7::new_class("OddMv2", parent = multivariate_distrib, package = NULL)
  g <- mvgaussian_distrib(2)
  gen <- mv_sigma
  S7::method(gen, Odd) <- function(distrib, theta) {
    matrix(c(2, 0.5, 0.5, 3), 2, 2)
  }
  d <- Odd(
    distrib_name = "odd", dimension = "multivariate", n_dim = 2L,
    bounds = c(-Inf, Inf), params = c("a", "b"),
    params_interpretation = c(a = "location", b = "location"), n_params = 2L,
    params_bounds = list(a = c(-Inf, Inf), b = c(-Inf, Inf)),
    link_params = list(
      a = linkfunctions7::identity_link(), b = linkfunctions7::identity_link()
    )
  )
  der <- mv_derived(d, list(a = 0, b = 0))
  expect_identical(
    names(der$value), c("sigma_v1_v1", "sigma_v2_v1", "sigma_v2_v2")
  )
  expect_equal(unname(der$value), c(2, 0.5, 3))
})

test_that("mv_summary carries the standard errors across and respects the domains", {
  d <- mvgaussian_distrib(2)
  true <- list(mu1 = 0, mu2 = 1, sigma_log_L1 = 0, sigma_log_L2 = 0,
               sigma_L2.1 = 0.7)
  set.seed(64)
  y <- distrib_rng(d, 1000, true)
  fit <- fit_distrib(d, y)
  tab <- mv_summary(fit)

  expect_identical(rownames(tab), c("sd_v1", "sd_v2", "cor_v1_v2"))
  expect_identical(names(tab)[1:2], c("Estimate", "Std. Error"))

  # the estimates are the sample ones, since the gaussian's MLE is
  s_hat <- crossprod(sweep(y, 2L, colMeans(y))) / nrow(y)
  expect_equal(tab["sd_v1", "Estimate"], sqrt(s_hat[1, 1]), tolerance = 1e-4)
  expect_equal(tab["cor_v1_v2", "Estimate"],
    s_hat[1, 2] / sqrt(s_hat[1, 1] * s_hat[2, 2]), tolerance = 1e-4)

  # the standard error against the delta method done here from scratch
  der <- mv_derived(d, as.list(coef(fit)))
  se <- sqrt(diag(der$jacobian %*% vcov(fit) %*% t(der$jacobian)))
  expect_equal(tab[["Std. Error"]], unname(se))

  # an interval built on the raw scale would routinely leave the domain; these
  # are built on log and on Fisher's z and mapped back
  expect_true(all(tab[1:2, 3] > 0))
  expect_true(all(abs(unlist(tab[3, 3:4])) < 1))
  expect_true(all(tab[, 3] < tab[, 1] & tab[, 1] < tab[, 4]))

  # a standard deviation's interval is not symmetric, which is the point of
  # building it on the log scale
  expect_false(isTRUE(all.equal(
    tab[1, 1] - tab[1, 3], tab[1, 4] - tab[1, 1], tolerance = 1e-6
  )))

  # the level is honoured
  wide <- mv_summary(fit, level = 0.99)
  expect_true(all(wide[, 3] < tab[, 3]))
  expect_true(all(wide[, 4] > tab[, 4]))
  expect_equal(wide[["Estimate"]], tab[["Estimate"]])
})

test_that("the correlation's standard error matches the known asymptotic one", {
  # For a bivariate gaussian the maximum likelihood correlation has asymptotic
  # variance (1 - rho^2)^2 / n. That is a statement about the model, computed
  # nowhere in the package, so it tests the whole delta-method chain at once.
  d <- mvgaussian_distrib(2)
  rho <- 0.6
  s <- matrix(c(1, rho, rho, 1), 2, 2)
  eta <- parameters7::param_free(d@param, s)
  true <- as.list(stats::setNames(c(0, 0, unname(eta)), d@params))

  n <- 4000
  set.seed(65)
  y <- distrib_rng(d, n, true)
  tab <- mv_summary(fit_distrib(d, y))
  expect_equal(
    tab["cor_v1_v2", "Std. Error"],
    sqrt((1 - tab["cor_v1_v2", "Estimate"]^2)^2 / n),
    tolerance = 0.05
  )

  # and a standard deviation's, which is sigma / sqrt(2n)
  expect_equal(
    tab["sd_v1", "Std. Error"], tab["sd_v1", "Estimate"] / sqrt(2 * n),
    tolerance = 0.05
  )
})

test_that("mv_summary refuses what it is not for", {
  set.seed(66)
  uni <- fit_distrib(gaussian1_distrib(), stats::rnorm(100))
  expect_error(mv_summary(uni), "multivariate fit")
  expect_error(mv_summary("not a fit"), "distrib_fit")

  d <- mvgaussian_distrib(2)
  set.seed(67)
  fit <- fit_distrib(d, distrib_rng(d, 200, generate_random_theta(d)))
  expect_error(mv_summary(fit, level = 1), "in \\(0, 1\\)")
})

test_that("print shows the interpretable blocks rather than only the coordinates", {
  d <- mvgaussian_distrib(2)
  set.seed(68)
  y <- distrib_rng(d, 400, list(mu1 = 0, mu2 = 0, sigma_log_L1 = 0,
                                sigma_log_L2 = 0, sigma_L2.1 = 0.3))
  out <- utils::capture.output(print(fit_distrib(d, y)))

  expect_true(any(grepl("^Standard deviations:", out)))
  expect_true(any(grepl("^Correlations:", out)))
  expect_true(any(grepl("^sd_v1", out)))
  expect_true(any(grepl("^cor_v1_v2", out)))
  expect_true(any(grepl("^Location:", out)))

  # the t names its diagonal quantities for what they are: the scale matrix is
  # not the covariance
  dt <- mvstudent_t_distrib(2)
  set.seed(69)
  yt <- distrib_rng(dt, 400, list(mu1 = 0, mu2 = 0, sigma_log_L1 = 0,
                                  sigma_log_L2 = 0, sigma_L2.1 = 0.3, nu = 6))
  out_t <- utils::capture.output(
    print(fit_distrib(dt, yt, method = optimizers7::newton()))
  )
  expect_true(any(grepl("^Scale standard deviations:", out_t)))
  expect_false(any(grepl("^Standard deviations:", out_t)))
})


test_that("the expected information of a multivariate fit is the right size", {
  # The link-scale branch of distrib_expected_hessian() built its first-order
  # term with rep(0, length(y)), and for a matrix response length(y) counts
  # ENTRIES. Recycled against the components it inflated every diagonal entry
  # of the information by a factor of p -- the first-order term of the order-2
  # chain rule appears only on the diagonal -- so every standard error of a
  # multivariate fit came out sqrt(p) times too small. The observed Hessian was
  # correct throughout, which is why a fit by Newton disagreed with a fit by
  # Fisher scoring on the same data.
  set.seed(81)
  for (p in 2:4) {
    d <- mvgaussian_distrib(p)
    th <- generate_random_theta(d)
    n <- 500
    y <- distrib_rng(d, n, th)
    s <- unname(mv_sigma(d, th))

    # the mean block of the expected information is n * Sigma^{-1}, exactly
    info <- -distributions7:::fit_hess_matrix(d, y, th, expected = TRUE)
    expect_equal(unname(info[seq_len(p), seq_len(p)]), n * solve(s),
      tolerance = 1e-8, label = paste("p =", p)
    )

    # and the expected and observed informations agree in the mean block,
    # which for a gaussian they do exactly
    obs <- -distributions7:::fit_hess_matrix(d, y, th, expected = FALSE)
    expect_equal(unname(obs[seq_len(p), seq_len(p)]),
      unname(info[seq_len(p), seq_len(p)]),
      tolerance = 1e-8, label = paste("p =", p)
    )
  }
})

test_that("the reported standard errors match the asymptotic ones", {
  # Two facts about the multivariate gaussian MLE, neither computed anywhere in
  # the package: se(mu_j) = sd_j / sqrt(n), and Var(S_jj) = 2 Sigma_jj^2 / n,
  # so se(sd_j) = sd_j / sqrt(2n). They hold whatever p and whatever the
  # correlations, which makes them a check on the whole chain from the
  # information to the delta method.
  set.seed(82)
  for (p in c(2L, 4L)) {
    d <- mvgaussian_distrib(p)
    n <- 2000
    y <- distrib_rng(d, n, generate_random_theta(d))
    fit <- fit_distrib(d, y)
    tab <- mv_summary(fit)

    sd1 <- tab["sd_v1", "Estimate"]
    expect_equal(tab["sd_v1", "Std. Error"], sd1 / sqrt(2 * n),
      tolerance = 1e-6, label = paste("sd, p =", p)
    )
    expect_equal(unname(fit@se[["mu1"]]), sd1 / sqrt(n),
      tolerance = 1e-6, label = paste("mu, p =", p)
    )
  }
})

test_that("the two parametrisations report the same uncertainty", {
  # They are the same model, so a derived quantity has the same standard error
  # under either. A disagreement means one of the two informations is wrong.
  set.seed(83)
  ds <- mvgaussian_distrib(3)
  th <- generate_random_theta(ds)
  y <- distrib_rng(ds, 1000, th)
  do <- mvgaussian_distrib(3, omega = parameters7::log_cholesky(3))

  a <- mv_summary(fit_distrib(ds, y))
  b <- mv_summary(fit_distrib(do, y))
  common <- rownames(a)
  expect_equal(a[common, "Estimate"], b[common, "Estimate"], tolerance = 1e-5)
  expect_equal(a[common, "Std. Error"], b[common, "Std. Error"],
    tolerance = 1e-4
  )
})


test_that("a structured matrix reports the quantities the family is about", {
  set.seed(77)
  d <- mvgaussian_distrib(6, sigma = parameters7::autoregressive(6, order = 2))
  th <- as.list(stats::setNames(
    c(rep(0, 6), log(3), atanh(0.7), atanh(-0.3)), d@params
  ))
  y <- distrib_rng(d, 400, th)
  fit <- fit_distrib(d, y)
  s <- mv_summary(fit)
  blk <- attr(s, "block")

  expect_true("Autoregressive structure" %in% blk)
  rows <- rownames(s)[blk == "Autoregressive structure"]
  expect_identical(rows, c("scale", "pacf1", "pacf2", "phi1", "phi2"))

  # The block is not a restatement of the covariance: the lag-one correlation
  # IS the first partial autocorrelation, and the two must agree to the digit,
  # while the coefficients appear nowhere else in the summary.
  expect_equal(s["pacf1", "Estimate"], s["cor_v1_v2", "Estimate"],
    tolerance = 1e-10)
  expect_equal(s["pacf1", "Std. Error"], s["cor_v1_v2", "Std. Error"],
    tolerance = 1e-10)
  expect_false(any(grepl("^phi", rownames(s)[blk != "Autoregressive structure"])))

  # an order-q autoregression has phi_q = r_q whatever the data
  expect_equal(s["phi2", "Estimate"], s["pacf2", "Estimate"], tolerance = 1e-10)

  # every interval stays in the set its quantity lives in, which is what
  # building it on the declared scale is for
  expect_gt(s["scale", "2.5%"], 0)
  expect_gt(s["pacf1", "2.5%"], -1)
  expect_lt(s["pacf1", "97.5%"], 1)
  expect_lt(s["pacf2", "97.5%"], 1)
})


test_that("the delta method behind the block agrees with a numerical Jacobian", {
  # The Jacobian comes from the Levinson-Durbin recursion in jets; the
  # reference differentiates the map from the free vector to the quantities,
  # which shares nothing with it.
  skip_if_not_installed("numDeriv")
  d <- mvgaussian_distrib(5, sigma = parameters7::autoregressive(5, order = 2))
  th <- as.list(stats::setNames(
    c(0.2, -0.1, 0.4, 0, 0.3, log(2), atanh(0.6), atanh(-0.25)), d@params
  ))
  der <- mv_derived(d, th)
  keep <- attr(der$block, "names")[der$block == "Autoregressive structure"]
  v <- unlist(th, use.names = FALSE)
  num <- numDeriv::jacobian(function(x) {
    mv_derived(d, as.list(stats::setNames(x, d@params)))$value[keep]
  }, v)
  expect_equal(unname(der$jacobian[keep, ]), num, tolerance = 1e-7)
})


test_that("a precision parametrisation says which matrix the block describes", {
  d <- mvgaussian_distrib(4, omega = parameters7::ar1(4))
  th <- as.list(stats::setNames(c(0, 0, 0, 0, log(2), atanh(0.5)), d@params))
  der <- mv_derived(d, th)
  expect_true("Autoregressive structure (precision)" %in% der$block)
  # and the quantities are those of Omega, not of its inverse
  expect_equal(unname(der$value[["rho"]]), 0.5, tolerance = 1e-12)
})


test_that("an unstructured matrix adds no block and is unchanged", {
  d <- mvgaussian_distrib(3)
  th <- as.list(stats::setNames(c(0, 0, 0, 0.1, -0.2, 0.05, 0.4, 0.3, -0.1),
                                d@params))
  der <- mv_derived(d, th)
  expect_setequal(unique(der$block),
                  c("Standard deviations", "Correlations"))
})
